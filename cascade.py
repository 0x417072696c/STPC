# -*- coding: utf-8 -*-
"""
cascade.py – оригинальное ядро + SASA + глубина (расстояние до центра масс).
"""
import numpy as np
import math
import os

J_NAT = 0.125
J_NON = 0.08
T = 300.0
kB = 8.617333262145e-5

def monomer_free_energy(nh, nd, J, N):
    if N <= 0: return 0.0
    return nh * (-J_NAT) + nd * (-J) - T * kB * (np.log(N) if N>0 else 0)

def dimer_free_energy(state1, state2, J_interface, nd_interface):
    F1 = monomer_free_energy(state1['nh'], state1['nd'], state1['J'], state1['N'])
    F2 = monomer_free_energy(state2['nh'], state2['nd'], state2['J'], state2['N'])
    F_int = nd_interface * (-J_interface)
    return F1 + F2 + F_int

def secondary_structure_J(sequence, pos, ss_map=None):
    if ss_map and 0 <= pos < len(ss_map):
        ss = ss_map[pos]
        if ss == 'E': return 0.125, 8
        elif ss == 'H': return 0.115, 12
        else: return 0.08, 10
    pred = predict_ss_chou_fasman(sequence)
    if pred and pos < len(pred):
        ss = pred[pos]
        if ss == 'E': return 0.125, 8
        elif ss == 'H': return 0.115, 12
        else: return 0.08, 10
    return 0.08, 10

def predict_ss_chou_fasman(sequence):
    return ['C'] * len(sequence)

def estimate_aggregation_rate(sequence, pos, wt, alt):
    return 0.0, 0.0

def buried_charge_penalty(aa, exposed):
    if not exposed and aa in ("K","R","D","E"): return 1.5
    return 0.0

def detect_disulfide_bonds(sg_coords, threshold=2.5):
    bonds = set()
    residues = list(sg_coords.items())
    for i in range(len(residues)):
        ri, ci = residues[i]
        for j in range(i+1, len(residues)):
            rj, cj = residues[j]
            if np.linalg.norm(np.array(ci)-np.array(cj)) <= threshold:
                bonds.add((min(ri,rj), max(ri,rj)))
    return bonds

def detect_metal_coordinations(metal_sites, pdb_file, threshold=3.0):
    coordinated = set()
    if not metal_sites: return coordinated
    atoms = []
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    try:
                        res_num = int(line[22:26].strip())
                        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                        atoms.append((res_num, (x,y,z)))
                    except: pass
    except: pass
    for _, mxyz in metal_sites:
        for rn, axyz in atoms:
            if np.linalg.norm(np.array(mxyz)-np.array(axyz)) < threshold:
                coordinated.add(rn)
    return coordinated

def detect_steric_clash(sequence, pos, alt_aa, ss_map, pdb_file):
    vols = {'A':88.6,'V':140.0,'I':166.7,'L':166.7,'F':189.9,'W':227.8,
            'M':162.9,'P':112.7,'G':60.1,'S':89.0,'T':116.1,'C':108.5,
            'N':114.1,'Q':143.8,'H':153.2,'E':138.4,'D':111.1,'K':168.6,'R':173.4,'Y':193.6}
    return alt_aa in vols and vols.get(alt_aa,100) > 150

def surface_charge_inversion(wt_aa, alt_aa, exposed):
    return exposed and wt_aa in {'E','D','R','K'} and alt_aa in {'V','I','L','F','W'}

def detect_dimer_interface(pdb_path):
    try:
        import pdb_parser
        return pdb_parser.extract_dimer_interface_residues(pdb_path)
    except:
        return set()

def calibrate_to_protherm(exp_ddg, calc_hb):
    if len(exp_ddg) < 2 or len(calc_hb) < 2: return None
    A = np.vstack([calc_hb, np.ones(len(calc_hb))]).T
    try:
        scale, intercept = np.linalg.lstsq(A, exp_ddg, rcond=None)[0]
        return scale, intercept
    except: return None

# ---------- SASA + ГЛУБИНА (расстояние до центра масс) ----------
try:
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    HAS_SASA = True
except ImportError:
    HAS_SASA = False

def compute_sasa(pdb_file, chain='A'):
    if not HAS_SASA: return {}
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('prot', pdb_file)
        sr = ShrakeRupley()
        sr.compute(structure, level="R")
        sasa_dict = {}
        for model in structure:
            for ch in model:
                if ch.id != chain: continue
                for res in ch:
                    if 'CA' in res:
                        sasa_dict[res.id[1]] = res.sasa
        return sasa_dict
    except: return {}

def compute_distance_to_center(pdb_file, chain='A'):
    """Расстояние от Cα до центра масс белка (Å)."""
    if not HAS_SASA: return {}
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('prot', pdb_file)
        coords = []
        res_ids = []
        for model in structure:
            for ch in model:
                if ch.id != chain: continue
                for res in ch:
                    if 'CA' in res:
                        coords.append(res['CA'].coord)
                        res_ids.append(res.id[1])
        if not coords: return {}
        center = np.mean(coords, axis=0)
        dist_dict = {}
        for rid, c in zip(res_ids, coords):
            dist_dict[rid] = np.linalg.norm(c - center)
        return dist_dict
    except: return {}

def is_trap_exposed(sequence, pos, ref, alt, ss_map=None, pdb_file=None, sasa_thr=50.0, dist_thr=12.0):
    """
    Открытая ловушка:
      - петля (ss_map[pos] == 'C'),
      - SASA > sasa_thr,
      - расстояние до центра масс > dist_thr (на поверхности).
    """
    if ss_map is None or pos >= len(ss_map):
        hydrophobic = set('AILMFWYV')
        hydrophilic = set('RNDCQEGHKPST')
        return ref in hydrophobic and alt in hydrophilic

    if ss_map[pos] != 'C':
        return False

    if pdb_file and os.path.exists(pdb_file) and HAS_SASA:
        sasa_dict = compute_sasa(pdb_file)
        dist_dict = compute_distance_to_center(pdb_file)
        if sasa_dict and dist_dict:
            sasa = sasa_dict.get(pos + 1, 0.0)
            dist = dist_dict.get(pos + 1, 0.0)
            return (sasa > sasa_thr) and (dist > dist_thr)

    return True  # без данных — петля считается открытой
def cubic_ddG(depth):
    """Кубический закон затухания: Y = 3.61·X³ - 10.25·X² + 8.47·X - 1.26"""
    x = depth
    return 3.61 * x**3 - 10.25 * x**2 + 8.47 * x - 1.26

