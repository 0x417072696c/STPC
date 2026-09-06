# -*- coding: utf-8 -*-
"""
AllyMind 10.2 FINAL – PDB parser с поддержкой target_chain (ускорение для крио-EM).
Конфиденциально. Валентин Лебедкин, 11.07.2026.
"""
import os, urllib.request, numpy as np

CACHE_DIR = 'pdb_cache'
os.makedirs(CACHE_DIR, exist_ok=True)

def download_pdb(pdb_id):
    """Скачивает PDB-файл. Поддерживает суффикс _Chain (напр. 8G54_A)."""
    base_id = pdb_id.split('_')[0].lower()
    path = os.path.join(CACHE_DIR, f'{base_id}.pdb')
    if os.path.exists(path):
        return path
    url = f'https://files.rcsb.org/download/{base_id}.pdb'
    try:
        urllib.request.urlretrieve(url, path)
        return path
    except Exception:
        return None

def extract_aligned_ss_map(sequence, pdb_file, chain='A'):
    """Извлекает вторичную структуру из PDB-файла (HELIX/SHEET) для указанной цепи."""
    try:
        ss_map = ['C'] * len(sequence)
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('HELIX'):
                    ch = line[19].strip()
                    if ch != chain: continue
                    try:
                        start = int(line[21:25].strip())
                        end = int(line[33:37].strip())
                    except ValueError:
                        continue
                    for r in range(start, end+1):
                        idx = r - 1
                        if 0 <= idx < len(sequence):
                            ss_map[idx] = 'H'
                elif line.startswith('SHEET'):
                    ch = line[21].strip()
                    if ch != chain: continue
                    try:
                        start = int(line[22:26].strip())
                        end = int(line[33:37].strip())
                    except ValueError:
                        continue
                    for r in range(start, end+1):
                        idx = r - 1
                        if 0 <= idx < len(sequence):
                            ss_map[idx] = 'E'
        return ss_map
    except Exception:
        return None

def extract_metal_sites(pdb_file, target_chain=None):
    """Возвращает список (имя_металла, (x,y,z))."""
    metals = []
    metal_names = {"ZN","CA","MG","FE","MN","CO","NI","CU","CD","HG","K","NA","PT","AU"}
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    atom = line[12:16].strip()
                    if atom in metal_names:
                        ch = line[21:22].strip()
                        if target_chain and ch != target_chain:
                            continue
                        try:
                            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                            metals.append((atom, (x, y, z)))
                        except: pass
    except: pass
    return metals

def extract_dimer_interface_residues(pdb_file, distance_cutoff=5.0, target_chain=None):
    """Находит остатки на интерфейсе между разными цепями (быстро)."""
    residues = {}
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") and line[13:15].strip() == "CA":
                    chain = line[21:22].strip()
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    residues.setdefault(chain, []).append((res_num, (x, y, z)))
    except: pass
    interface_residues = set()
    chains = list(residues.keys())
    for i in range(len(chains)):
        for j in range(i+1, len(chains)):
            if target_chain and chains[i] != target_chain and chains[j] != target_chain:
                continue  # проверяем только пары с участием target_chain
            for r1, c1 in residues[chains[i]]:
                for r2, c2 in residues[chains[j]]:
                    if np.linalg.norm(np.array(c1) - np.array(c2)) < distance_cutoff:
                        interface_residues.add(r1)
                        interface_residues.add(r2)
    return interface_residues
