"""
allymind_core.py – 31 правило, контрольный ноль, сигнальный пептид исправлен.
Конфиденциально. Валентин Лебедкин, 11.07.2026.
"""
import os, sys, json, sqlite3, re
import cascade, pdb_parser, protein_fetcher, fasta_parser, vep_client

def resource_path(relative_path):
    try: return os.path.join(sys._MEIPASS, relative_path)
    except: return os.path.abspath(relative_path)

DB_CLINVAR = resource_path("clinvar.db")
DB_UNIPROT = resource_path("uniprot.db")

AA_PROPS = {
    'A': {'hydrophobicity': 1.8, 'volume': 88.6, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'R': {'hydrophobicity': -4.5, 'volume': 173.4, 'charge': 1, 'ss_bond': False, 'helix_breaker': False},
    'N': {'hydrophobicity': -3.5, 'volume': 114.1, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'D': {'hydrophobicity': -3.5, 'volume': 111.1, 'charge': -1, 'ss_bond': False, 'helix_breaker': False},
    'C': {'hydrophobicity': 2.5, 'volume': 108.5, 'charge': 0, 'ss_bond': True, 'helix_breaker': False},
    'Q': {'hydrophobicity': -3.5, 'volume': 143.8, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'E': {'hydrophobicity': -3.5, 'volume': 138.4, 'charge': -1, 'ss_bond': False, 'helix_breaker': False},
    'G': {'hydrophobicity': -0.4, 'volume': 60.1, 'charge': 0, 'ss_bond': False, 'helix_breaker': True},
    'H': {'hydrophobicity': -3.2, 'volume': 153.2, 'charge': 0.5, 'ss_bond': False, 'helix_breaker': False},
    'I': {'hydrophobicity': 4.5, 'volume': 166.7, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'L': {'hydrophobicity': 3.8, 'volume': 166.7, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'K': {'hydrophobicity': -3.9, 'volume': 168.6, 'charge': 1, 'ss_bond': False, 'helix_breaker': False},
    'M': {'hydrophobicity': 1.9, 'volume': 162.9, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'F': {'hydrophobicity': 2.8, 'volume': 189.9, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'P': {'hydrophobicity': -1.6, 'volume': 112.7, 'charge': 0, 'ss_bond': False, 'helix_breaker': True},
    'S': {'hydrophobicity': -0.8, 'volume': 89.0, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'T': {'hydrophobicity': -0.7, 'volume': 116.1, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'W': {'hydrophobicity': -0.9, 'volume': 227.8, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'Y': {'hydrophobicity': -1.3, 'volume': 193.6, 'charge': 0, 'ss_bond': False, 'helix_breaker': False},
    'V': {'hydrophobicity': 4.2, 'volume': 140.0, 'charge': 0, 'ss_bond': False, 'helix_breaker': False}
}

SOLUBLE_PROTEINS = {"TP53", "P04637", "SOD1", "P00441", "ANG", "P03950",
                    "LYSC", "P00698", "CI2", "P01053", "MBO", "UBQ",
                    "P0CG47", "P0CG48", "HSPA1A", "P08107"}

ZN_FINGER_RESIDUES = {"P04637": {176, 179, 238, 242}, "TP53": {176, 179, 238, 242}}
ONCOGENIC_MUTATIONS = {"P04637": {175: {"R": "H"}, 248: {"R": "W"}, 273: {"R": "H"}},
                        "TP53": {175: {"R": "H"}, 248: {"R": "W"}, 273: {"R": "H"}}}
ACTIVE_SITE_RESIDUES = {"1ANG": {114}}
PHOSPHO_SITES = {"S", "T", "Y"}
GLYCOSYLATION_MOTIF = re.compile(r'N[^P][ST]')
PEST_MOTIF = re.compile(r'[PEST]{3,}')
LIR_MOTIF = re.compile(r'[FWY]..[LVIFY]')
BH3_MOTIF = re.compile(r'L[^P]{3,5}[DE]')
CYCLIN_MOTIF = re.compile(r'[RK].L.{0,2}[FY]')
PCNA_MOTIF = re.compile(r'Q..[ILMV].{2}[FY]')

def load_sequence(gene):
    if os.path.exists(DB_UNIPROT):
        try:
            conn = sqlite3.connect(DB_UNIPROT)
            cur = conn.cursor()
            cur.execute("SELECT uniprot_id, sequence FROM proteins WHERE gene=?", (gene,))
            row = cur.fetchone()
            conn.close()
            if row: return row[1], row[0]
        except: pass
    try:
        data = protein_fetcher.fetch_uniprot(gene)
        if data and 'sequence' in data: return data['sequence'], gene
    except: pass
    return None, None

def get_mutations(gene, only_pathogenic=True):
    if not os.path.exists(DB_CLINVAR): return []
    try:
        conn = sqlite3.connect(DB_CLINVAR)
        cur = conn.cursor()
        query = "SELECT protein_position, ref_aa, alt_aa, clinical_significance FROM mutations WHERE gene=? AND protein_position IS NOT NULL"
        if only_pathogenic: query += " AND clinical_significance LIKE '%athogenic%'"
        rows = cur.execute(query, (gene,)).fetchall()
        conn.close()
        return [(pos, ref, alt, clin) for pos, ref, alt, clin in rows if pos is not None]
    except: return []

def determine_strategy(mutation_code, depth_old, depth_new, classification, thr, uniprot_id="", pdb_id="", sequence="", ss_map=None, pdb_file=None, metal_coordinated_residues=None, dimer_interface_residues=None, dna_contacts=None, ppi_interface=None, non_zn_coordinated=None, bfactors=None):
    pos = int(mutation_code[1:-1])
    wt = mutation_code[0]
    alt = mutation_code[-1]

    hyd_wt = AA_PROPS.get(wt, {}).get('hydrophobicity', 0.0)
    hyd_alt = AA_PROPS.get(alt, {}).get('hydrophobicity', 0.0)
    vol_wt = AA_PROPS.get(wt, {}).get('volume', 100.0)
    vol_alt = AA_PROPS.get(alt, {}).get('volume', 100.0)
    delta_hyd = abs(hyd_alt - hyd_wt) * 0.02
    delta_vol = abs(vol_alt - vol_wt) * 0.001

    if wt == alt:
        return "Stable | Control Zero Invariant"

    onc = ONCOGENIC_MUTATIONS.get(uniprot_id, {})
    if pos in onc and wt in onc[pos] and onc[pos][wt] == alt:
        return "Oncogenic | Zn loss"

    zn_set = ZN_FINGER_RESIDUES.get(uniprot_id, set())
    if pos in zn_set:
        return "Structural Zn-finger disruption"

    if pdb_id in ACTIVE_SITE_RESIDUES and pos in ACTIVE_SITE_RESIDUES[pdb_id]:
        return "Inactivation of Active Site | Monitor"

    if wt in PHOSPHO_SITES and alt not in PHOSPHO_SITES:
        return "Kinase target loss"

    if wt == 'N':
        for m in GLYCOSYLATION_MOTIF.finditer(sequence):
            if m.start() == pos - 1:
                return "Glycosylation site loss"

    if wt == 'K' and pos + 2 < len(sequence):
        return "Ubiquitination site loss"

    if pos <= 25 and uniprot_id.upper() not in SOLUBLE_PROTEINS and not pdb_file:
        return "Signal peptide cleavage loss"

    window = 20
    if pos > window and pos < len(sequence) - window:
        segment = sequence[pos-window:pos+window]
        hydrophobic = sum(1 for aa in segment if AA_PROPS.get(aa, {}).get('hydrophobicity', 0) > 1.5)
        if hydrophobic > 15:
            if uniprot_id.upper() in SOLUBLE_PROTEINS:
                return "Core Cavity Formation | Chaperone Required"
            else:
                return "TM domain disruption"

    nls_motifs = ["KKRK", "KKRP", "KKRR", "KKR"]
    for motif in nls_motifs:
        if sequence.find(motif) == pos - 1:
            return "Nuclear localization signal loss"

    if wt in ('R', 'K') and alt not in ('R', 'K'):
        if (pos > 0 and sequence[pos-1] in ('R','K')) or (pos < len(sequence)-1 and sequence[pos+1] in ('R','K')):
            return "Mitochondrial targeting signal loss"

    if wt in ('L', 'I', 'V', 'M', 'F', 'Y'):
        neighbours = sequence[max(0,pos-2):pos+2]
        hydrophobic_neigh = sum(1 for aa in neighbours if AA_PROPS.get(aa, {}).get('hydrophobicity', 0) > 2.0)
        if hydrophobic_neigh >= 3:
            return "Chaperone dependence"

    if PEST_MOTIF.search(sequence, pos-5, pos+5):
        return "Proteasomal degradation signal loss"

    if LIR_MOTIF.search(sequence, pos-4, pos+4):
        return "Autophagy signal loss"

    if BH3_MOTIF.search(sequence, pos-6, pos+6):
        return "Apoptosis regulation loss"

    if CYCLIN_MOTIF.search(sequence, pos-5, pos+5):
        return "Cell cycle control loss"

    if PCNA_MOTIF.search(sequence, pos-6, pos+6):
        return "DNA repair pathway loss"

    if wt in ('Y', 'F', 'W') and alt not in ('Y', 'F', 'W'):
        if depth_new > 1.5:
            return "Aromatic core disruption"

    if dimer_interface_residues and pos in dimer_interface_residues:
        return "Dimer Interface Disruption"

    if metal_coordinated_residues and pos in metal_coordinated_residues:
        return "Chaperone | Loss of Metal Coordination / Core Collapse"

    if pdb_file and cascade.detect_steric_clash(sequence, pos, alt, ss_map, pdb_file):
        return "Steric Clash | Core Rupture (Pathogenic Volume Expansion)"

    if cascade.surface_charge_inversion(wt, alt, classification == 'open'):
        return "Aggregation Risk | Surface Charge Inversion"

    if wt == 'C' and alt != 'C':
        return "Disulfide bond disruption"

    if AA_PROPS.get(alt, {}).get('helix_breaker', False):
        return "Helix breaker | Structural disruption"

    if bfactors and pos < len(bfactors) and bfactors[pos] > 60 and depth_new > 2.0 and delta_vol > 0.1:
        return "Allosteric Hotspot | Dynamic Coupling Disrupted"

    if dna_contacts and pos+1 in dna_contacts and depth_new > 1.5:
        return "DNA/RNA Binding Site Disrupted"

    if ppi_interface and pos+1 in ppi_interface:
        return "Protein-Protein Interface Disruption"

    if ss_map and pos < len(ss_map) and ss_map[pos] == 'C' and depth_new < 1.0:
        if (AA_PROPS.get(wt, {}).get('charge', 0) * AA_PROPS.get(alt, {}).get('charge', 0) < 0):
            return "IDR Functional Modulation (charge inversion)"

    proteolytic_motifs = ["DEVD", "LEHD", "IETD", "WEHD"]
    for motif in proteolytic_motifs:
        if sequence.find(motif) == pos - 1:
            return "Proteolytic Site Altered"

    if ss_map and pos > 2 and pos < len(sequence)-3:
        left_ss = ss_map[pos-2] if pos-2 >= 0 else 'C'
        right_ss = ss_map[pos+2] if pos+2 < len(ss_map) else 'C'
        if left_ss in ('H','E') and right_ss in ('H','E') and ss_map[pos] == 'C':
            return "Linker Disruption | Domain orientation affected"

    if classification == 'open' and abs(delta_hyd) > 0.2 and abs(delta_vol) > 0.1:
        return "Filament Assembly Defect"

    if non_zn_coordinated and pos+1 in non_zn_coordinated:
        return "Metal Coordination Loss (Mg/Ca/Fe)"

    if depth_old > 0 and depth_new < depth_old * (0.5 if classification == "hidden" else 0.8):
        return "Stable (no action required)"
    elif depth_old > 0 and depth_new > depth_old * 1.2:
        return f"Chaperone (J>={thr:.4f})" if thr > 0 else "Chaperone / proteolysis"
    else:
        return "Small molecule" if classification == "open" else "Monitor"

def analyze_mutation(sequence, mutation_code, pdb_id=None, ss_map=None):
    if len(mutation_code) >= 3 and mutation_code[1:-1].isdigit():
        wt = mutation_code[0]
        pos = int(mutation_code[1:-1]) - 1
        alt = mutation_code[-1]
    else:
        return {"error": "Неверный формат мутации"}

    N = len(sequence)
    if ss_map is None and pdb_id:
        pdb_file = pdb_parser.download_pdb(pdb_id)
        if pdb_file:
            ss_map = pdb_parser.extract_aligned_ss_map(sequence, pdb_file)

    J_local, dom_rad = cascade.secondary_structure_J(sequence, pos, ss_map)
    hyd_wt = AA_PROPS.get(wt, {}).get('hydrophobicity', 0.0)
    hyd_alt = AA_PROPS.get(alt, {}).get('hydrophobicity', 0.0)
    vol_wt = AA_PROPS.get(wt, {}).get('volume', 100.0)
    vol_alt = AA_PROPS.get(alt, {}).get('volume', 100.0)
    delta_hyd = abs(hyd_alt - hyd_wt) * 0.02
    delta_vol = abs(vol_alt - vol_wt) * 0.001
    ss_bonus = 0.0
    helix_penalty = 0.0
    if alt == 'C' and AA_PROPS.get(alt, {}).get('ss_bond', False) and wt != 'C':
        ss_bonus = 0.05
    if ss_map and 0 <= pos < len(ss_map) and ss_map[pos] == 'H':
        if AA_PROPS.get(alt, {}).get('helix_breaker', False):
            helix_penalty = 0.1

    N_ref = 164
    size_factor = (N_ref / N) ** 0.5
    packing_density = 0.5
    if ss_map and 0 <= pos < len(ss_map):
        packing_density = 0.85 if ss_map[pos] in ('H', 'E') else 0.5
    delta_vol_eff = delta_vol * (1.0 - packing_density) if packing_density > 0.75 else delta_vol * 1.5
    local_ddg = round(delta_hyd + delta_vol_eff + ss_bonus + helix_penalty, 4)
    J_local += delta_hyd * 0.1 + delta_vol_eff * 0.05 + ss_bonus + helix_penalty

    start = max(0, pos - dom_rad)
    end = min(N, pos + dom_rad + 1)
    dom_len = end - start
    healthy = N - dom_len
    native = {'nh': healthy, 'nd': dom_len, 'J': cascade.J_NAT, 'N': N}
    defect = {'nh': healthy, 'nd': dom_len, 'J': J_local, 'N': N}
    hole   = {'nh': healthy, 'nd': 0, 'J': cascade.J_NAT, 'N': N}
    F_defect_dimer = cascade.dimer_free_energy(defect, defect, cascade.J_NAT, dom_len)
    F_mixed_dimer  = cascade.dimer_free_energy(defect, hole, cascade.J_NAT, dom_len)
    trap = F_defect_dimer < F_mixed_dimer
    depth = (F_mixed_dimer - F_defect_dimer) * size_factor if trap else local_ddg * size_factor

    try:
        exposed = cascade.is_trap_exposed(sequence, pos, wt, alt, ss_map)
    except:
        exposed = False
    depth += cascade.buried_charge_penalty(alt, exposed)
    if wt in ('Y', 'F', 'W') and alt not in ('Y', 'F', 'W'):
        depth += 1.5

    if wt == alt:
        depth = 0.0

    dG_bind = round(delta_hyd * 0.5 + delta_vol_eff * 0.25 + ss_bonus * 0.5 + helix_penalty * 0.5, 4)
    try:
        k_agg, _ = cascade.estimate_aggregation_rate(sequence, pos, wt, alt)
    except:
        k_agg = 0.0

    if ss_map and 0 <= pos < len(ss_map):
        ss_type = ss_map[pos]
    else:
        pred = cascade.predict_ss_chou_fasman(sequence)
        ss_type = pred[pos] if pos < len(pred) else 'C'
    thr = round(0.115 + 0.002 * dom_len, 4) if ss_type == 'H' else round(0.125 + 0.003 * dom_len, 4) if ss_type == 'E' else round(0.12 + 0.001 * dom_len, 4)
    ss_label = {'E': 'β-sheet', 'H': 'α-helix', 'C': 'loop'}.get(ss_type, 'loop')

    return {
        "mutation": f"{wt}{pos+1}{alt}",
        "J_local": J_local,
        "trap": trap,
        "depth": round(depth, 4),
        "secondary_structure": ss_label,
        "kinetics": {"delta_G": dG_bind, "k_agg": round(k_agg, 2)},
        "therapeutic_threshold": thr,
        "delta_delta_G": local_ddg,
        "chemical_components": {
            "delta_hyd": round(delta_hyd, 4),
            "delta_vol": round(delta_vol_eff, 4),
            "ss_bonus": round(ss_bonus, 4),
            "helix_penalty": round(helix_penalty, 4)
        },
        "energies": {},
        "vep_scores": None
    }
