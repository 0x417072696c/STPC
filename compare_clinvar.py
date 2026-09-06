"""
compare_clinvar.py KAPPA — только патогенные белки, κ = средняя глубина / J_NAT.
"""
import allymind_core, cascade, pdb_parser
import numpy as np, os

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")
WINDOW = 5

targets = [
    ("P00441", "2C9V", "SOD1 (БАС)"),
    ("Q13148", "4IOO", "TDP-43 (БАС)"),
    ("P10636", "8G54_A", "Тау (Альцгеймер)"),
    ("P37840", "1XQ8", "α-синуклеин (Паркинсон)"),
    ("P04156", "1A2V", "Прион (БКЯ)"),
    ("P38398", "1OQZ", "Паркин (ювен. Паркинсон)"),
]

def find_core_domain(ss_map):
    best_start, best_end = 0, 0
    current_start = None
    for i, s in enumerate(ss_map):
        if s in ('H', 'E'):
            if current_start is None: current_start = i
        else:
            if current_start is not None:
                if i - current_start > best_end - best_start:
                    best_start, best_end = current_start, i
                current_start = None
    if current_start is not None:
        if len(ss_map) - current_start > best_end - best_start:
            best_start, best_end = current_start, len(ss_map)
    if best_end - best_start == 0:
        L = len(ss_map); best_start, best_end = int(L*0.33), int(L*0.66)
    return max(0, best_start-WINDOW), min(len(ss_map), best_end+WINDOW)

print(f"{'Белок':<30} {'Зона (а.о.)':<12} {'Средняя глубина':<16} {'κ = глубина / J_NAT':<20}")
print("-" * 80)

kappa_list = []

for uniprot, pdb_id, name in targets:
    seq, uid = allymind_core.load_sequence(uniprot)
    if not seq: continue
    pdb_file = pdb_parser.download_pdb(pdb_id)
    ss_map = None
    if pdb_file and os.path.exists(pdb_file):
        try: ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file, chain='A')
        except: pass
    if ss_map is None: continue

    zone_start, zone_end = find_core_domain(ss_map)
    depths = []

    for pos in range(zone_start+1, zone_end+1):
        wt = seq[pos-1]
        for alt in AMINO_ACIDS:
            if alt == wt: continue
            report = allymind_core.analyze_mutation(seq, f"{wt}{pos}{alt}", ss_map=ss_map)
            depth = report.get('depth',0)
            # Учитываем только ловушки с глубиной > 0.5 эВ (потенциально опасные)
            if depth > 0.5:
                depths.append(depth)

    if depths:
        mean_depth = np.mean(depths)
    else:
        mean_depth = 0.0

    kappa = mean_depth / cascade.J_NAT
    kappa_list.append(kappa)
    print(f"{name:<30} {zone_end-zone_start:<12} {mean_depth:<16.4f} {kappa:<20.4f}")

if kappa_list:
    final_kappa = np.mean(kappa_list)
    print(f"\nСредний κ по всем патогенным белкам = {final_kappa:.4f}")
    print(f"Теоретическое κ = 8/3 = {8/3:.4f}")