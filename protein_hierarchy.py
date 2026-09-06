"""
Иерархия белков по когерентной уязвимости.
Сравнивает здоровые и патогенные белки.
"""
import allymind_core, cascade, pdb_parser
import numpy as np

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")

# Набор белков: (UniProt, PDB, название, категория)
proteins = [
    ("P00441", "1HLN", "SOD1 (БАС)", "патогенный"),
    ("Q13148", "4IOO", "TDP-43 (БАС)", "патогенный"),
    ("P10636", "5XLG_A", "Тау (Альцгеймер)", "патогенный"),
    ("P0CG47", "1UBQ", "Убиквитин", "здоровый"),
    ("P42212", "1GFL", "GFP", "здоровый"),
    ("P00698", "1A2U", "Лизоцим", "здоровый"),
]

def analyze_protein(uniprot, pdb_id):
    seq, uid = allymind_core.load_sequence(uniprot)
    if not seq: return None
    pdb_file = pdb_parser.download_pdb(pdb_id)
    ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file, chain='A') if pdb_file else None

    total_ddG = 0.0
    total_depth = 0.0
    critical = 0
    count_open = 0

    for pos in range(1, len(seq)+1):
        wt = seq[pos-1]
        for alt in AMINO_ACIDS:
            if alt == wt: continue
            mut = f"{wt}{pos}{alt}"
            report = allymind_core.analyze_mutation(seq, mut, ss_map=ss_map)
            depth = report.get('depth', 0)
            ddG = report.get('kinetics', {}).get('delta_G', 0) or 0
            exposed = cascade.is_trap_exposed(seq, pos-1, wt, alt, ss_map)
            if exposed and depth > 0.5:
                count_open += 1
                total_ddG += ddG
                total_depth += depth
                if ddG > 0.05:
                    critical += 1

    if count_open == 0:
        return None

    avg_depth = total_depth / count_open
    sum_ddG = total_ddG
    return {
        "uniprot": uniprot,
        "name": uid,
        "length": len(seq),
        "open_traps": count_open,
        "critical": critical,
        "sum_ddG": f"{sum_ddG:.4f}",
        "avg_depth": f"{avg_depth:.4f}"
    }

if __name__ == "__main__":
    print(f"{'Белок':<20} {'Длина':<6} {'Откр. лов.':<10} {'Крит. поз.':<10} {'Сумм. ΔΔG':<12} {'Сред. глуб.':<12}")
    print("-" * 75)
    results = []
    for uniprot, pdb, name, category in proteins:
        print(f"Анализ {name}...")
        r = analyze_protein(uniprot, pdb)
        if r:
            r["category"] = category
            r["display_name"] = name
            results.append(r)
            print(f"{name:<20} {r['length']:<6} {r['open_traps']:<10} {r['critical']:<10} {r['sum_ddG']:<12} {r['avg_depth']:<12}")
        else:
            print(f"{name:<20} — данные отсутствуют")
    
    print("\nИерархия по интегральной уязвимости (суммарный ΔΔG):")
    results.sort(key=lambda x: float(x['sum_ddG']), reverse=True)
    for i, r in enumerate(results):
        print(f"  {i+1}. {r['display_name']} ({r['category']}) — ∑ΔΔG = {r['sum_ddG']} эВ, критических позиций: {r['critical']}")