"""
Скрининг MAPT с балансным фильтром.
Только открытые ловушки с depth > 0.5 эВ.
"""
import allymind_core, cascade, pdb_parser
import csv, time

UNIPROT = "P10636"
PDB_ID = "5XLG_A"
OUTPUT = "mapt_balanced.csv"

ALL_AA = list("ARNDCQEGHILKMFPSTWYV")

print(f"Загрузка {UNIPROT}...")
seq, uid = allymind_core.load_sequence(UNIPROT)
if not seq:
    print("Не удалось загрузить последовательность.")
    exit()

print(f"Длина белка: {len(seq)} а.о.")
print(f"Загрузка PDB {PDB_ID}...")
pdb_file = pdb_parser.download_pdb(PDB_ID)
ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file) if pdb_file else None

results = []
total = 0
start_time = time.time()

for pos in range(1, len(seq)+1):
    wt = seq[pos-1]
    for alt in ALL_AA:
        if alt == wt:
            continue
        total += 1
        mut = f"{wt}{pos}{alt}"
        report = allymind_core.analyze_mutation(seq, mut, ss_map=ss_map)
        depth = report['depth']
        exposed = cascade.is_trap_exposed(seq, pos-1, wt, alt, ss_map)
        
        if exposed and depth > 0.5:
            ddG = report.get('kinetics', {}).get('delta_G', 0) or 0
            thr = report.get('therapeutic_threshold', 0) or 0
            ss = report.get('secondary_structure', '?') or '?'
            strategy = "Small molecule" if depth < 1.0 else "Chaperone"
            results.append((depth, pos, wt, alt, ddG, thr, ss, strategy))

results.sort(reverse=True)
elapsed = time.time() - start_time

with open(OUTPUT, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["pos", "wt", "alt", "depth", "ddG", "threshold_J", "structure", "strategy"])
    for depth, pos, wt, alt, ddG, thr, ss, strategy in results:
        writer.writerow([pos, wt, alt, f"{depth:.4f}", f"{ddG:.4f}", f"{thr:.4f}", ss, strategy])

print(f"Готово. Проверено {total} мутаций за {elapsed:.1f} секунд.")
print(f"Найдено {len(results)} открытых ловушек с глубиной > 0.5 эВ. Результат сохранён в {OUTPUT}.")
if results:
    print("\nТоп-5 самых глубоких ловушек:")
    for i, (depth, pos, wt, alt, ddG, thr, ss, strategy) in enumerate(results[:5]):
        print(f"  {i+1}. {wt}{pos}{alt} — глубина {depth:.4f} эВ, структура {ss}, стратегия {strategy}")