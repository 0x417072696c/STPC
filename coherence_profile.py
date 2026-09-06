"""
Когерентный профиль белка v2 — по глубине ловушки (depth), а не по J_local.
Выделяет реальные патогенные позиции.
"""
import allymind_core, cascade, pdb_parser
import numpy as np
import matplotlib.pyplot as plt
import csv, time, sys

UNIPROT = sys.argv[1] if len(sys.argv) > 1 else "P00441"
PDB_ID  = sys.argv[2] if len(sys.argv) > 2 else "1HLN"
OUT_CSV = f"{UNIPROT}_depth_profile.csv"
OUT_PNG = f"{UNIPROT}_depth_profile.png"

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")

print(f"Загрузка {UNIPROT}...")
seq, uid = allymind_core.load_sequence(UNIPROT)
if not seq:
    print("Ошибка загрузки последовательности.")
    sys.exit(1)

print(f"Длина белка: {len(seq)} а.о.")
print(f"Загрузка PDB {PDB_ID}...")
pdb_file = pdb_parser.download_pdb(PDB_ID)
ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file, chain='A') if pdb_file else None

positions, depth_max, depth_mean, structures = [], [], [], []
print(f"Анализ {len(seq)} позиций × 19 замен...")
start_time = time.time()

for pos in range(1, len(seq) + 1):
    wt = seq[pos - 1]
    depths = []
    for alt in AMINO_ACIDS:
        if alt == wt:
            continue
        report = allymind_core.analyze_mutation(seq, f"{wt}{pos}{alt}", ss_map=ss_map)
        depths.append(report.get('depth', 0.0))
    
    positions.append(pos)
    depth_max.append(np.max(depths) if depths else 0.0)
    depth_mean.append(np.mean(depths) if depths else 0.0)
    structures.append(ss_map[pos - 1] if ss_map else 'C')

elapsed = time.time() - start_time

# Сохраняем CSV
with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["pos", "wt", "structure", "depth_max", "depth_mean"])
    for i, pos in enumerate(positions):
        writer.writerow([pos, seq[pos-1], structures[i],
                         f"{depth_max[i]:.4f}", f"{depth_mean[i]:.4f}"])

# График
fig, ax = plt.subplots(figsize=(14, 6))
pos_arr = np.array(positions)
ax.plot(pos_arr, depth_max, 'r-', linewidth=1.2, label='Макс. глубина ловушки (худшая замена)', alpha=0.9)
ax.fill_between(pos_arr, 0, depth_max, alpha=0.2, color='red')
ax.plot(pos_arr, depth_mean, 'orange', linewidth=0.8, label='Средняя глубина ловушки', alpha=0.7)
# Подпишем топ-5 позиций
top5_idx = np.argsort(depth_max)[-5:]
for i in top5_idx:
    ax.annotate(f"{seq[positions[i]-1]}{positions[i]}",
                (positions[i], depth_max[i]),
                textcoords="offset points", xytext=(0, 10), ha='center',
                fontsize=8, color='darkred', fontweight='bold')

ax.set_xlabel('Позиция в белке')
ax.set_ylabel('Глубина ловушки (эВ)')
ax.set_title(f'Профиль уязвимости белка: {UNIPROT} ({uid})\nПроверено {len(seq)*19} мутаций за {elapsed:.1f} с')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=150)
plt.show()

print(f"Готово. Профиль сохранён в {OUT_CSV} и {OUT_PNG}.")
print("Топ-5 самых уязвимых позиций (макс. глубина ловушки):")
for i in top5_idx[::-1]:
    print(f"  Поз. {positions[i]} ({seq[positions[i]-1]}): глубина = {depth_max[i]:.4f} эВ, структура = {structures[i]}")