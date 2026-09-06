"""
sars2_screening.py – скрининг спайк-белка SARS-CoV-2 (P0DTC2).
Предсказывает неизвестные уязвимые позиции.
"""
import allymind_core, cascade, pdb_parser
import numpy as np, matplotlib.pyplot as plt, csv, time

UNIPROT = "P0DTC2"
PDB_ID = "6VXX"  # закрытая конформация спайка
OUTPUT_CSV = "sars2_profile.csv"
OUTPUT_PNG = "sars2_vulnerability_profile.png"

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")

print(f"Загрузка {UNIPROT}...")
seq, uid = allymind_core.load_sequence(UNIPROT)
if not seq:
    print("Ошибка загрузки последовательности.")
    exit()

print(f"Длина белка: {len(seq)} а.о.")
print(f"Загрузка PDB {PDB_ID}...")
pdb_file = pdb_parser.download_pdb(PDB_ID)
ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file, chain='A') if pdb_file else None

positions, depth_max, depth_mean, worst_mut, structures = [], [], [], [], []
total_mutations = 0
start_time = time.time()

print(f"Анализ {len(seq)} позиций × 19 замен...")
for pos in range(1, len(seq)+1):
    wt = seq[pos-1]
    depths = []
    for alt in AMINO_ACIDS:
        if alt == wt:
            continue
        total_mutations += 1
        report = allymind_core.analyze_mutation(seq, f"{wt}{pos}{alt}", ss_map=ss_map)
        depths.append(report.get('depth', 0.0))

    positions.append(pos)
    depth_max.append(np.max(depths) if depths else 0.0)
    depth_mean.append(np.mean(depths) if depths else 0.0)
    worst_idx = np.argmax(depths) if depths else 0
    worst_alt = AMINO_ACIDS[worst_idx] if depths else wt
    worst_mut.append(f"{wt}{pos}{worst_alt}")
    structures.append(ss_map[pos-1] if ss_map else 'C')

elapsed = time.time() - start_time

# Сохраняем CSV
with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["pos", "wt", "worst_mutation", "structure", "depth_max", "depth_mean"])
    for i, pos in enumerate(positions):
        writer.writerow([pos, seq[pos-1], worst_mut[i], structures[i],
                         f"{depth_max[i]:.4f}", f"{depth_mean[i]:.4f}"])

# Строим график
fig, ax = plt.subplots(figsize=(14, 6))
pos_arr = np.array(positions)
ax.plot(pos_arr, depth_max, 'r-', linewidth=1.2, label='Максимальная глубина ловушки', alpha=0.9)
ax.fill_between(pos_arr, 0, depth_max, alpha=0.15, color='red')
ax.plot(pos_arr, depth_mean, 'orange', linewidth=0.8, label='Средняя глубина ловушки', alpha=0.7)

# Подписываем топ-10 позиций
top10_idx = np.argsort(depth_max)[-10:]
for i in top10_idx:
    ax.annotate(f"{seq[positions[i]-1]}{positions[i]}",
                (positions[i], depth_max[i]),
                textcoords="offset points", xytext=(0, 10), ha='center',
                fontsize=8, color='darkred', fontweight='bold')

ax.set_xlabel('Позиция в белке')
ax.set_ylabel('Глубина ловушки (эВ)')
ax.set_title(f'Профиль уязвимости спайк-белка SARS-CoV-2 ({uid})\nПроверено {total_mutations} мутаций за {elapsed:.1f} с')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=150)
plt.show()

print(f"Готово. Профиль сохранён в {OUTPUT_CSV} и {OUTPUT_PNG}.")
print(f"Проверено {total_mutations} мутаций за {elapsed:.1f} секунд.")
print("Топ-10 самых уязвимых позиций (макс. глубина ловушки):")
for i in top10_idx[::-1]:
    print(f"  Поз. {positions[i]} ({seq[positions[i]-1]}): глубина = {depth_max[i]:.4f} эВ, структура = {structures[i]}")