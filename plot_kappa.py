"""
plot_kappa.py – κ vs длина белка + разрешение PDB.
"""
import matplotlib.pyplot as plt
import numpy as np

# Последние значения κ (средняя глубина / J_NAT) из расчёта доменов
data = [
    # (имя, длина, κ, тип_структуры)
    ("SOD1 (БАС)", 154, 19.0454, "X-ray"),
    ("TDP-43 (БАС)", 414, 18.1964, "X-ray"),
    ("Тау (Альцгеймер)", 377, 13.2126, "Cryo-EM"),
    ("α-синуклеин (Паркинсон)", 499, 17.7999, "NMR"),
    ("Прион (БКЯ)", 253, 19.9893, "X-ray"),
    ("Паркин (ювен. Паркинсон)", 1863, 11.0055, "X-ray"),
    ("Убиквитин", 229, 0.0, "X-ray"),
    ("GFP", 238, 0.0, "X-ray"),
    ("Лизоцим", 147, 0.0, "X-ray"),
    ("Цитохром C", 105, 0.0, "X-ray"),
    ("Миоглобин", 154, 0.0, "X-ray"),
    ("Трипсин", 247, 0.0, "X-ray"),
]

names = [d[0] for d in data]
lengths = [d[1] for d in data]
kappas = [d[2] for d in data]
structures = [d[3] for d in data]

colors = {'X-ray': '#3daee9', 'Cryo-EM': '#ff5555', 'NMR': '#ffaa00'}
color_list = [colors[s] for s in structures]

fig, ax = plt.subplots(figsize=(10,6))
scatter = ax.scatter(lengths, kappas, c=color_list, s=120, edgecolors='white', linewidth=1.5, zorder=5)

# Подписи точек
for i, name in enumerate(names):
    offset = 10 if i % 2 == 0 else -15
    ax.annotate(name, (lengths[i], kappas[i]), textcoords="offset points",
                xytext=(0, offset), ha='center', fontsize=8)

# Горизонтальная линия теоретического κ = 8+5+3 = 16
ax.axhline(y=16, color='gold', linestyle='--', linewidth=2, alpha=0.8, label='κ = 8 + 5 + 3 = 16')

# Легенда по типам структур
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors['X-ray'], label='X-ray'),
                   Patch(facecolor=colors['Cryo-EM'], label='Cryo-EM'),
                   Patch(facecolor=colors['NMR'], label='NMR')]
ax.legend(handles=legend_elements, title='PDB Structure Type', loc='upper right')

ax.set_xlabel('Protein Length (residues)')
ax.set_ylabel('κ = mean depth / J_NAT')
ax.set_title('Dependence of κ on protein length and structure resolution\nD=67 generation of the Crystal')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('kappa_vs_length.png', dpi=150)
plt.show()

print("График сохранён в kappa_vs_length.png")