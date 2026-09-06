"""
neuro_batch.py – анализ TDP-43, Тау, Альфа-синуклеина
Использование: python neuro_batch.py
"""
import allymind_core, cascade, pdb_parser

mutations = [
    ("Q13148", "4IOO", "A315T"),   # TDP-43
    ("P10636", "8G54", "P301L"),   # Тау
    ("P37840", "1XQ8", "A53T"),    # Альфа-синуклеин
]

print(f"{'Белок':<20} {'Мутация':<8} {'Глубина (HB)':<14} {'ΔΔG (HB)':<10} {'Порог J':<10} {'Структура':<12} {'Класс':<10} {'Стратегия':<25}")
print("-" * 120)

for uniprot, pdb_id, mut in mutations:
    seq, uid = allymind_core.load_sequence(uniprot)
    if not seq:
        print(f"{uniprot:<20} {'—':<8} {'не загружен':<14}")
        continue

    pdb_file = pdb_parser.download_pdb(pdb_id)
    ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file) if pdb_file else None

    wt, pos, alt = mut[0], int(mut[1:-1]) - 1, mut[-1]
    report = allymind_core.analyze_mutation(seq, mut, ss_map=ss_map)
    depth = report['depth']
    ddG = report.get('kinetics', {}).get('delta_G', 0) or 0
    thr = report.get('therapeutic_threshold', 0) or 0
    ss = report.get('secondary_structure', 'loop') or 'loop'
    exposed = cascade.is_trap_exposed(seq, pos, wt, alt, ss_map)
    classification = "open" if exposed else "hidden"

    # Простейшая стратегия
    if classification == "open":
        strategy = "Small molecule" if depth < 1.0 else "Chaperone"
    else:
        strategy = "Monitor" if depth < thr else "Chaperone (hidden)"

    print(f"{uid:<20} {mut:<8} {depth:<14.3f} {ddG:<10.4f} {thr:<10.4f} {ss:<12} {classification:<10} {strategy:<25}")