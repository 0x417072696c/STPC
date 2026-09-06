"""
double_mutant.py – анализ двойной мутации (A4V+I113T)
Использование: python double_mutant.py P00441 1HLN A4V I113T
"""

import sys
import allymind_core
import cascade
import pdb_parser

def analyze_double(seq, mut1, mut2, ss_map):
    """Считает две мутации и суммирует эффекты."""
    report1 = allymind_core.analyze_mutation(seq, mut1, ss_map=ss_map)
    report2 = allymind_core.analyze_mutation(seq, mut2, ss_map=ss_map)

    combined = {
        "mutation": f"{report1['mutation']} + {report2['mutation']}",
        "depth": round(report1['depth'] + report2['depth'], 4),
        "J_local": round((report1['J_local'] + report2['J_local']) / 2, 4),
        "trap": report1['trap'] or report2['trap'],
        "ΔΔG": round(report1['kinetics']['delta_G'] + report2['kinetics']['delta_G'], 4),
    }
    return report1, report2, combined


if __name__ == "__main__":
    uniprot = sys.argv[1] if len(sys.argv) > 1 else "P00441"
    pdb_id = sys.argv[2] if len(sys.argv) > 2 else "1HLN"
    mut1 = sys.argv[3] if len(sys.argv) > 3 else "A4V"
    mut2 = sys.argv[4] if len(sys.argv) > 4 else "I113T"

    seq, uid = allymind_core.load_sequence(uniprot)
    if not seq:
        print(f"Не удалось загрузить {uniprot}")
        sys.exit(1)

    pdb_file = pdb_parser.download_pdb(pdb_id)
    if pdb_file:
        ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file)
    else:
        ss_map = None

    r1, r2, combined = analyze_double(seq, mut1, mut2, ss_map)

    print(f"SOD1 ({uid}), длина {len(seq)} а.о., PDB {pdb_id}")
    print(f"Мутации: {mut1} + {mut2}")
    print(f"  {mut1}: глубина {r1['depth']:.3f} эВ, ΔΔG {r1['kinetics']['delta_G']:.4f} эВ, структура {r1['secondary_structure']}")
    print(f"  {mut2}: глубина {r2['depth']:.3f} эВ, ΔΔG {r2['kinetics']['delta_G']:.4f} эВ, структура {r2['secondary_structure']}")
    print(f"  Совместно: глубина {combined['depth']:.3f} эВ, ΔΔG {combined['ΔΔG']:.4f} эВ, ловушка: {'да' if combined['trap'] else 'нет'}")