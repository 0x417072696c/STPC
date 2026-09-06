# -*- coding: utf-8 -*-
"""
AllyMind v0.1 — Парсер FASTA и предсказатель доменов (fasta_parser.py)
Полностью автономный. Читает FASTA, определяет тип домена.
Конфиденциально. Валентин Лебедкин, 02.07.2026.
"""

import re

# Правила для предсказания вторичной структуры
ALPHA_FAVOR = set('LEQAMKHR')
BETA_FAVOR  = set('VILFYW')
HYDROPHOBIC = set('VILFYWMA')
CYSTEINE    = set('C')
CHARGED     = set('KRDE')

DOMAIN_ENERGIES = {
    'alpha': 0.110,
    'beta': 0.125,
    'hydrophobic': 0.09,
    'disulfide': 0.20,
    'salt_bridge': 0.15,
    'mixed': 0.10,
    'unknown': 0.08
}


def read_fasta(filename):
    """
    Читает FASTA-файл. Пробует UTF-8, если не получается — UTF-16.
    Возвращает список кортежей (идентификатор, последовательность).
    """
    records = []
    current_id = None
    current_seq = []
    
    # Пробуем UTF-8, потом UTF-16
    encodings = ['utf-8', 'utf-16-le', 'utf-16-be']
    lines = None
    for enc in encodings:
        try:
            with open(filename, 'r', encoding=enc) as f:
                lines = f.readlines()
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
    
    if lines is None:
        raise ValueError(f"Не удалось прочитать файл {filename} в кодировках {encodings}")
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_id is not None:
                records.append((current_id, ''.join(current_seq)))
            current_id = line[1:].split()[0]  # первый идентификатор
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_id is not None:
        records.append((current_id, ''.join(current_seq)))
    
    return records


def classify_domain_type(sequence):
    """Определяет тип домена по аминокислотному составу."""
    seq_set = set(sequence)
    alpha_count = len(seq_set & ALPHA_FAVOR)
    beta_count  = len(seq_set & BETA_FAVOR)
    hydro_count = len(seq_set & HYDROPHOBIC)
    cys_count   = len(seq_set & CYSTEINE)
    charged_count = len(seq_set & CHARGED)
    
    total = len(sequence)
    if total == 0:
        return 'unknown'
    
    if cys_count >= 2:
        return 'disulfide'
    if charged_count >= total * 0.3:
        return 'salt_bridge'
    if hydro_count >= total * 0.5:
        return 'hydrophobic'
    if alpha_count > beta_count:
        return 'alpha'
    elif beta_count > alpha_count:
        return 'beta'
    else:
        return 'mixed'


def extract_domains(sequence, window_size=20, step=5):
    """
    Сканирует последовательность скользящим окном и возвращает список доменов.
    Каждый домен — словарь с ключами: start, end, sequence, type, energy.
    """
    domains = []
    L = len(sequence)
    for start in range(0, L - window_size + 1, step):
        end = start + window_size
        subseq = sequence[start:end]
        dtype = classify_domain_type(subseq)
        energy = DOMAIN_ENERGIES.get(dtype, 0.08)
        domains.append({
            'start': start,
            'end': end,
            'sequence': subseq,
            'type': dtype,
            'energy': energy
        })
    return domains


def get_domain_summary(domains):
    """Краткая текстовая сводка по доменам."""
    lines = []
    for d in domains:
        lines.append(f"  [{d['start']:4d}-{d['end']:4d}] {d['type']:15s} J={d['energy']:.3f} эВ  seq={d['sequence'][:20]}...")
    return "\n".join(lines)


# ------------------------------------------------------------
# Тест при запуске модуля напрямую
# ------------------------------------------------------------
if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
        records = read_fasta(fasta_file)
        for rec_id, seq in records:
            print(f"\n=== {rec_id} (длина {len(seq)}) ===")
            domains = extract_domains(seq, window_size=20, step=10)
            print(get_domain_summary(domains))
    else:
        # Встроенный тест
        test_seq = "MANLGCWMLVLFVATWSDLGLCKKRPKPGG"
        print(f"Тестовая последовательность: {test_seq}")
        domains = extract_domains(test_seq, window_size=10, step=5)
        print(get_domain_summary(domains))
