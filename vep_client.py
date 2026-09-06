# -*- coding: utf-8 -*-
"""
AllyMind 4.0 – Онлайн-клиент Ensembl VEP.
Возвращает оценки Polyphen, SIFT, CADD для мутации.
Конфиденциально. Валентин Лебедкин, 03.07.2026.
"""
import json, urllib.request, urllib.error

VEP_URL = "https://rest.ensembl.org/vep/human/region/{}:{}-{}/{}?content-type=application/json"

def get_vep_scores(chrom, pos, ref, alt):
    """
    Получает оценки Polyphen, SIFT и CADD из Ensembl VEP.
    Аргументы:
        chrom – хромосома (например, '21')
        pos – позиция в геноме (1-based)
        ref – референсный аллель
        alt – альтернативный аллель
    Возвращает словарь с ключами polyphen_score, polyphen_pred, sift_score, sift_pred, cadd_phred.
    Возвращает None, если запрос не удался.
    """
    url = VEP_URL.format(chrom, pos, pos, f"{ref}/{alt}")
    req = urllib.request.Request(url, headers={
        'User-Agent': 'AllyMind/4.0',
        'Accept': 'application/json'
    })
    try:
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode('utf-8'))
    except (urllib.error.URLError, json.JSONDecodeError) as e:
        print(f"VEP error: {e}")
        return None

    if not data or 'transcript_consequences' not in data[0]:
        return None

    # Берём первый канонический транскрипт
    for tc in data[0]['transcript_consequences']:
        polyphen_score = tc.get('polyphen_score')
        polyphen_pred = tc.get('polyphen_prediction')
        sift_score = tc.get('sift_score')
        sift_pred = tc.get('sift_prediction')
        cadd_phred = tc.get('cadd_phred')

        if polyphen_score is not None or sift_score is not None:
            return {
                'polyphen_score': polyphen_score,
                'polyphen_pred': polyphen_pred,
                'sift_score': sift_score,
                'sift_pred': sift_pred,
                'cadd_phred': cadd_phred
            }
    return None
