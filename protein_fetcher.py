# -*- coding: utf-8 -*-
"""
protein_fetcher.py – загрузка последовательностей из UniProt (AllyMind 9.0).
Конфиденциально. Валентин Лебедкин, 05.07.2026.
"""
import os
import urllib.request
import urllib.error
import json

CACHE_DIR = "cache"

def _cache_path(identifier):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, f"{identifier}.json")

def fetch_uniprot(gene_or_uniprot):
    """
    Загружает информацию о белке из UniProt REST API.
    Возвращает словарь с ключами 'sequence', 'uniprot_id', 'gene_symbol' или None при ошибке.
    """
    identifier = gene_or_uniprot.strip().upper()
    cache_file = _cache_path(identifier)

    # пробуем загрузить из кэша
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
                if data and 'sequence' in data:
                    return data
        except Exception:
            pass

    # запрос к UniProt
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{identifier}+OR+accession:{identifier}&fields=accession,sequence,protein_name&format=json&size=1"
    try:
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'AllyMind/9.0')
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())
            if data.get('results') and len(data['results']) > 0:
                result = data['results'][0]
                entry = {
                    'sequence': result.get('sequence', {}).get('value', ''),
                    'uniprot_id': result.get('primaryAccession', identifier),
                    'gene_symbol': gene_or_uniprot
                }
                # кэшируем
                with open(cache_file, 'w', encoding='utf-8') as f:
                    json.dump(entry, f)
                return entry
    except urllib.error.HTTPError as e:
        print(f"HTTP Error {e.code}: {e.reason}")
    except Exception as e:
        print(f"Ошибка при загрузке {identifier}: {e}")

    return None
