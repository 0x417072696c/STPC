import os
import requests
from tqdm import tqdm

def download_uniprot():
    # Надежная HTTP-ссылка на актуальный релиз Swiss-Prot FASTA
    url = "https://uniprot.org"
    
    # Явно задаем имя файла в текущей папке
    local_filename = "uniprot_sprot.fasta.gz"
    
    print(f"Начало скачивания: {local_filename}")
    
    # Потоковый запрос
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    # Размер файла для прогресс-бара
    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024 * 1024  # Чтение блоками по 1 МБ
    
    with open(local_filename, 'wb') as file, tqdm(
        desc=local_filename,
        total=total_size,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as progress_bar:
        for data in response.iter_content(block_size):
            size = file.write(data)
            progress_bar.update(size)
            
    print(f"\nГотово! Файл сохранен в: {os.path.abspath(local_filename)}")

if __name__ == "__main__":
    download_uniprot()