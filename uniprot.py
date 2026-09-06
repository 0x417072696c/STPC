import gzip
import xml.etree.ElementTree as ET

file_path = r"C:\Users\Colorful X15 AT23\Desktop\СТПК\Levi\AllyMind 10_2\uniprot_sprot.xml.gz"
print("Запуск исправленного парсера UniProt...")

try:
    with gzip.open(file_path, 'rb') as f:
        # Используем iterparse для экономии RAM
        context = ET.iterparse(f, events=('end',))
        count = 0
        
        for event, elem in context:
            # Очищаем имя тега от префиксов (пространств имен)
            tag = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
            
            if tag == 'entry':
                count += 1
                uid = "Unknown"
                sequence = ""
                
                # Ищем accession и sequence внутри текущего entry без привязки к урлам
                for child in elem:
                    child_tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                    if child_tag == 'accession' and uid == "Unknown":
                        uid = child.text
                    elif child_tag == 'sequence':
                        sequence = (child.text or "").replace('\n', '').replace(' ', '')
                
                # Показываем первые 10 белков для теста структуры
                if count <= 10:
                    print(f"[{count}] UID: {uid} | Длина: {len(sequence)} aa")
                    print(f"Последовательность: {sequence[:40]}...")
                    print("-" * 30)
                elif count == 11:
                    print("...база успешно определена, чтение продолжается потоком...")
                    break # Остановим цикл на 11, чтобы просто проверить, что всё починилось
                
                elem.clear() # Мгновенно чистим RAM
                
    print("Проверка структуры завершена успешно!")
except FileNotFoundError:
    print(f"Ошибка: Файл не найден по пути {file_path}")
except Exception as e:
    print(f"Произошла ошибка при парсинге: {e}")
