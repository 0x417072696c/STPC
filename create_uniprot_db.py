import gzip, sqlite3, xml.etree.ElementTree as ET, sys

file_path = r"C:\Users\Colorful X15 AT23\Desktop\СТПК\Levi\AllyMind 10_2\uniprot_sprot.xml.gz"
db_path = "uniprot.db"
print("Создание базы данных Swiss-Prot...")
conn = sqlite3.connect(db_path)
conn.execute("CREATE TABLE IF NOT EXISTS proteins (gene TEXT, uniprot_id TEXT, sequence TEXT)")
conn.execute("DELETE FROM proteins")
conn.commit()

inserted = 0
with gzip.open(file_path, 'rb') as f:
    context = ET.iterparse(f, events=('end',))
    for event, elem in context:
        tag = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
        if tag == 'entry':
            uid = None
            seq = None
            gene = None
            for child in elem:
                child_tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                if child_tag == 'accession' and uid is None:
                    uid = child.text
                elif child_tag == 'sequence':
                    seq = (child.text or "").replace('\n', '').replace(' ', '')
                elif child_tag == 'gene':
                    for gchild in child:
                        gtag = gchild.tag.split('}')[-1] if '}' in gchild.tag else gchild.tag
                        if gtag == 'name' and gchild.attrib.get('type') == 'primary':
                            gene = gchild.text
                            break
            if uid and seq:
                conn.execute("INSERT INTO proteins (gene, uniprot_id, sequence) VALUES (?,?,?)",
                             (gene or uid, uid, seq))
                inserted += 1
                if inserted % 10000 == 0:
                    print(f"Обработано: {inserted}")
                    conn.commit()
            elem.clear()

conn.commit()
conn.close()
print(f"Готово! Записей: {inserted}")
