"""
expand_screening_clean.py – массовый скрининг с фильтрацией генов.
Пропускает гены с ;, пробелами, длиной > 50 символов, без PDB.
"""
import allymind_core, cascade, pdb_parser
import numpy as np, sqlite3, os, csv, time
import matplotlib.pyplot as plt

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")
WINDOW = 5
J_NAT = cascade.J_NAT

# 1. Получаем список генов из ClinVar с фильтрацией
print("Loading filtered gene list from ClinVar...")
conn = sqlite3.connect(allymind_core.DB_CLINVAR)
all_genes = [row[0] for row in conn.execute(
    "SELECT DISTINCT gene FROM mutations WHERE gene != '' AND gene IS NOT NULL"
).fetchall()]
conn.close()

# Фильтруем гены
clean_genes = []
for g in all_genes:
    if ';' in g or ' ' in g or len(g) > 50:
        continue
    clean_genes.append(g)
print(f"Total unique genes: {len(all_genes)} → filtered: {len(clean_genes)}")

results = []
start_time = time.time()

for i, gene in enumerate(clean_genes):
    print(f"\rProcessing {i+1}/{len(clean_genes)}: {gene:<20}", end="", flush=True)

    # Загружаем последовательность
    seq, uniprot = allymind_core.load_sequence(gene)
    if not seq:
        continue

    # Ищем PDB-структуру (только если есть, иначе пропускаем)
    pdb_file = pdb_parser.download_pdb(gene.lower())
    if not pdb_file:
        continue

    # Извлекаем ss_map и кэшируем SASA/расстояние до центра
    ss_map = None
    sasa_dict = {}
    dist_dict = {}
    try:
        ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file, chain='A')
        sasa_dict = cascade.compute_sasa(pdb_file)
        dist_dict = cascade.compute_distance_to_center(pdb_file)
    except:
        pass

    if ss_map is None:
        continue

    # Находим зону влияния домена
    def find_core_domain(ss_map):
        best_start, best_end = 0, 0
        current_start = None
        for i, s in enumerate(ss_map):
            if s in ('H', 'E'):
                if current_start is None:
                    current_start = i
            else:
                if current_start is not None:
                    if i - current_start > best_end - best_start:
                        best_start, best_end = current_start, i
                    current_start = None
        if current_start is not None:
            if len(ss_map) - current_start > best_end - best_start:
                best_start, best_end = current_start, len(ss_map)
        if best_end - best_start == 0:
            L = len(ss_map)
            best_start, best_end = int(L * 0.33), int(L * 0.66)
        return max(0, best_start - WINDOW), min(len(ss_map), best_end + WINDOW)

    zone_start, zone_end = find_core_domain(ss_map)
    depths = []

    for pos in range(zone_start + 1, zone_end + 1):
        wt = seq[pos - 1]
        for alt in AMINO_ACIDS:
            if alt == wt:
                continue
            report = allymind_core.analyze_mutation(seq, f"{wt}{pos}{alt}", ss_map=ss_map)
            depth = report.get('depth', 0)
            ddG = cascade.cubic_ddG(depth)

            # Проверка открытости (петля + SASA + расстояние до центра)
            exposed = False
            if ss_map[pos - 1] == 'C':
                if sasa_dict and dist_dict:
                    sasa = sasa_dict.get(pos, 0.0)
                    dist = dist_dict.get(pos, 0.0)
                    exposed = (sasa > 70.0) and (dist > 12.0)
                else:
                    exposed = True

            if exposed and depth > 0.8 and ddG > 0.1:
                depths.append(depth)

    if depths:
        mean_depth = np.mean(depths)
        kappa = mean_depth / J_NAT
        results.append((gene, uniprot, len(seq), kappa, len(depths)))
    else:
        results.append((gene, uniprot, len(seq), 0.0, 0))

# 2. Сохраняем в CSV
with open("kappa_clean_proteins.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["gene", "uniprot", "length", "kappa", "num_traps"])
    for row in results:
        writer.writerow(row)

elapsed = time.time() - start_time
print(f"\nProcessed {len(results)} proteins in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes).")

# 3. Строим график
if results:
    lengths = [r[2] for r in results]
    kappas = [r[3] for r in results]
    genes = [r[0] for r in results]

    fig, ax = plt.subplots(figsize=(14, 8))
    scatter = ax.scatter(lengths, kappas, c='#3daee9', alpha=0.6, s=30)
    ax.axhline(y=16, color='gold', linestyle='--', linewidth=2, alpha=0.8, label='κ = 16 (8+5+3)')
    ax.set_xlabel('Protein Length (residues)')
    ax.set_ylabel('κ = mean depth / J_NAT')
    ax.set_title(f'Pan-Protein Vulnerability Map: {len(results)} proteins')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('kappa_clean_proteins.png', dpi=150)
    plt.show()
    print("Graph saved to kappa_clean_proteins.png")