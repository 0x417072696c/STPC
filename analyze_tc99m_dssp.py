import allymind_core, cascade, pdb_parser
# Встроенная последовательность HSA
hsa_seq = (
    'DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL'
)

print("Скачиваю PDB 1AO6...")
pdb_file = pdb_parser.download_pdb('1AO6')
if pdb_file:
    print("Извлекаю вторичную структуру через DSSP...")
    ss_map = pdb_parser.extract_aligned_ss_map(hsa_seq, pdb_file)
    if ss_map:
        site = 'KVPQVSTPTLVEVSR'
        pos = hsa_seq.find(site)
        print("Структура сайта Tc-99m (DSSP):")
        for i, aa in enumerate(site):
            res_pos = pos + i
            ss = ss_map[res_pos] if res_pos < len(ss_map) else '?'
            print(f"{res_pos+1}:{aa} -> {ss}")
        print()
        print("=== Tc-99m site scan WITH DSSP ===")
        print(f'{"Pos.":<6} {"Orig":<4} {"Depth(HB)":<10} {"ddG(HB)":<10} {"Structure":<12} {"J_thr":<10}')
        for i, aa in enumerate(site):
            res_pos = pos + i
            mut = f'{aa}{res_pos+1}A'
            report = allymind_core.analyze_mutation(hsa_seq, mut, ss_map=ss_map)
            depth = report.get('depth', 0)
            ddG = report.get('delta_delta_G', 0)
            ss = report.get('secondary_structure', '?')
            thr = report.get('therapeutic_threshold', 0)
            print(f'{res_pos+1:<6} {aa:<4} {depth:<10.4f} {ddG:<10.4f} {ss:<12} {thr:<10.4f}')
    else:
        print("DSSP не сработал. Проверь, что mkdssp.exe лежит в текущей папке.")
else:
    print("Не удалось скачать PDB 1AO6.")
