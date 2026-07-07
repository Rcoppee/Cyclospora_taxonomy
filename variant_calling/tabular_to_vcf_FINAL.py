#!/usr/bin/env python3
import csv, sys, os, ast

if len(sys.argv) < 3:
    print("Usage: python tabular_to_vcf.py input.tsv output.vcf")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

def safe_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0

def parse_observations(obs_str):
    try:
        return ast.literal_eval(obs_str)
    except Exception:
        return {}

sample_name = os.path.splitext(os.path.basename(input_file))[0]
deletion_memory = {}

with open(input_file, newline='') as infile, open(output_file, 'w') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write("##source=tabular_to_vcf_refined\n")
    outfile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
    outfile.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant Type\">\n")
    outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    # Ajout du champ VF dans l'en-tête
    outfile.write("##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant allele frequency (%)\">\n")
    outfile.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

    for row in reader:
        chrom = row.get("Chromosome", ".")
        pos_str = row.get("Position", ".")
        try:
            pos = int(pos_str)
        except:
            continue
            
        ref_base = row.get("RefBase", ".")
        dp_str = row.get("Depth", "0")
        dp = int(dp_str) if dp_str.isdigit() else 0
        vtype = row.get("VariantType", "REF").upper()
        vf_col = safe_float(row.get("VariantFreq(%)", "0"))
        variant_seq = row.get("VariantSeq", "").strip()
        obs = parse_observations(row.get("Observations", "{}"))

        # Nettoyage mémoire délétions
        deletion_memory = {k: v for k, v in deletion_memory.items() if k[0] == chrom and k[1] >= pos}
        inherited_count = deletion_memory.get((chrom, pos), 0)

        # Valeurs VCF par défaut
        ref = ref_base
        alt = "."
        vf = vf_col
        gt = "./."

        # --- LOGIQUE PAR TYPE ---

        # 1. RÉFÉRENCE
        if vtype == "REF":
            alt = "."
            gt = "0/0"
            vf = vf_col

        # 2. INSERTIONS (Corrigé pour éviter la duplication de l'ancre)
        elif vtype == "INS":
            ins_seq = ""
            # On récupère la séquence brute
            if variant_seq and variant_seq not in (".", "*"):
                ins_seq = variant_seq
            else:
                ins_keys = [k.lstrip('+') for k in obs.keys() if isinstance(k, str) and k.startswith('+')]
                if ins_keys:
                    ins_seq = max(ins_keys, key=len)
            
            if ins_seq:
                # Si la séquence insérée commence déjà par la base de REF, 
                # on l'utilise telle quelle pour l'ALT.
                # Sinon, on ajoute la base de REF devant (standard VCF).
                if ins_seq.upper().startswith(ref.upper()):
                    alt = ins_seq.upper()
                else:
                    alt = ref.upper() + ins_seq.upper()
            else:
                alt = "<INS>"
            
            vf = vf_col if vf_col > 0 else 0.0
            gt = "1/1" if vf > 80 else "0/1"

        # 3. DÉLÉTIONS (avec mémoire pour les single-nucleotides)
        elif vtype in ("DEL", "DEL_DOWNSTREAM"):
            del_candidates = {}
            for k, v in obs.items():
                if isinstance(k, str) and k.startswith('-'):
                    del_candidates[k.lstrip('-')] = int(v)
                elif isinstance(k, str) and all(c.upper() in "ATGCN" for c in k):
                    del_candidates[k] = int(v)

            if del_candidates:
                del_count = max(del_candidates.values())
                vf = (del_count / dp * 100) if dp > 0 else vf_col
            elif inherited_count > 0:
                vf = (inherited_count / dp * 100) if dp > 0 else vf_col
            else:
                vf = vf_col

            alt = "*"
            gt = "1/1" if vf > 80 else "0/1"

        # 4. SNP / SUBSTITUTIONS
        elif vtype in ("SNP", "SUB"):
            alt = variant_seq if variant_seq not in (".", "", "*") else "."
            vf = vf_col
            gt = "1/1" if vf > 80 else "0/1"

        # --- MISE À JOUR MÉMOIRE (pour les DEL suivantes) ---
        for k, v in obs.items():
            if isinstance(k, str) and k.startswith('-'):
                del_seq = k.lstrip('-')
                for i in range(1, len(del_seq) + 1):
                    deletion_memory[(chrom, pos + i)] = int(v)

        # ÉCRITURE
        if dp == 0 and vtype != "REF": gt = "./."
        
        info = f"DP={dp};TYPE={vtype}"
        format_field = "GT:VF"
        sample_field = f"{gt}:{vf:.2f}"
        
        outfile.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\t{format_field}\t{sample_field}\n")