#!/usr/bin/env python3
import re
import sys
from collections import Counter

def parse_pileup_variants_stable(pileup_file, output_file, min_depth=5, indel_freq_threshold=0.5, sep='\t'):
    """
    Parse pileup robustly and output variants:
    - Treat insertions/deletions correctly per position
    - Downstream positions of dominant deletions are marked DEL_downstream
    - Filter indels <50% frequency
    - Maintain all positions, no skips
    """
    indel_re = re.compile(r'([+\-])(\d+)([ACGTNacgtn]+)')
    deletion_map = {}  # (chrom,pos) -> '*' for downstream positions

    with open(pileup_file, 'r') as inf, open(output_file, 'w') as outf:
        outf.write(sep.join([
            "Chromosome","Position","RefBase","Depth","VariantType",
            "VariantSeq","VariantFreq(%)","Observations","AnchorType"
        ]) + "\n")

        for line in inf:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[2].upper()
            try:
                depth = int(fields[3])
            except ValueError:
                continue
            bases = fields[4]

            key = (chrom, pos)

            # 1) Vérifier si cette position est planifiée DEL downstream
            if key in deletion_map:
                outf.write(sep.join([chrom, str(pos), ref, str(depth),
                                     "DEL", "*", "", str(deletion_map[key]), "DEL_downstream"]) + "\n")
                continue

            # 2) Vérifier profondeur minimale
            if depth < min_depth:
                outf.write(sep.join([chrom, str(pos), ref, str(depth),
                                     "LOW_DEPTH", "N", "", "", "LOW_DEPTH"]) + "\n")
                continue

            # 3) Nettoyer la chaîne de bases
            bases_clean = re.sub(r'\^.', '', bases).replace('$', '')

            counts = Counter()
            insertion_counts = Counter()
            deletion_counts = Counter()

            # 4) Parcours des bases
            i = 0
            L = len(bases_clean)
            while i < L:
                ch = bases_clean[i]

                if ch in "+-":
                    m = indel_re.match(bases_clean[i:])
                    if m:
                        sign, num_s, seq = m.groups()
                        num = int(num_s)
                        seq = seq.upper()
                        if sign == '+':
                            insertion_counts[seq] += 1
                            counts[f"+{seq}"] += 1
                        else:
                            deletion_counts[seq] += 1
                            counts[f"-{seq}"] += 1
                        i += 1 + len(num_s) + num
                        continue
                    else:
                        i += 1
                        continue
                elif ch in ".,":  # base référence
                    counts[ref] += 1
                    i += 1
                    continue
                uc = ch.upper()
                if uc in ('A','C','G','T','N'):
                    counts[uc] += 1
                    i += 1
                    continue
                i += 1

            # 5) Construire Observations
            obs_dict = {k:v for k,v in counts.items() if v>0}

            variant_type = "REF"
            var_seq = ref
            variant_freq = 0.0
            anchor_type = "REF"

            # a) Vérifier insertions
            dominant_ins = None
            if insertion_counts:
                ins_seq, ins_count = insertion_counts.most_common(1)[0]
                if ins_count / depth >= indel_freq_threshold:
                    dominant_ins = ins_seq
                    variant_type = "INS"
                    var_seq = ref + ins_seq
                    variant_freq = (ins_count / depth) * 100
                    anchor_type = "INS_anchor"

            # b) Vérifier délétions
            dominant_del = None
            if deletion_counts:
                del_seq, del_count = deletion_counts.most_common(1)[0]
                if del_count / depth >= indel_freq_threshold:
                    dominant_del = del_seq
                    anchor_type = "DEL_anchor"
                    # planifier positions downstream
                    for j in range(1, len(del_seq)+1):
                        deletion_map[(chrom, pos + j)] = str({del_seq: del_count})
                    # La position actuelle reste REF
                else:
                    dominant_del = None

            # c) Si pas d’insertion dominante, choisir base dominante
            if variant_type != "INS":
                base_counts = {k:v for k,v in counts.items() if not k.startswith(('+','-'))}
                if base_counts:
                    top_base, top_count = max(base_counts.items(), key=lambda x:x[1])
                    var_seq = top_base
                    variant_freq = (top_count / depth) * 100
                    if top_base == ref:
                        variant_type = "REF"
                        anchor_type = "REF"
                    else:
                        variant_type = "SNP"
                        anchor_type = "SNP"

            freq_s = f"{variant_freq:.2f}" if variant_type in ["SNP","REF","INS"] else ""
            outf.write(sep.join([chrom, str(pos), ref, str(depth),
                                 variant_type, var_seq, freq_s, str(obs_dict), anchor_type]) + "\n")

    print(f"[OK] Fichier écrit : {output_file}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: parse_pileup_variants_stable.py input.pileup output.tsv [min_depth]")
        sys.exit(1)
    inp = sys.argv[1]
    outp = sys.argv[2]
    md = int(sys.argv[3]) if len(sys.argv) >= 4 else 5
    parse_pileup_variants_stable(inp, outp, min_depth=md)
