#!/bin/bash

# ------------------------------------------------------------
# Script pour convertir tous les fichiers .tsv en .vcf
# en utilisant le script Python tabular_to_vcf.py
# ------------------------------------------------------------

# ==== PARAMÈTRES À MODIFIER ====
REF="ccayetanensis.fasta"              # ton fichier de référence FASTA
SCRIPT="tabular_to_vcf_FINAL.py"         # ton script Python de conversion
INPUT_DIR="."                      # dossier contenant les fichiers .tsv
OUTPUT_DIR="vcf"                   # dossier de sortie
# ================================

# Vérifie que le fichier de référence existe
if [ ! -f "$REF" ]; then
    echo "❌ Erreur : fichier de référence introuvable ($REF)"
    exit 1
fi

# Vérifie que le script Python existe
if [ ! -f "$SCRIPT" ]; then
    echo "❌ Erreur : script Python introuvable ($SCRIPT)"
    exit 1
fi

# Crée le dossier de sortie s'il n'existe pas
mkdir -p "$OUTPUT_DIR"

# Boucle sur tous les fichiers TSV du dossier
for file in "$INPUT_DIR"/*.tsv; do
    [ -e "$file" ] || continue  # saute si aucun fichier trouvé

    base=$(basename "$file" .tsv)
    output="${OUTPUT_DIR}/${base}.vcf"

    echo "➡️  Conversion de $file → $output"
    python3 "$SCRIPT" "$file" "$output" "$REF"

    if [ $? -eq 0 ]; then
        echo "✅ Fichier généré : $output"
    else
        echo "⚠️  Erreur lors du traitement de $file"
    fi
done

echo "🎉 Conversion terminée. Les fichiers VCF sont dans $OUTPUT_DIR/"
