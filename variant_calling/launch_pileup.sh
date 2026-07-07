#!/bin/bash

# Arrête le script si une commande échoue
set -e

# Boucle sur tous les fichiers .sorted.bam du dossier courant
for f in *.sorted.bam; do
    # Récupère le nom de base sans extension
    base=$(basename "$f" .sorted.bam)

    echo "------------------------------------------"
    echo "🔹 Processing file: $f"
    echo "------------------------------------------"

    # Étape 1 : Génération du pileup avec samtools
    echo "➡️  Creating pileup: ${base}.pileup"
    for ctg in $(cut -f1 ccayetanensis.fasta.fai); do
    	samtools mpileup -A -B -x -a -f ccayetanensis.fasta -r ${ctg} "$f" >> "${base}.pileup"
    done

    # Étape 2 : Nettoyage avec ton script clean.sh
    echo "➡️  Cleaning pileup: ${base}.pileup -> ${base}_cleaned.pileup"
    ./clean.sh "${base}.pileup" "${base}_cleaned.pileup"

    # (Optionnel) Supprimer le pileup brut si tu veux économiser de la place
    # rm "${base}.pileup"

    echo "✅ Done with ${base}.sorted.bam"
    echo
done

