#!/bin/bash

# Arrête le script si une commande échoue
set -e

for f in *_cleaned.pileup; do
    base=$(basename "$f" _cleaned.pileup)

    echo "------------------------------------------"
    echo "🔹 Processing file: $f"
    echo "------------------------------------------"

    python pileup_to_diversity2.py "$f" "${base}.tsv"

    echo "✅ Created: ${base}.tsv"
    echo
done
