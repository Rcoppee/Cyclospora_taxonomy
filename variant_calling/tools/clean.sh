#!/bin/bash
# Usage: ./clean.sh input.pileup output.pileup

INPUT=$1
OUTPUT=$2

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: $0 input.pileup output.pileup"
    exit 1
fi

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "$INPUT" > "$OUTPUT"

echo "[OK] Fichier filtré : $OUTPUT"
