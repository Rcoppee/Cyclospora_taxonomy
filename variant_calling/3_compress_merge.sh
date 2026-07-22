#!/bin/bash

# --- Vérification basique ---
if [ "$(ls *.vcf 2>/dev/null | wc -l)" -eq 0 ]; then
    echo "❌ Aucun fichier .vcf trouvé dans le dossier courant."
    exit 1
fi

echo "🔹 Étape 1 : compression et indexation de chaque VCF"
for file in *.vcf; do
    base=$(basename "$file")
    echo "  → Traitement de $base"
    bgzip -c "$file" > "${file}.gz"
    tabix -p vcf "${file}.gz"
done

echo "🔹 Étape 2 : fusion avec bcftools"
bcftools merge *.vcf.gz -Oz -o merged.vcf.gz

echo "🔹 Étape 3 : indexation du fichier fusionné"
tabix -p vcf merged.vcf.gz

echo "✅ Terminé ! Fichier final : merged.vcf.gz"
