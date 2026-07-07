#! /bin/bash

#Script Romain Coppee
#Creation data: 09/02/2020
#Last modification: 09/06/2026

# Spécifiez le chemin du répertoire contenant vos fichiers
repertoire="$PWD"

# Assurez-vous que le répertoire existe
if [ ! -d "$repertoire" ]; then
    echo "Le répertoire spécifié n'existe pas. Veuillez spécifier le chemin correct."
    exit 1
fi

# Parcourez les fichiers correspondants au motif dans le répertoire
for fichier in "$repertoire"/*_R1.fastq.gz; do
    # Vérifiez si le fichier existe avant de procéder
    if [ -e "$fichier" ]; then
        # Obtenez le nom de base du fichier (sans le suffixe _1.fastq.gz)
        current_name=$(basename "$fichier" "_R1.fastq.gz")
        
        echo "$current_name"
        bwa mem -t 12 ccayetanensis.fasta "${current_name}_R1.fastq.gz" "${current_name}_R2.fastq.gz" > "${current_name}.sam"
        samtools view -@ 12 -b -S "${current_name}.sam" > "${current_name}.bam"
        rm "${current_name}.sam"
        samtools sort -@ 12 "${current_name}.bam" -o "${current_name}.sorted.bam"
        rm "${current_name}.bam"
        samtools index "${current_name}.sorted.bam"
    fi
done

