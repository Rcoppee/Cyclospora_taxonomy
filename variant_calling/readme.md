# Variant Calling

This directory contains the scripts used to identify genome-wide genetic variants (SNPs and indels) from mapped *Cyclospora cayetanensis* genomes and to generate a merged VCF file for downstream population genomic analyses.

---

## Workflow

The variant calling workflow consists of three successive steps:

1. Generation and cleaning of pileup files
2. Variant calling and VCF generation
3. Compression and merging of individual VCF files

The scripts must be executed sequentially.

---

## Files

### `1_launch_pileup.sh`

This script generates a pileup file for each mapped BAM file (`*_sorted.bam`) using **Samtools**.

After pileup generation, each file is automatically processed using the script:

```text
tools/clean.sh
```

which formats the pileup for downstream analyses.

### Input

- `*_sorted.bam`

### Output

- `*_cleaned.pileup`

### Usage

Run the script from the directory containing the BAM files:

```bash
bash 1_launch_pileup.sh
```

---

### `2_variant_calling.sh`

This script identifies genome-wide SNPs and indels from the cleaned pileup files.

Two auxiliary scripts located in the `tools/` directory are automatically called:

```text
tools/diversity_launch.sh
```

Identifies genetic variants (SNPs and indels).

```text
tools/convert_file.sh
```

Converts the detected variants into standard VCF format.

### Input

- `*_cleaned.pileup`

### Output

- Individual VCF files for each sample.

### Usage

```bash
bash 2_variant_calling.sh
```

---

### `3_compress_merge.sh`

This script compresses all individual VCF files and merges them into a single multi-sample VCF using **BCFtools**.

### Usage

```bash
bash 3_compress_merge.sh
```

---

## Variant filtering

Following variant calling, additional filtering was performed using **VCFtools** and **BCFtools** to generate the final dataset used for population genomic analyses.

The following criteria were applied:

- only polymorphic positions (at least one sample carrying an alternative allele) were retained;
- variants were required to have at least **75% interpretable genotypes** across samples;
- all positions containing at least one heterozygous genotype were removed to simplify downstream population structure analyses.

The resulting filtered VCF file was used for all phylogenetic and population genomic analyses presented in the manuscript.

---

## Software requirements

The workflow requires:

- Samtools
- BCFtools
- VCFtools

The auxiliary scripts located in the `tools/` directory must remain in their original location for the workflow to execute correctly.
