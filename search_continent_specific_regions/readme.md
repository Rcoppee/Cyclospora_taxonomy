# Search for continent-discriminant genomic regions

This directory contains the scripts used to identify genomic regions showing the strongest differentiation between predefined geographic groups (e.g. America, Asia, Africa). The approach scans a multi-sample VCF file and reports genomic regions containing variants that maximize discrimination between the selected populations.

## Requirements

- Python 3 (>=3.8)
- Standard Python libraries (no additional dependencies required unless specified in the script)

## Input files

Two input files are required:

### 1. Multi-sample VCF

A merged, uncompressed VCF file generated from the variant calling pipeline.

Example:

```text
data.vcf
```

### 2. Group assignment file

A tab-delimited text file containing one sample per line together with its geographic group.

Example:

```text
Sample_001    America
Sample_002    America
Sample_003    Asia
Sample_004    Africa
...
```

Only samples listed in this file will be included in the analysis.

## Running the analysis

Execute the script using:

```bash
python3 find_discriminant_regions.py data.vcf groups.txt results.tsv
```

where:

- `data.vcf` is the merged VCF file;
- `groups.txt` contains the population assignment of each sample;
- `results.tsv` is the output file listing the candidate discriminant regions.

## Output

The output is a tab-separated (`.tsv`) file containing the genomic regions showing the highest discriminatory power between the predefined geographic groups.

For each candidate region, the output includes genomic coordinates together with summary statistics allowing prioritization of markers for downstream molecular assay development (Melting qPCR or Sanger sequencing).

## Notes

- The script is designed to compare any user-defined populations and is therefore not restricted to continental analyses.
- Group labels are provided by the user and may represent continents, countries, lineages, or any other predefined sample classification.
- The resulting candidate regions should be manually inspected before primer or probe design.
