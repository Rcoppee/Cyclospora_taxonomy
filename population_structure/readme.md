# Population structure analyses

This directory contains the R scripts used to perform the genome-wide population genomic analyses described in the manuscript.

The scripts reproduce all analyses used to investigate the genetic structure of *Cyclospora cayetanensis*, including phylogenetic reconstruction, multivariate analyses, and population assignment analyses.

## Analyses performed

The workflow includes:

- Neighbor-Joining (NJ) phylogenetic reconstruction;
- Principal Component Analysis (PCA);
- Discriminant Analysis of Principal Components (DAPC);
- Population assignment (posterior membership probabilities).
  
The scripts also generate the figures presented in the manuscript.

---

## Running the analyses

Open the corresponding script in **RStudio** (or any R environment) and execute it.

```R
source("population_structure.R")
```

The scripts are fully annotated and can be run section by section.

---

## Choice of the analysis

The same workflow was used to investigate population structure at different hierarchical levels:

- CYCLONE lineages;
- Continents;
- Countries.

These analyses are controlled within the script by commenting or uncommenting the corresponding code sections.

To select the desired analysis, simply add or remove the `#` symbols as indicated in the script comments.

Only one analysis should be activated at a time.

---

## Input data

The scripts require:

- the filtered multi-sample VCF file generated in the `variant_calling` workflow;
- metadata files describing sample group assignments (lineages, continents or countries).

Input file paths can be modified directly within the script.

---

## Output

The workflow generates:

- Neighbor-Joining phylogenetic trees;
- PCA plots;
- DAPC plots;
- Posterior membership probability plots.

These outputs correspond to the analyses presented in the manuscript.

---

## Required R packages

The scripts require several R packages, including:

- **vcfR**
- **adegenet**
- **ape**
- **poppr**
- **igraph**
- **ggplot2**
- **dplyr**
- **tidyr**
