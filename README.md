# Campylotropis species complex: analysis code and supporting data

Analysis scripts and supporting data for **Genomic, ecological, and morphological insights into speciation with introgression in the *Campylotropis macrocarpa* complex**.

[中文说明](README_zh.md) · [Analysis guide](docs/ANALYSIS_GUIDE.md) · [Methods](docs/METHODS.md) · [Data and Code Availability](docs/DATA_AND_CODE_AVAILABILITY.md)

The repository contains scripts, processed inputs and selected results for population genomics, demographic inference, morphological associations and genotype–environment analyses. Analyses can be run individually; workflows starting from raw sequencing reads require external data and software.

## Repository contents

| Analysis | Main scripts | Inputs and results |
|---|---|---|
| Population structure | `admixture.sh`, `plot_admixture.R` | Q matrices (K = 1–20), sample order and ancestry plot |
| Phylogeny | `prepare_svd_input.py`, `svdquartets.nex` | Population partitions and trees |
| Demography | `psmc.sh`, `summarize_fsc.py` | PSMC outputs, six fastsimcoal2 models, spectra and summaries |
| Genomic differentiation and allele sharing | `species_dinvestigate.sh`, `prepare_species_fd.py`, `window_features.R`, `reference_features.R` | Dsuite outputs, 10-kb statistics, annotations and reference-sensitivity analyses |
| ADP-binding/NBS overlap | `table_s9.py` | Gene counts, enrichment statistics and Table S9 |
| Morphology and climate | `morphology_rda_pca.R` | Eight traits in 192 individuals and climatic data |
| Genotype–environment associations | `genomic_prda.R`, `baypass.sh`, `gea_overlap.py`, `gea_window_statistics.R` | pRDA/BayPass candidates, shared windows and genomic-statistic comparisons |
| GWAS and candidate region | `gemma.sh`, `process_gemma.R`, local-affinity and plotting scripts | Phenotypes, glandular-hair associations, regional genotypes and MYB annotation |
| Hybrid ancestry | `hybrid_ancestry.R`, `plot_hybrid_ancestry.R` | Marker matrices, reference panels and ancestry estimates |
| Spatial relationships | `mantel.R`, `mmrr.R` | Geographic, environmental and genetic distances |

- `scripts/`: analysis and plotting scripts.
- `data/`: processed inputs, annotations, sample metadata and model files.
- `results/`: numerical summaries and selected outputs.
- `docs/`: methods, dependencies and usage instructions.
- `vendor/`: third-party code with attribution and license texts.
- `CONTENTS.tsv`: file descriptions, sizes and SHA-256 checksums.

## Analysis settings

- Morphological RDA and genomic pRDA use **bio3, bio7, bio10, bio13, bio15 and srad_07** from the same climatic dataset.
- Morphological RDA includes **192 individuals from 28 populations**: R² = 34.76%, adjusted R² = 32.64%; RDA1 and RDA2 account for 20.98% and 10.45% of total variation.
- Genomic pRDA conditions on genetic PC1–PC3. Candidate SNPs exceed 3.5 loading standard deviations and have an absolute genotype–environment correlation >0.5.
- GEMMA includes genomic relatedness and genetic PC1–PC3 as covariates.
- BayPass tests climatic associations across the filtered genome-wide SNP dataset.
- Dsuite Dinvestigate uses **50 usable SNPs per window and a 25-SNP step**. Species-level configurations are **MY1–CYU–CPO** and **CY32–CPO–CYU**. Available configuration means are averaged within the 10-kb reporting grid.
- ADMIXTURE uses intergenic SNPs filtered at MAF 0.01, call rate 0.90 and minimum spacing 1 kb. The included ancestry plot is `maf01ms09_re.pdf`.

See [Methods](docs/METHODS.md) for parameter definitions and [Introgression analyses](docs/INTROGRESSION.md) for window aggregation, reference sensitivity and figure scripts.

## Getting started

Run the following from the repository root:

```bash
Rscript scripts/morphology_rda_pca.R
python scripts/gea_overlap.py
python scripts/summarize_fsc.py
```

These analyses use bundled inputs. Other scripts require additional inputs or longer computation; see the [Analysis guide](docs/ANALYSIS_GUIDE.md). R dependencies include vegan, data.table, R.utils, ggplot2, tidyr, readr, dplyr and clusterProfiler. Python dependencies include numpy, pandas and pypdf. Software versions are listed in [software_versions.tsv](docs/software_versions.tsv).

## Data access

- Reference genome and associated sequencing data: **CNGBdb CNP0005557**.
- Whole-genome resequencing data: **CNGBdb CNP0007152**.
- Code and processed data: [Fengyaa/campylotropis-species-complex](https://github.com/Fengyaa/campylotropis-species-complex).

External software and databases retain their respective licenses and citation requirements; see [Third-party notices](THIRD_PARTY_NOTICES.md).
