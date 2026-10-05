# Campylotropis species complex: analysis code and supporting data

Selected scripts and compact data supporting the manuscript on differentiation, introgression, morphology and climate in the *Campylotropis macrocarpa* species complex. [中文说明](README_zh.md)

This is an analysis archive, not an end-to-end reproduction pipeline. It includes the necessary analysis scripts, compact inputs and selected retained results. Large genotypes and upstream software are external. There is no universal runner, CI workflow or generated work directory.

## Contents

| Analysis | Main scripts | Data and retained results |
|---|---|---|
| Population structure | `admixture.sh`, `plot_admixture.R` | Final maf01ms90thin1k Q matrices (K=1–20), sample order and maf01ms09_re.pdf |
| Phylogeny | `prepare_svd_input.py`, `svdquartets.nex` | Population partitions and retained trees |
| Demography | `psmc.sh`, `summarize_fsc.py` | PSMC outputs, six FSC models, spectra and summaries |
| Genomic windows | `dsuite.sh`, `species_dinvestigate.sh`, `prepare_species_fd.py`, `window_features.R`, `reference_features.R`, plotting scripts | MY1/CY32 reference inputs, updated 10-kb statistics, annotations, NBS/GO results and reference sensitivity |
| ADP-binding/NBS overlap | `table_s9.py` | Unique gene counts, GO FDR and updated Table S9 |
| Morphology/climate | `morphology_rda_pca.R` | 192 complete individuals and shared climate tables |
| GEA | `genomic_prda.R`, `baypass.sh`, `gea_overlap.py`, `gea_window_statistics.R` | Final pRDA/BayPass candidates, shared windows and comparison summaries |
| GWAS/candidate region | `gemma.sh`, `process_gemma.R`, candidate/local-affinity scripts | Phenotypes, GL variants, local genotype matrices and MYB annotation |
| Hybrid ancestry | `hybrid_ancestry.R`, `plot_hybrid_ancestry.R` | Marker matrices, reference information and estimates |
| Spatial relationships | `mantel.R`, `mmrr.R` | Geographic/environmental/genetic matrices and summaries |

All scripts are in `scripts/`, inputs in `data/`, and selected/generated outputs in `results/`. `CONTENTS.tsv` lists file sources and SHA-256 checksums. See [analysis guide](docs/ANALYSIS_GUIDE.md) for dependencies and commands. The [species-reference update](docs/SPECIES_FD_UPDATE.md) documents the added methods, Table S9 and updated Fig. 3, Fig. S7 and Fig. 5d.

## Final settings

- Morphology RDA and genomic pRDA use the same corrected climate dataset: **bio3, bio7, bio10, bio13, bio15 and srad_07**. The final author clarification retained **bio13**.
- Morphology RDA: **192 individuals / 28 populations**; R²=0.3475757, adjusted R²=0.3264160; RDA1/RDA2=20.98108%/10.44732% of total variation.
- Genomic pRDA: genetic PC1–PC3 conditioning; final candidates at 3.5 SD and |r|>0.5.
- GEMMA: relatedness plus genetic PC covariates; **no environmental fixed effects**.
- Final BayPass association scan: **all filtered genome-wide SNPs**, not the intergenic structure panel.
- Dinvestigate: **50 usable SNPs, step 25**. Final species-reference configurations are **MY1–CYU–CPO** and **CY32–CPO–CYU**; the 10-kb descriptive summary averages available configuration means. Component estimates and common-support sensitivity are retained.
- ADMIXTURE: only the version underlying **maf01ms09_re.pdf**; alternative MAF panels are omitted.

## Example use

```bash
Rscript scripts/morphology_rda_pca.R
python scripts/gea_overlap.py
python scripts/summarize_fsc.py
```

These short analyses use bundled data. Other scripts may need external inputs or longer computations. R dependencies include vegan, data.table, R.utils, ggplot2, tidyr, readr, dplyr and clusterProfiler; Python analyses use numpy, pandas and pypdf. Upstream tools are named in each command script. Recorded prior software versions are in `docs/software_versions.tsv`, not a locked environment.

## Data and citation

Reference genome and associated sequencing: **CNGBdb CNP0005557**. Whole-genome resequencing: **CNGBdb CNP0007152**. See [Data and Code Availability](docs/DATA_AND_CODE_AVAILABILITY.md).

This local archive has not been uploaded. Add the real repository URL, manuscript title/authors, code/data license and release DOI when available. Third-party software retains its original terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). [Method notes](docs/METHODS.md) distinguish author confirmations from historical file provenance.
