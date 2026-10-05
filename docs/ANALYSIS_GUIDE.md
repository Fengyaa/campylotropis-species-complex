# Analysis guide

Scripts use repository-relative paths unless arguments are required. Outputs go to `results/` or a requested new directory. These are independent analyses and command records, not a complete study pipeline.

## Structure and trees

`bash scripts/admixture.sh INPUT.vcf.gz INTERGENIC.pos OUTDIR` is cleaned from the supplied commands: MAF 0.01 throughout, consistent filenames, and one CV-enabled fit per K. The original ran two fits per K. The cleaned version is not described as the exact historical execution; stochastic new fits need not reproduce archived Q values. CV logs were not recovered, so no CV values are invented. `Rscript scripts/plot_admixture.R` uses the final K=2–8 matrices; K=1–20 matrices and the original selected PDF are retained.

Convert the final fourfold-degenerate, thin500 VCF with `convert_vcf_to_nexus.rb` from the [species-tree tutorial](https://github.com/millanek/tutorials/tree/master/bayesian_species_tree_inference). Run `python scripts/prepare_svd_input.py input.nex numbered.nex`, inspect the exported name/number map, then use `scripts/svdquartets.nex` in PAUP*. The helper checks 259 taxa and outgroup positions but assumes the supplied population ordering. Original conversion software and the full alignment remain external; mixed exploratory filenames have been removed from the cleaned commands.

## Demography

`bash scripts/psmc.sh SAMPLE_LIST REF.fa GVCF_DIR VCF2FQ.py PSMC_DIR OUTDIR` preserves the supplied GATK → vcf2fq → fq2psmcfa → PSMC sequence and plotting parameters. GATK 4.2.5.0, Python, Perl and [PSMC](https://github.com/lh3/psmc) are external. The original `vcf2fq.py` and per-individual gVCF sample list were not found locally; no substitute is invented. Samples run sequentially to avoid imposing an unspecified concurrency level on 28-GB GATK jobs. Supplied μ/g values do not independently establish the calibration of old figure exports.

`python scripts/summarize_fsc.py` summarizes six retained models. Their TPL/EST, likelihoods and AIC are in `data/demography/fsc`; the best-model parameter and 95% bootstrap summaries are retained. Identical observed spectra are stored once with original-name mappings. Full optimization is not required to inspect these materials.

## Genomic windows

`scripts/dsuite.sh` is a parameterized template for [Dsuite](https://github.com/millanek/Dsuite) using author-confirmed `-w 50,25`. It needs an external VCF, SETS, rooted tree and trios; names must agree and outgroups must be coded as required by Dsuite. The final species configurations are MY1–CYU–CPO and CY32–CPO–CYU; their separate SETS/trio files are bundled in `data/introgression/`. `species_dinvestigate.sh` runs these two local-window configurations without requiring a tree or a new Dtrios/Fbranch run.

The compact 10-kb table retains diversity, differentiation, divergence, XP-CLR and fd. Original full-genotype inputs and exact popgenWindows/XP-CLR commands remain external. Species-reference fd summaries have been recalculated from the two supplied 50/25 outputs; sympatric fd is unchanged. See [added analysis instructions](SPECIES_FD_UPDATE.md) for start-position assignment, filtering, configuration averaging, Table S9 and reference sensitivity.

`Rscript scripts/window_features.R` runs the retained threshold/coverage/NBS and GO/KEGG calculations, omitting only an unnecessary comparison to an older run. Dependencies: data.table, R.utils, clusterProfiler, ggplot2. Plot with `plot_window_features.R` or `plot_genomic_landscapes.R`. Putative NBS IDs are annotation outputs; source BLAST/UniProt inputs are external.

## Traits and GEA

`Rscript scripts/morphology_rda_pca.R` uses vegan for the final RDA/PCA and checks agreement with the shared six-variable climate table. It writes tables. The archive includes no alternative bio12 model.

For each chromosome dosage file, first create an output directory, then run:

```bash
Rscript scripts/genomic_prda.R chromosome.raw data/climate/all_indv_env.corrected.tsv --pcs 3 --axes 6 --id-column IID --out-prefix results/prda/chr1
```

The scan script retains 3.5 and 4.0 SD as options from the earlier script; only the final 3.5 SD/|r|>0.5 candidate set is packaged. `python scripts/gea_overlap.py` summarizes its overlap with BayPass. `Rscript scripts/gea_window_statistics.R` applies the previous spatial tests to this final set only.

`scripts/baypass.sh` is a BayPass auxiliary-model template, not an original execution log. Supply the full filtered SNP count file, climatic PC matrix, original Omega matrix and population count. It does not reconstruct unavailable MCMC settings or identify the covariance-estimation panel. The climate source is shared with RDA; the original BayPass PC matrix must retain its population order. Final candidates use BF>20 decibans and per-covariate upper-tail rank probability <0.001.

## GWAS, ancestry and candidate region

`scripts/gemma.sh` is a template for relatedness + genetic PC1–PC3, without climate. Supply matched BED/BIM/FAM, kinship, phenotype matrix and no-header covariates (intercept + PCs). The 210-row phenotype/FAM panel is retained; missing phenotypes are excluded per trait, leaving 195 for GL. The labeled phenotype TSV must be converted to GEMMA's numeric input format; it is not itself a no-header phenotype matrix. [GEMMA manual](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.tex).

`process_gemma.R` applies genome-wide corrections to all chromosome association files, and `plot_gemma.R` plots full test outputs. Both print argument usage. Compact significant/regional GL results are included; full scans require external input.

`Rscript scripts/hybrid_ancestry.R` uses selected-marker matrices and parental references; `plot_hybrid_ancestry.R` plots estimates. This two-lineage model does not establish categorical hybrid generations.

Run `python scripts/local_affinity.py`, then `python scripts/local_affinity_sensitivity.py`, then `Rscript scripts/plot_local_affinity.R`. The compact quality-masked chromosome 8/6 caches and reference order are retained, including original qualification flags. `plot_chr8_candidate.R` and `plot_MYB_domains.R` use retained regional/protein data. Full extraction VCFs remain external.

## Spatial relationships

Run `Rscript scripts/mantel.R` before `Rscript scripts/mmrr.R`. They use the retained diploid WC84 matrices, linearization FST/(1−FST), environmental/geographic distances, population subsets and 9,999 permutations. R dependencies: vegan and the included PopGenReport lgrMMRR with its original GPL notices.

## Checks

The curation checks syntax, sample/data consistency and short numerical summaries. Full upstream VCF scans and demographic fitting were not rerun. The species-reference update reused completed 2,000-replicate feature bootstraps and 9,999-shift GEA comparisons; portable window summaries and Table S9 were checked directly against those completed results. Software versions from earlier local checks are listed in `software_versions.tsv`; they are not a locked environment.
