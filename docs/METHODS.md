# Analysis methods and parameters

| Analysis | Description |
|---|---|
| Climatic predictors | Morphological RDA and genomic pRDA use bio3, bio7, bio10, bio13, bio15 and July solar radiation (srad_07). |
| Climate sources | 19 WorldClim v1.4 bioclimatic variables and 36 WorldClim v2.1 monthly solar radiation, wind speed and water vapor pressure variables, at 2.5-arc-minute resolution. |
| Morphology | Eight traits in 192 individuals with complete measurements from 28 populations, drawn from the 195-individual glandular-hair GWAS cohort. R² and adjusted R² are reported separately. |
| Genomic pRDA | 257 individuals; genetic PC1–PC3 conditioning; six constrained axes; candidate loading threshold 3.5 standard deviations and absolute marginal genotype–environment correlation >0.5. |
| BayPass | Climatic association testing across filtered genome-wide SNPs, using a population covariance matrix estimated from 204,672 intergenic SNPs thinned to a minimum spacing of 1 kb. Candidates require BF >20 decibans and per-covariate upper-tail rank probability <0.001. |
| GEMMA | Genomic relatedness and genetic PC1–PC3 covariates. The glandular-hair analysis includes 195 individuals and 3,469,695 SNP tests. |
| Dinvestigate | MY1–CYU–CPO and CY32–CPO–CYU configurations; 50 usable SNPs per window and a 25-SNP step; reference individuals excluded from their respective pooled groups. Native estimates are assigned by start coordinate to 1-based inclusive 10-kb intervals; available configuration means are averaged equally. |
| ADMIXTURE | Intergenic SNPs with MAF ≥0.01, call rate ≥0.90 and minimum spacing 1 kb. Q matrices are named `257.snps.hardfiltered.maf01ms90thin1k.K.Q`; the ancestry plot is `maf01ms09_re.pdf`. |
| SVDquartets | Fourfold-degenerate SNPs thinned to a minimum spacing of 500 bp; 259 taxa; outgroups at positions 131 and 244; 100 bootstrap replicates. The `SPECIES` partition defines populations. |
| PSMC | GATK all-sites calls → vcf2fq.py → fq2psmcfa -q20; PSMC parameters `-N25 -t15 -r5 -p '4+25*2+4+6'`. Plotting uses μ = 7 × 10⁻⁹ per site per generation and a generation time of 2 years. |
| Gene overlap | At least 1 bp overlap with the full annotated gene span, including exons and introns; each gene counted once per GO candidate set. Annotation intervals are merged before calculating coverage. |
| Window comparisons | Upper-tail cutoffs of 10%, 5%, 2.5% and 1%; genomic-statistic windows require ≥100 analyzed sites and finite values. Local fd summaries retain 0 < fd < 1. GEA comparisons use a 22-statistic BH correction family, 9,999 spatial shifts and 1,000 block-bootstrap replicates. |

## Genomic annotations and enrichment

Gene-body and repeat coverage are measured as the fraction of each window covered by the respective annotation. NBS density is expressed as the number of distinct overlapping putative NBS genes per megabase. Candidate-to-background ratios are evaluated with chromosome-stratified block bootstrapping.

GO enrichment uses one-sided hypergeometric tests. The reference universe comprises GO-annotated genes overlapping any eligible window in the corresponding comparison. Terms represented by 10–500 genes are tested, including zero-hit terms in the Benjamini–Hochberg correction. Full test tables and gene-hit IDs are included for Table S9 export.

## Inputs and interpretation

The repository contains processed inputs and selected results. Large genotype files and upstream tools are external; command scripts document their required arguments. Demographic model templates, parameter definitions, likelihood summaries and bootstrap summaries are included in `data/demography/fsc/`.

[Introgression analyses](INTROGRESSION.md) describes the species-reference summaries and their component fields. Local fd characterizes allele sharing and does not by itself determine the evolutionary origin of variation at an individual locus.
