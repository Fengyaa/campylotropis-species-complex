# Final settings and provenance

Author confirmations on 5 October 2026 supersede the earlier archive notes.

| Topic | Final description |
|---|---|
| Climate | Morphology RDA and genomic pRDA share bio3, bio7, bio10, **bio13**, bio15 and July solar radiation. The proposed bio12 change was explicitly corrected to bio13 by the author. |
| Climate sources | 19 WorldClim 1.4 bioclimatic variables plus 36 WorldClim 2.1 monthly solar/wind/vapor variables, 2.5 arc-minute resolution. Shared population/individual values are preserved. |
| Morphology | 192 complete individuals selected from the 195 GL-GWAS cohort; eight traits and 28 populations. R² and adjusted R² are separately reported. |
| Genomic pRDA | 257 individuals, conditioning on PC1–PC3; six constrained axes; final loading threshold 3.5 SD and absolute marginal correlation >0.5. |
| BayPass | Association testing uses the full filtered genome-wide SNP dataset. This clarification does not independently identify the panel used to estimate the supplied covariance matrix. |
| GEMMA | Relatedness and genetic PC covariates; no environmental fixed effects. GL uses 195 non-missing phenotypes and 3,469,695 tests. |
| Dinvestigate | Author confirms 50 usable SNPs/25 step for both scales. The historical species source was named Cma_Cyu_Cpo_localFstats_100_50.txt; retained values are not altered or relabeled as a newly executed run. |
| ADMIXTURE | structurePlot.R reads 257.snps.hardfiltered.maf01ms90thin1k.K.Q and writes the selected maf01ms09_re.pdf. Prepared with intergenic positions, MAF 0.01, max-missing 0.9 and thinning 1 kb. Alternative MAF panels omitted. |
| SVDquartets | Fourfold-degenerate thin500 alignment; 259 positions; supplied population partitions; outgroups 131/244; 100 bootstrap replicates. SPECIES is the partition name for populations. |
| PSMC | GATK 4.2.5.0 all-sites calls, original vcf2fq.py, fq2psmcfa -q20; -N25 -t15 -r5 -p 4+25*2+4+6. Supplied plot command uses μ=8.17e-8 and g=10. |
| Gene overlap | At least 1 bp of the full annotated gene span, including introns/exons; distinct genes within each GO candidate set. Overlapping intervals are merged for coverage. |
| Window comparisons | Four upper-tail cutoffs: 10%, 5%, 2.5%, 1%; eligible sites >=100; finite values; 0<fd<1. Final GEA comparisons retain the 22-statistic BH family, 9,999 shifts and 1,000 bootstrap replicates. |

No full demographic fitting, new VCF scan or rescaling of manuscript demographic dates was performed during this curation. Old differently calibrated PSMC plot exports were omitted. FSC input/model/likelihood summaries are retained without inventing missing optimization commands.

The selected archive excludes DILS, SEM, balancing-selection exploration, old phenotype/RDA versions, other ADMIXTURE panels, alternate GEA cutoffs and the earlier overall reproduction framework. Large raw/genotype files and original conversion/configuration inputs are supplied externally where noted in the guide.
