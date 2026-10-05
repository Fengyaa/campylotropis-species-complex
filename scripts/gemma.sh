#!/usr/bin/env bash
# Template: relatedness + genetic PC1-PC3, no environmental fixed effects.
# Covariates: no-header matrix with intercept + PC1-PC3, in BED/FAM sample order.
# BED input must meet final GWAS filters (MAF >0.05).
set -euo pipefail
[[ $# -eq 6 ]] || { echo "Usage: $0 BED_PREFIX PHENOTYPES KINSHIP PC_COVARIATES TRAIT_COLUMN OUTDIR" >&2; exit 2; }
mkdir "$6"
gemma -bfile "$1" -p "$2" -k "$3" -c "$4" -n "$5" -lmm 1 -outdir "$6" -o association
