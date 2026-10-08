#!/usr/bin/env bash
# BayPass auxiliary-model command template.
# GENOTYPE_COUNTS must contain ALL final filtered genome-wide SNPs, not the structure panel.
set -euo pipefail
[[ $# -eq 5 ]] || { echo "Usage: $0 ALL_SNP_COUNTS CLIMATE_PC_MATRIX OMEGA NPOP OUTDIR" >&2; exit 2; }
[[ "$4" =~ ^[1-9][0-9]*$ ]] || { echo 'NPOP must be a positive integer' >&2; exit 2; }
mkdir "$5"
g_baypass -npop "$4" -gfile "$1" -efile "$2" -scalecov -omegafile "$3" -auxmodel -outprefix "$5/aux"
