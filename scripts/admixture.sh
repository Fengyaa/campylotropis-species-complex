#!/usr/bin/env bash
# Cleaned author commands; final MAF 0.01, call rate 0.90, intergenic, 1-kb spacing.
# Usage: bash scripts/admixture.sh INPUT.vcf.gz INTERGENIC.pos OUTDIR
set -euo pipefail
[[ $# -eq 3 ]] || { echo "Usage: $0 INPUT.vcf.gz INTERGENIC.pos OUTDIR" >&2; exit 2; }
input=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
positions=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
mkdir "$3"
out=$(cd "$3" && pwd)
vcftools --gzvcf "$input" --positions "$positions" --maf 0.01 --max-missing 0.9 --thin 1000 \
 --remove-indv 2014-128-7A --remove-indv Gln2014-8-3A --remove-indv xubo266-7A \
 --recode --out "$out/259.maf01ms90thin1k"
vcftools --vcf "$out/259.maf01ms90thin1k.recode.vcf" --remove-indv XB_DR_C --remove-indv xubo1407 \
 --recode --out "$out/257.maf01ms90thin1k"
vcftools --vcf "$out/257.maf01ms90thin1k.recode.vcf" --plink --out "$out/257.maf01ms90thin1k"
plink --file "$out/257.maf01ms90thin1k" --allow-extra-chr --recode12 --out "$out/final_maf01ms90thin1k"
awk '{print $2}' "$out/final_maf01ms90thin1k.ped" > "$out/sample_order.txt"
mkdir "$out/admixture"
cd "$out/admixture"
# Single CV-enabled fit per K avoids overwriting an earlier non-CV fit.
for k in $(seq 1 20); do
 admixture --cv ../final_maf01ms90thin1k.ped "$k" | tee "admixture_${k}.log"
done
