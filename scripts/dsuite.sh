#!/usr/bin/env bash
# Dsuite command template: 50 usable SNPs per window and a 25-SNP step.
set -euo pipefail
[[ $# -eq 5 ]] || { echo "Usage: $0 INPUT.vcf.gz SETS.txt ROOTED_TREE.nwk TRIOS.tsv OUTDIR" >&2; exit 2; }
absolute() { printf '%s/%s\n' "$(cd "$(dirname "$1")" && pwd)" "$(basename "$1")"; }
vcf=$(absolute "$1");sets=$(absolute "$2");tree=$(absolute "$3");trios=$(absolute "$4")
mkdir "$5"
cd "$5"
Dsuite Dtrios -t "$tree" -o dtrios "$vcf" "$sets"
Dsuite Fbranch "$tree" dtrios_tree.txt > fbranch.tsv
Dsuite Dinvestigate -w 50,25 "$vcf" "$sets" "$trios"
