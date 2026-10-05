#!/usr/bin/env bash
# Optional upstream rerun; processed 50/25 outputs are already bundled.
set -euo pipefail
[[ $# -eq 2 ]] || { echo "Usage: $0 INPUT.vcf.gz NEW_OUTDIR" >&2; exit 2; }
root=$(cd "$(dirname "$0")/.." && pwd)
vcf=$(cd "$(dirname "$1")" && printf '%s/%s' "$PWD" "$(basename "$1")")
mkdir "$2"
out=$(cd "$2" && pwd)
mkdir "$out/MY1" "$out/CY32"
cd "$out/MY1"
Dsuite Dinvestigate -w 50,25 "$vcf" "$root/data/introgression/sp_set_Cyu.txt" "$root/data/introgression/sp_set_Cyu_trio.tsv"
cd "$out/CY32"
Dsuite Dinvestigate -w 50,25 "$vcf" "$root/data/introgression/sp_set_Cpo.txt" "$root/data/introgression/sp_set_Cpo_trio.tsv"
