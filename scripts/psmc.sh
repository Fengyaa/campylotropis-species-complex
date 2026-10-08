#!/usr/bin/env bash
# Infer PSMC trajectories from per-individual gVCFs.
# Requires an external compatible vcf2fq.py converter.
# Usage: bash scripts/psmc.sh SAMPLE_LIST REF.fa GVCF_DIR VCF2FQ.py PSMC_DIR OUTDIR
set -euo pipefail
[[ $# -eq 6 ]] || { echo "Usage: $0 SAMPLE_LIST REF.fa GVCF_DIR VCF2FQ.py PSMC_DIR OUTDIR" >&2; exit 2; }
samples=$1;reference=$2;gvcfs=$3;converter=$4;psmc_dir=$5;out=$6
[[ -f "$samples" && -f "$reference" && -f "$converter" ]] || { echo 'Missing input/converter' >&2; exit 2; }
mkdir "$out"
while IFS= read -r sample || [[ -n "$sample" ]]; do
 sample=${sample%$'\r'}
 [[ -n "$sample" ]] || continue
 gatk --java-options '-XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -Xmx28g' GenotypeGVCFs \
  -R "$reference" -V "$gvcfs/${sample}_MD.g.vcf.gz" -O "$out/$sample.raw.vcf.gz" -all-sites
 gzip -dc "$out/$sample.raw.vcf.gz" | python "$converter" | gzip > "$out/$sample.fq.gz"
 "$psmc_dir/utils/fq2psmcfa" -q20 "$out/$sample.fq.gz" > "$out/$sample.psmcfa"
 "$psmc_dir/psmc" -N25 -t15 -r5 -p '4+25*2+4+6' -o "$out/$sample.psmc" "$out/$sample.psmcfa"
 perl "$psmc_dir/utils/psmc_plot.pl" -u 7e-9 -g 2 -R -p "$out/$sample" "$out/$sample.psmc"
done < "$samples"
