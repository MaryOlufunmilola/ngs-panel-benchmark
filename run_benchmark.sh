#!/usr/bin/env bash
# Variant calling benchmark for a targeted panel on GIAB HG002.
# Builds a panel of 50 genes, extracts reads over it, measures coverage and uniformity, calls variants,
# compares them with the GIAB truth set, and repeats at lower depths. Run from the repository root.
set -euo pipefail
cd "$(dirname "$0")"
source config/settings.sh
mkdir -p data results logs

for tool in samtools bcftools mosdepth bedtools python; do
  command -v "$tool" >/dev/null || { echo "Missing tool: $tool (see environment.yml)"; exit 1; }
done
for f in "$REF" "$REF.fai" "$GTF" "$TRUTH_VCF" "$TRUTH_BED"; do
  [ -e "$f" ] || { echo "Missing input: $f (see README, Data)"; exit 1; }
done
USE_HAPPY=false
HAPPY_IMAGE="${HAPPY_IMAGE:-quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0}"; HAPPY_CMD="${HAPPY_CMD:-hap.py}"
if command -v docker >/dev/null && docker image inspect "$HAPPY_IMAGE" >/dev/null 2>&1; then USE_HAPPY=true; fi
echo "hap.py available: $USE_HAPPY (otherwise exact match comparison is used)"

echo "==> 1. Panel BED from GTF"
python src/cli.py make-panel --gtf "$GTF" --out data/panel.bed | tee logs/panel.log

echo "==> 2. Reads over the panel"
if [ ! -s data/hg002_panel.bam ]; then
  samtools view -@ "$THREADS" -b -M -L data/panel.bed -T "$REF" -o data/hg002_panel.bam "$READS"
  samtools index data/hg002_panel.bam
fi

echo "==> 3. Evaluation regions (panel within GIAB high confidence regions) and truth subset"
bedtools intersect -a data/panel.bed -b "$TRUTH_BED" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > data/eval.bed
bcftools view -R data/eval.bed "$TRUTH_VCF" -Ou | bcftools norm -f "$REF" -m -any -Oz -o data/truth_panel.vcf.gz
bcftools index -f -t data/truth_panel.vcf.gz
echo "Evaluation regions: $(awk '{s+=$3-$2} END {print s}' data/eval.bed) bp"

for FRAC in $FRACTIONS; do
  D="results/frac_${FRAC}"; mkdir -p "$D"
  echo "==> 4. Depth fraction $FRAC"
  if [ "$FRAC" = "1.0" ]; then
    BAM=data/hg002_panel.bam
  else
    BAM="$D/reads.bam"
    samtools view -@ "$THREADS" -b -s "${SEED}${FRAC#0}" -o "$BAM" data/hg002_panel.bam   # e.g. -s 42.25
    samtools index "$BAM"
  fi

  mosdepth -t "$THREADS" --by data/panel.bed "$D/cov" "$BAM"
  python src/cli.py coverage --perbase "$D/cov.per-base.bed.gz" --out-prefix "$D/coverage" > "$D/coverage.log"

  bcftools mpileup -f "$REF" -R data/panel.bed -a AD,DP -q 20 -Q 20 -d 1000 -L 1000 -Ou "$BAM" \
    | bcftools call -mv -Ou | bcftools norm -f "$REF" -m -any -Oz -o "$D/calls.vcf.gz"
  bcftools index -f -t "$D/calls.vcf.gz"

  if $USE_HAPPY; then
    docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w "$HAPPY_IMAGE" "$HAPPY_CMD" \
      data/truth_panel.vcf.gz "$D/calls.vcf.gz" -f data/eval.bed -r "$REF" -o "$D/happy" > "$D/happy.log" 2>&1
  else
    python src/cli.py compare --truth data/truth_panel.vcf.gz --query "$D/calls.vcf.gz" \
      --bed data/eval.bed --out "$D/compare.csv"
  fi
done

echo "==> 5. Titration summary and plots"
python src/cli.py collect --results results --fractions $FRACTIONS | tee logs/titration.log

if [ -n "$CLINVAR_VCF" ] && [ -e "$CLINVAR_VCF" ]; then
  echo "==> 6. ClinVar annotation of calls at full depth"
  # ClinVar names chromosomes 1, 2, ...; the GRCh38 analysis set uses chr1, chr2, ...
  bcftools view -h "$CLINVAR_VCF" | grep -q '^##contig=<ID=chr' && CV="$CLINVAR_VCF" || {
    for c in $(seq 1 22) X Y; do echo "$c chr$c"; done > data/chr_map.txt
    echo "MT chrM" >> data/chr_map.txt
    bcftools annotate --rename-chrs data/chr_map.txt "$CLINVAR_VCF" -Oz -o data/clinvar_chr.vcf.gz
    bcftools index -f -t data/clinvar_chr.vcf.gz; CV=data/clinvar_chr.vcf.gz; }
  bcftools annotate -a "$CV" -c INFO/CLNSIG,INFO/CLNDN results/frac_1.0/calls.vcf.gz -Oz -o results/calls_clinvar.vcf.gz
  bcftools query -i 'INFO/CLNSIG!="."' -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
    results/calls_clinvar.vcf.gz > results/clinvar_matches.tsv || true
  echo "Calls with a ClinVar record: $(wc -l < results/clinvar_matches.tsv)"
fi

echo "Done. See results/titration_metrics.csv, results/*.png, and logs/."
