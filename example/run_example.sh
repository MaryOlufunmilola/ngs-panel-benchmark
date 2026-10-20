#!/usr/bin/env bash
# Run the benchmark's Python steps (panel BED, coverage, call comparison,
# titration summary) on tiny synthetic inputs and compare the result with
# example/expected/titration_metrics.csv. Needs only Python with pandas and
# matplotlib; no reference genome, reads, samtools, or Docker. Takes seconds.
#
# Usage: bash example/run_example.sh
set -euo pipefail
EX="$(cd "$(dirname "$0")" && pwd)"
CLI="$EX/../src/cli.py"
IN="$EX/toy_inputs"
WORK="$EX/work"
rm -rf "$IN" "$WORK"; mkdir -p "$WORK/data" "$WORK/results"
export MPLBACKEND=Agg

python "$EX/make_toy_inputs.py" "$IN"

echo "==> 1. Panel BED from GTF"
python "$CLI" make-panel --gtf "$IN/toy.gtf.gz" --genes "$IN/genes.txt" --out "$WORK/data/panel.bed"

for FRAC in 1.0 0.5 0.25 0.1; do
  D="$WORK/results/frac_${FRAC}"; mkdir -p "$D"
  echo "==> Depth fraction $FRAC: coverage and comparison"
  python "$CLI" coverage --perbase "$IN/frac_${FRAC}/cov.per-base.bed.gz" --panel "$WORK/data/panel.bed" \
    --out-prefix "$D/coverage" > "$D/coverage.log"
  python "$CLI" compare --truth "$IN/truth.vcf.gz" --query "$IN/frac_${FRAC}/calls.vcf.gz" \
    --bed "$WORK/data/panel.bed" --out "$D/compare.csv" > /dev/null
done

echo "==> Titration summary"
python "$CLI" collect --results "$WORK/results" --fractions 1.0 0.5 0.25 0.1

echo "==> Comparing with example/expected/titration_metrics.csv"
python - "$WORK/results/titration_metrics.csv" "$EX/expected/titration_metrics.csv" <<'EOF'
import sys
import pandas as pd
got, exp = pd.read_csv(sys.argv[1]), pd.read_csv(sys.argv[2])
pd.testing.assert_frame_equal(got.reset_index(drop=True), exp.reset_index(drop=True), check_exact=False, rtol=1e-6)
print("Output matches expected/titration_metrics.csv. Plots: example/work/results/*.png")
EOF
