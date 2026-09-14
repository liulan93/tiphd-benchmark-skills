#!/usr/bin/env bash
# ============================================================================
# TiPhD Benchmark — run the evaluate.R stage for ALL (or selected) tools.
#
# Usage (run from the project root that contains ./data and ./results):
#   conda activate <r-env>     # must expose Rscript >= 4.4 with Seurat 5
#   bash run_evaluate_all.sh
#   bash run_evaluate_all.sh music scab
#
# Produces five *_final_metrics.csv per tool (one per cancer), when the
# corresponding pair outputs exist. evaluate.R locates its own config via
# __file__ but reads data/results from cwd (TIPHD_DATA_DIR / TIPHD_OUT_DIR),
# so it is invoked by path WITHOUT changing directory.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ALL=(music scpp scstar2 scab scipac scissor scpas pipet \
     scad scdeal scsurv sctrend sidish scper statescope tirank)
TOOLS=("$@")
[ "$#" -eq 0 ] && TOOLS=("${ALL[@]}")

RS="${TIPHD_RSCRIPT:-Rscript}"
echo "=== evaluate_all using: $($RS --version 2>&1) ==="
echo "    cwd: $(pwd)  (expecting ./data and ./results here)"

for tool in "${TOOLS[@]}"; do
  [ -f "$HERE/$tool/evaluate.R" ] || { echo "SKIP $tool (no evaluate.R)"; continue; }
  echo "----- [$tool] evaluate -----"
  "$RS" "$HERE/$tool/evaluate.R" || echo "  evaluate reported failures for $tool"
done
echo "=== evaluate_all complete ==="
