#!/usr/bin/env bash
# ============================================================================
# TiPhD Benchmark — run the evaluate.R stage for ALL (or selected) tools.
#
# Usage:
#   conda activate <r-env>     # must expose Rscript >= 4.4 with Seurat 5
#   bash run_evaluate_all.sh
#   bash run_evaluate_all.sh music scab
#
# Produces five *_final_metrics.csv per tool (one per cancer), when the
# corresponding pair outputs exist.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

ALL=(music scpp scstar2 scab scipac scissor scpas pipet \
     scad scdeal scsurv sctrend sidish scper statescope tirank)
TOOLS=("$@")
[ "$#" -eq 0 ] && TOOLS=("${ALL[@]}")

RS="${TIPHD_RSCRIPT:-Rscript}"
echo "=== evaluate_all using: $($RS --version 2>&1) ==="

for tool in "${TOOLS[@]}"; do
  [ -f "$tool/evaluate.R" ] || { echo "SKIP $tool (no evaluate.R)"; continue; }
  echo "----- [$tool] evaluate -----"
  ( cd "$tool" && "$RS" evaluate.R ) || echo "  evaluate reported failures for $tool"
done
echo "=== evaluate_all complete ==="
