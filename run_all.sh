#!/usr/bin/env bash
# ============================================================================
# TiPhD Benchmark — run ALL algorithms (batch stage) with the right interpreter
#
# Usage:
#   conda activate <your-r-env-with-rscript>   # must expose `Rscript` on PATH
#   bash run_all.sh                            # run every tool's batch_run.py
#   bash run_all.sh music scab tirank          # run only the named tools
#   R_ONLY=1   bash run_all.sh                # R tools only
#   PY_ONLY=1  bash run_all.sh                # Python tools only
#
# After the batch stage finishes, run the evaluation for each tool:
#   bash run_evaluate_all.sh
#
# Environment selection (override as needed):
#   TORCH_ENV  conda env name for tiphd-torch  (default: tiphd-torch)
#   STATS_ENV  conda env name for tiphd-stats (default: tiphd-stats)
#   PY310_ENV  conda env name for tiphd-py310  (default: tiphd-py310)
#   The `Rscript` used for R tools and evaluate.R is taken from PATH
#   (or TIPHD_RSCRIPT). Activate an R >= 4.4 + Seurat 5 environment first.
# ============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# tool -> interpreter class: R (system Rscript) | torch | stats | py310
declare -A ENV_OF=(
  [music]=R [scpp]=R [scstar2]=R [scab]=R [scipac]=R [scissor]=R [scpas]=R [pipet]=R
  [scad]=torch [scdeal]=torch [scsurv]=torch [sctrend]=torch [sidish]=torch
  [tirank]=stats
  [statescope]=py310
  [scper]=torch
)
R_TOOLS=(music scpp scstar2 scab scipac scissor scpas pipet)
TORCH_TOOLS=(scad scdeal scsurv sctrend sidish scper)
STATS_TOOLS=(tirank)
PY310_TOOLS=(statescope)

# Resolve the python executable of a conda env by name (works even if not active)
env_python() {
  local name="$1"
  local cand=""
  for base in "${CONDA_PREFIX%/envs/*}" "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" /opt/conda; do
    [ -z "$base" ] && continue
    [ -x "$base/envs/$name/bin/python" ] && cand="$base/envs/$name/bin/python" && break
    [ -x "$base/bin/python" ] && [ "$name" = "base" ] && cand="$base/bin/python" && break
  done
  echo "$cand"
}

# Build the ordered tool list from args or default order
if [ "$#" -gt 0 ]; then
  TOOLS=("$@")
else
  TOOLS=("${R_TOOLS[@]}" "${TORCH_TOOLS[@]}" "${STATS_TOOLS[@]}" "${PY310_TOOLS[@]}")
fi
[ "${R_ONLY:-0}" = 1 ]   && TOOLS=("${R_TOOLS[@]}")
[ "${PY_ONLY:-0}" = 1 ]  && TOOLS=("${TORCH_TOOLS[@]}" "${STATS_TOOLS[@]}" "${PY310_TOOLS[@]}")

TORCH_PY="$(env_python "${TORCH_ENV:-tiphd-torch}")"
STATS_PY="$(env_python "${STATS_ENV:-tiphd-stats}")"
PY310_PY="$(env_python "${PY310_ENV:-tiphd-py310}")"

echo "=== TiPhD run_all: ${#TOOLS[@]} tools ==="
echo "Rscript: ${TIPHD_RSCRIPT:-$(command -v Rscript || echo 'NOT ON PATH')}"
echo "torch python : ${TORCH_PY:-<env not found>}"
echo "stats python : ${STATS_PY:-<env not found>}"
echo "py310 python : ${PY310_PY:-<env not found>}"
echo

failed=()
for tool in "${TOOLS[@]}"; do
  [ -f "$tool/batch_run.py" ] || { echo "SKIP $tool (no batch_run.py)"; continue; }
  cls="${ENV_OF[$tool]:-R}"
  echo "----- [$tool] class=$cls -----"
  case "$cls" in
    R)
      "${TIPHD_RSCRIPT:-Rscript}" --version >/dev/null 2>&1 || { echo "  FAIL: Rscript not found"; failed+=("$tool"); continue; }
      ( cd "$tool" && python batch_run.py ) || failed+=("$tool")
      ;;
    torch)
      [ -z "$TORCH_PY" ] && { echo "  FAIL: tiphd-torch env not found"; failed+=("$tool"); continue; }
      ( cd "$tool" && "$TORCH_PY" batch_run.py ) || failed+=("$tool")
      ;;
    stats)
      [ -z "$STATS_PY" ] && { echo "  FAIL: tiphd-stats env not found"; failed+=("$tool"); continue; }
      ( cd "$tool" && TIRANK_GPU="${TIRANK_GPU:-0}" "$STATS_PY" batch_run.py ) || failed+=("$tool")
      ;;
    py310)
      [ -z "$PY310_PY" ] && { echo "  FAIL: tiphd-py310 env not found"; failed+=("$tool"); continue; }
      ( cd "$tool" && "$PY310_PY" batch_run.py ) || failed+=("$tool")
      ;;
  esac
done

echo
echo "=== run_all complete ==="
if [ "${#failed[@]}" -gt 0 ]; then
  echo "Tools with failures (rerun is idempotent — existing outputs are skipped):"
  printf '  %s\n' "${failed[@]}"
  exit 1
fi
echo "All selected tools finished their batch stage."
