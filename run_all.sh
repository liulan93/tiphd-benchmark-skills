#!/usr/bin/env bash
# ============================================================================
# TiPhD Benchmark — run ALL algorithms (batch stage) with the right interpreter
#
# Usage:
#   # Activate an R >= 4.4 + Seurat 5 environment first so `Rscript` is on
#   # PATH (or point TIPHD_RSCRIPT at it). Then:
#   bash run_all.sh                       # run every tool's batch_run.py
#   bash run_all.sh music scab tirank     # run only the named tools
#   R_ONLY=1  bash run_all.sh             # R tools only
#   PY_ONLY=1 bash run_all.sh             # Python tools only
#
# After the batch stage, evaluate all tools:
#   bash run_evaluate_all.sh
#
# Environment / interpreter selection (override any of these):
#   TIPHD_RSCRIPT  path to Rscript (default: Rscript on PATH)
#   TORCH_PY       python for tiphd-torch tools (default: auto-resolve)
#   STATS_PY       python for tiphd-stats (TiRank)
#   PY310_PY       python for tiphd-py310 (Statescope)
#   TORCH_ENV / STATS_ENV / PY310_ENV   conda env names (default tiphd-*)
#   TIRANK_GPU=1   run TiRank on GPU (its default env may ship CPU torch)
#
# Auto-resolution order for each conda-env python:
#   1. the explicit *_PY variable
#   2. $CONDA_PREFIX/../../envs/<name>/bin/python if a conda is active
#   3. `conda run -n <name> which python` (works with any conda install)
#   4. common conda locations: ~/miniconda3, ~/anaconda3, ~/miniforge3, /opt/conda
# ============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# Map a tool name to its interpreter class:
#   R (system Rscript) | torch | stats | py310
tool_class() {
  case "$1" in
    music|scpp|scstar2|scab|scipac|scissor|scpas|pipet) echo R ;;
    scad|scdeal|scsurv|sctrend|sidish|scper) echo torch ;;
    tirank) echo stats ;;
    statescope) echo py310 ;;
    *) echo R ;;
  esac
}
R_TOOLS=(music scpp scstar2 scab scipac scissor scpas pipet)
TORCH_TOOLS=(scad scdeal scsurv sctrend sidish scper)
STATS_TOOLS=(tirank)
PY310_TOOLS=(statescope)

# Resolve the python executable of a named conda env, robust to non-standard
# conda install locations.
env_python() {
  local name="$1" explicit="$2"
  [ -n "$explicit" ] && { [ -x "$explicit" ] && { echo "$explicit"; return 0; } || return 1; }
  # 2. relative to an active conda prefix
  if [ -n "${CONDA_PREFIX:-}" ]; then
    local base="${CONDA_PREFIX%/envs/*}"
    [ -x "$base/envs/$name/bin/python" ] && { echo "$base/envs/$name/bin/python"; return 0; }
    [ "$CONDA_PREFIX" = "$base" ] && [ -x "$base/bin/python" ] && { echo "$base/bin/python"; return 0; }
  fi
  # 3. ask conda itself (works regardless of install path)
  if command -v conda >/dev/null 2>&1; then
    local p
    p="$(conda run -n "$name" bash -c 'command -v python' 2>/dev/null | tail -1)"
    [ -n "$p" ] && [ -x "$p" ] && { echo "$p"; return 0; }
  fi
  # 4. common locations
  local base
  for base in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" /opt/conda; do
    [ -x "$base/envs/$name/bin/python" ] && { echo "$base/envs/$name/bin/python"; return 0; }
  done
  return 1
}

# Build the ordered tool list from args or default order
if [ "$#" -gt 0 ]; then
  TOOLS=("$@")
else
  TOOLS=("${R_TOOLS[@]}" "${TORCH_TOOLS[@]}" "${STATS_TOOLS[@]}" "${PY310_TOOLS[@]}")
fi
[ "${R_ONLY:-0}" = 1 ]  && TOOLS=("${R_TOOLS[@]}")
[ "${PY_ONLY:-0}" = 1 ] && TOOLS=("${TORCH_TOOLS[@]}" "${STATS_TOOLS[@]}" "${PY310_TOOLS[@]}")

RS="${TIPHD_RSCRIPT:-Rscript}"
TORCH_PY="$(env_python "${TORCH_ENV:-tiphd-torch}" "${TORCH_PY:-}" 2>/dev/null || true)"
STATS_PY="$(env_python "${STATS_ENV:-tiphd-stats}" "${STATS_PY:-}" 2>/dev/null || true)"
PY310_PY="$(env_python "${PY310_ENV:-tiphd-py310}" "${PY310_PY:-}" 2>/dev/null || true)"

echo "=== TiPhD run_all: ${#TOOLS[@]} tools ==="
echo "Rscript     : $($RS --version 2>&1 | head -1 || echo 'NOT FOUND')"
echo "torch python: ${TORCH_PY:-<env not found>}"
echo "stats python: ${STATS_PY:-<env not found>}"
echo "py310 python: ${PY310_PY:-<env not found>}"
echo

failed=()
for tool in "${TOOLS[@]}"; do
  [ -f "$HERE/$tool/batch_run.py" ] || { echo "SKIP $tool (no batch_run.py)"; continue; }
  cls="$(tool_class "$tool")"
  echo "----- [$tool] class=$cls -----"
  # IMPORTANT: invoke batch_run.py by path from the CURRENT working directory
  # (do NOT cd into the tool dir). batch_run.py locates its own scripts via
  # __file__, while data/ and results/ resolve from cwd ($TIPHD_DATA_DIR /
  # $TIPHD_OUT_DIR, default ./data and ./results under the project root).
  case "$cls" in
    R)
      if ! $RS --version >/dev/null 2>&1; then
        echo "  FAIL: Rscript not found (activate an R 4.4 + Seurat 5 env or set TIPHD_RSCRIPT)"
        failed+=("$tool"); continue
      fi
      python "$HERE/$tool/batch_run.py" || failed+=("$tool")
      ;;
    torch)
      if [ -z "$TORCH_PY" ]; then
        echo "  FAIL: tiphd-torch env not found (set TORCH_PY or create the env)"; failed+=("$tool"); continue
      fi
      "$TORCH_PY" "$HERE/$tool/batch_run.py" || failed+=("$tool")
      ;;
    stats)
      if [ -z "$STATS_PY" ]; then
        echo "  FAIL: tiphd-stats env not found (set STATS_PY or create the env)"; failed+=("$tool"); continue
      fi
      TIRANK_GPU="${TIRANK_GPU:-0}" "$STATS_PY" "$HERE/$tool/batch_run.py" || failed+=("$tool")
      ;;
    py310)
      if [ -z "$PY310_PY" ]; then
        echo "  FAIL: tiphd-py310 env not found (set PY310_PY or create the env)"; failed+=("$tool"); continue
      fi
      "$PY310_PY" "$HERE/$tool/batch_run.py" || failed+=("$tool")
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
