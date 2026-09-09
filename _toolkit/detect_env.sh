#!/usr/bin/env bash
# detect_env.sh — 探测 R / Python / conda 环境（供 setup-env 与各算法 skill 前置调用）
# 用法: bash detect_env.sh
# 退出码: 0 = 全部就绪；非 0 = 有缺失（仍会输出探测结果）

set -u
fail=0

echo "=== tiphd-torch/stats/R env probe ==="

# R
if command -v R >/dev/null 2>&1; then
  R_VER=$(R --version 2>/dev/null | head -1)
  echo "R       : OK  ($R_VER)"
else
  echo "R       : MISSING (https://cran.r-project.org/)"
  fail=1
fi

# Rscript + R package versions
if command -v Rscript >/dev/null 2>&1; then
  RS_VER=$(Rscript --version 2>&1 | head -1)
  echo "Rscript : OK  ($RS_VER)"
  if [ -n "${TIPHD_RSCRIPT:-}" ]; then
    echo "         : TIPHD_RSCRIPT=$TIPHD_RSCRIPT (覆盖默认)"
  fi
  # 检查 R 主版本号（需要 ≥ 4.4）
  R_MAJOR=$(Rscript -e 'cat(R.version$major, ".", R.version$minor, sep="")' 2>/dev/null)
  R_MAJOR_NUM=$(echo "$R_MAJOR" | cut -d. -f1)
  R_MINOR_NUM=$(echo "$R_MAJOR" | cut -d. -f2)
  if [ "$R_MAJOR_NUM" -lt 4 ] || { [ "$R_MAJOR_NUM" -eq 4 ] && [ "$R_MINOR_NUM" -lt 4 ]; }; then
    echo "         : WARNING R $R_MAJOR < 4.4 — Seurat v5 数据需要 R ≥ 4.4.0"
    fail=1
  fi
  # 检查 Seurat 版本
  SEURAT_VER=$(Rscript -e 'cat(if(requireNamespace("Seurat",quietly=TRUE)) as.character(packageVersion("Seurat")) else "MISSING")' 2>/dev/null)
  if [ "$SEURAT_VER" = "MISSING" ]; then
    echo "Seurat  : MISSING (install.packages('Seurat'))"
    fail=1
  else
    SEURAT_MAJOR=$(echo "$SEURAT_VER" | cut -d. -f1)
    if [ "$SEURAT_MAJOR" -lt 5 ] 2>/dev/null; then
      echo "Seurat  : $SEURAT_VER (WARNING < 5.0 — 无法正确读取 Assay5 .rds 数据)"
    else
      echo "Seurat  : $SEURAT_VER (OK)"
    fi
  fi
  # 检查 Matrix 版本（Seurat 5 需要 ≥ 1.6.4）
  MATRIX_VER=$(Rscript -e 'cat(if(requireNamespace("Matrix",quietly=TRUE)) as.character(packageVersion("Matrix")) else "MISSING")' 2>/dev/null)
  echo "Matrix  : $MATRIX_VER"
else
  echo "Rscript : MISSING (Windows 用户可设 export TIPHD_RSCRIPT=/path/to/Rscript.exe)"
  fail=1
fi

# Windows-only: Rtools44
if uname -s 2>/dev/null | grep -qiE "mingw|msys|cygwin"; then
  if command -v gcc >/dev/null 2>&1; then
    GCC_VER=$(gcc --version 2>/dev/null | head -1)
    echo "gcc     : OK  ($GCC_VER)  [Rtools44]"
  else
    echo "gcc     : MISSING  [需 Rtools44 才能编译 scPAS / Scissor 的 C++ 代码]"
    fail=1
  fi
fi

# Python
if command -v python >/dev/null 2>&1; then
  PY_VER=$(python --version 2>&1)
  echo "python  : OK  ($PY_VER)"
else
  echo "python  : MISSING"
  fail=1
fi

# pip
if command -v pip >/dev/null 2>&1; then
  echo "pip     : OK"
else
  echo "pip     : MISSING"
  fail=1
fi

# conda
if command -v conda >/dev/null 2>&1; then
  CONDA_VER=$(conda --version 2>&1)
  echo "conda   : OK  ($CONDA_VER)"
  for env in tiphd-torch tiphd-stats; do
    if conda env list 2>/dev/null | grep -qE "^\s*${env}\s"; then
      echo "env $env: OK  (exists)"
    else
      echo "env $env: MISSING (会从 environment*.yml 创建)"
    fi
  done
else
  echo "conda   : MISSING (optional, but recommended for Python deps)"
fi

# modelscope (setup-data 用)
if python -c "import modelscope" 2>/dev/null; then
  echo "modelscope: OK"
else
  echo "modelscope: MISSING (run: pip install modelscope)"
  fail=1
fi

echo "=== exit: $fail ==="
exit $fail
