#!/usr/bin/env bash
# install_python_packages.sh — pip install 自定义 Python 包（幂等）
# 用法:
#   conda activate tiphd-torch  && bash install_python_packages.sh   # 深度学习类包
#   conda activate tiphd-stats  && bash install_python_packages.sh   # 统计/图模型类包
#
# 脚本根据当前激活的 conda 环境名自动选择装哪些包。
#
# 说明：
#   SIDISH 的 setup.py 锁定了严格的版本号（部分依赖要求 Python 3.10+），
#   在 Python 3.9 环境下直接 `pip install -e` 会解析依赖失败。这里改用
#   --no-deps 安装包本身，再手动安装它实际 import 所需的运行时依赖。

set -e
HERE="$(cd "$(dirname "$0")" && pwd)"

CURRENT_ENV="${CONDA_DEFAULT_ENV:-unknown}"
echo "=== Installing custom Python packages (editable, idempotent) ==="
echo "    conda env: $CURRENT_ENV"

TORCH_PKGS=(SIDISH scTREND scSurv scPER)
STATS_PKGS=(TiRank)

# 普通自定义包：直接 editable 安装
install_pkg() {
  local pkg="$1"
  local import_name="$2"
  import_name="${import_name:-$pkg}"
  if python -c "import $import_name" 2>/dev/null; then
    echo "  OK   $pkg already installed, skipping"
    return 0
  fi
  local SRC_DIR="$HERE/third_party/python/$pkg"
  if [ ! -d "$SRC_DIR" ]; then
    echo "  SKIP $pkg: source dir not found ($SRC_DIR)"
    return 0
  fi
  echo "  ...  installing $pkg from $SRC_DIR"
  pip install -e "$SRC_DIR"
}

# SIDISH 特殊处理：--no-deps 安装包体，再补装运行时依赖
install_sidish() {
  if python -c "from SIDISH.SIDISH import SIDISH" 2>/dev/null; then
    echo "  OK   SIDISH already installed, skipping"
    return 0
  fi
  local SRC_DIR="$HERE/third_party/python/SIDISH"
  echo "  ...  installing SIDISH (--no-deps) + runtime deps"
  # SIDISH 实际 import 用到的第三方库
  pip install pyro-ppl torchvision torch_geometric shap bioinfokit \
              imbalanced-learn lifelines statsmodels 2>/dev/null || true
  pip install -e "$SRC_DIR" --no-deps
}

case "$CURRENT_ENV" in
  tiphd-torch)
    echo "    → installing deep-learning packages: ${TORCH_PKGS[*]}"
    install_sidish
    for pkg in scTREND scSurv scPER; do install_pkg "$pkg"; done
    ;;
  tiphd-stats)
    echo "    → installing stats/graph packages: ${STATS_PKGS[*]}"
    # TiRank 没有标准 setup 安装，依赖通过 PYTHONPATH 指向源码目录
    install_pkg TiRank
    install_pkg scPER
    ;;
  *)
    echo "    → unknown env '$CURRENT_ENV', installing all packages"
    install_sidish
    for pkg in scTREND scSurv scPER; do install_pkg "$pkg"; done
    install_pkg TiRank
    ;;
esac

echo "=== Done ==="
echo "提示: TiRank 通过 batch_run.py 自动把 third_party/python/TiRank 加入"
echo "      PYTHONPATH，无需 pip 安装；其依赖 lifelines/optuna/leidenalg/"
echo "      python-igraph 由 tiphd-stats conda 环境提供。"
