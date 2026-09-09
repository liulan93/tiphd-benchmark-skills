"""
download_data.py — 从 ModelScope 下载 TiPhD 测试数据集到本地
用法:
    python download_data.py                          # 默认 <cwd>/data
    python download_data.py --target-dir /path      # 自定义目录
幂等: 若 <target>/.tiphd_data_ready 已存在则跳过。

源: https://www.modelscope.cn/datasets/fansailing/TiPhD_test_data/files
依赖: pip install modelscope
"""
import argparse
import os
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description="Download TiPhD test dataset from ModelScope")
    parser.add_argument(
        "--target-dir",
        default=os.environ.get("TIPHD_DATA_DIR") or os.path.join(os.getcwd(), "data"),
        help="Target directory (default: $TIPHD_DATA_DIR or <cwd>/data)",
    )
    parser.add_argument(
        "--repo-id",
        default="fansailing/TiPhD_test_data",
        help="ModelScope dataset repo id",
    )
    parser.add_argument(
        "--revision",
        default="master",
        help="ModelScope repo revision (branch / tag)",
    )
    args = parser.parse_args()

    target = os.path.abspath(args.target_dir)
    os.makedirs(target, exist_ok=True)

    marker = os.path.join(target, ".tiphd_data_ready")
    if os.path.exists(marker):
        print(f"[tiphd] data already ready at {target}  (marker: {marker})")
        print("       delete the marker to force re-download.")
        return 0

    try:
        from modelscope import snapshot_download
    except ImportError:
        print("[tiphd] ERROR: modelscope not installed. Run: pip install modelscope",
              file=sys.stderr)
        return 1

    print(f"[tiphd] downloading {args.repo_id} (rev={args.revision}) -> {target}")
    try:
        path = snapshot_download(
            repo_id=args.repo_id,
            repo_type="dataset",
            local_dir=target,
            revision=args.revision,
        )
        print(f"[tiphd] downloaded to: {path}")
    except Exception as e:
        print(f"[tiphd] ERROR: download failed: {e}", file=sys.stderr)
        return 1

    # 完整性验证：5 癌种核心文件 + 5 金标准
    expected = [
        # AML
        "scRNA/AML/GSE116256.rds",
        "scRNA/AML/GSE116256.h5ad",
        "bulk/AML/exp/TCGA_exp.csv",
        "bulk/AML/exp/wave12_exp.csv",
        "bulk/AML/exp/wave34_exp.csv",
        "bulk/AML/clinical/TCGA_clinical.csv",
        # CRC
        "scRNA/CRC/GSE144735.rds",
        "scRNA/CRC/GSE132465.rds",
        "bulk/CRC/exp/GSE14333_exp.csv",
        "bulk/CRC/exp/GSE17536_exp.csv",
        "bulk/CRC/exp/GSE33113_exp.csv",
        "bulk/CRC/exp/GSE37892_exp.csv",
        "bulk/CRC/exp/GSE39582_exp.csv",
        # HCC
        "scRNA/HCC/GSE149614.h5ad",
        "scRNA/HCC/GSE151530.h5ad",
        "bulk/HCC/exp/GSE116174_exp.csv",
        # LUAD
        "scRNA/LUAD/GSE127465.h5ad",
        "scRNA/LUAD/GSE131907.h5ad",
        "scRNA/LUAD/GSE148071.h5ad",
        # GC
        "scRNA/GC/GSE183904.h5ad",
        "bulk/GC/exp/GSETCGA_exp.csv",
        # gold standard
        "gold_standard_AML.csv",
        "gold_standard_all_CRC.csv",
        "gold_standard_all_HCC.csv",
        "gold_standard_all_LUDA.csv",  # NOTE: 原始拼写（不是 LUAD）
        "gold_standard_all_GC.csv",
    ]
    missing = [p for p in expected if not os.path.exists(os.path.join(target, p))]
    if missing:
        print(f"[tiphd] WARNING: {len(missing)} expected files missing:", file=sys.stderr)
        for p in missing[:10]:
            print(f"          - {p}", file=sys.stderr)
        if len(missing) > 10:
            print(f"          ... and {len(missing) - 10} more", file=sys.stderr)
        return 2

    with open(marker, "w", encoding="utf-8") as f:
        f.write(f"ready at {os.path.getmtime(target)}\n")
    print(f"[tiphd] OK: all {len(expected)} expected files present at {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
