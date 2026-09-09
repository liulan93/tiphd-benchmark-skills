"""
TiPhD Benchmark — shared Python configuration
---------------------------------------------
所有算法的 batch_run.py 与 Python 类 run_*.py 都 import 本文件，
因此数据路径、5 个癌种参数、金标准名称映射只在这里维护一份。

约定：
  - 单细胞数据【不子采样】，读取完整细胞。
  - 所有路径从本文件位置推导（相对路径），不写死绝对路径。
  - 运行脚本按 (cancer, sc_name, bulk_name) 三个命令行参数，单配对独立运行。
"""
import os
import sys

# ---- 数据与输出路径：环境变量优先，缺省 <cwd>/data 与 <cwd>/results ----
DATA_DIR = os.environ.get("TIPHD_DATA_DIR") or os.path.join(os.getcwd(), "data")
OUT_BASE = os.environ.get("TIPHD_OUT_DIR")  or os.path.join(os.getcwd(), "results")
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUT_BASE, exist_ok=True)

# CONFIG_DIR 保留以向后兼容（不再驱动路径）
CONFIG_DIR = os.path.dirname(os.path.abspath(__file__))   # skills/_toolkit/

# ---- Rscript / 特定 Python 环境（环境变量覆盖） ----
# Rscript: 缺省 PATH 中的 "Rscript"；Windows 用户可设 TIPHD_RSCRIPT 为绝对路径
RSCRIPT = os.environ.get("TIPHD_RSCRIPT", "Rscript")

# scPER 的 ADAE 步骤需 tensorflow；不同机器 conda env 路径不同，环境变量覆盖
# 留空字符串 = 用 sys.executable（即当前 conda env）
SCPER_PYTHON = os.environ.get("TIPHD_SCPER_PYTHON", "")

# scPER 的 MAGIC 插补（Rmagic/reticulate）需 magic-impute；环境变量覆盖
MAGIC_PYTHON = os.environ.get("TIPHD_MAGIC_PYTHON", "")

# ---- 5 个癌种统一配置（与 config.R 保持一致） ----
CANCERS = {
    "AML": dict(ct_col="CellType", tcol="OS_time", scol="OS_status",
                gold="gold_standard_AML.csv",
                nmap={"HSC-like": "HSC-like", "Prog-like": "Prog-like", "Mono-like": "Mono-like"},
                scs=["GSE116256"], bulks=["TCGA", "wave12", "wave34"]),
    "CRC": dict(ct_col="Cell_subtype", tcol="Time", scol="Status",
                gold="gold_standard_all_CRC.csv",
                nmap={"T follicular helper cells (Tfh)": "T follicular helper cells",
                      "IgA+ Plasma cells": "IgA+ Plasma", "SPP1+": "SPP1+A",
                      "CD8+ T cells": "CD8+ T cells", "Regulatory T cells": "Regulatory T cells",
                      "Myofibroblasts": "Myofibroblasts", "Stalk-like ECs": "Stalk-like ECs",
                      "T helper 17 cells": "T helper 17 cells", "Tip-like ECs": "Tip-like ECs"},
                scs=["GSE144735", "GSE132465"],
                bulks=["GSE14333", "GSE17536", "GSE33113", "GSE37892", "GSE39582"]),
    "HCC": dict(ct_col="subtype", tcol="OS_time", scol="OS_status",
                gold="gold_standard_all_HCC.csv",
                nmap={"CD4+ Tex": "CD4_Tex", "CD4+ Th": "CD4_Th", "Treg": "Treg",
                      "cDC1": "cDC1_CLEC9A", "cDC2": "cDC2_CD1C", "CD8+ Tex": "CD8_Tex",
                      "FOLR2+ TAMs (TAM1)": "TAM_FOLR2", "IgA+ B cells": "IgA+_Plasma_B_cells",
                      "LAMP3+ DCs": "LAMP3_DC", "NK cell": "NK_cells",
                      "PLVAP+ ECs": "EC_PLVAP",
                      "SPP1+ tumor-associated macrophages": "TAM_SPP1"},
                scs=["GSE149614", "GSE151530"],
                bulks=["GSE116174", "GSE14520", "GSE76427"]),
    "LUAD": dict(ct_col="Cell_subtype", tcol="OS_time", scol="OS_status",
                 gold="gold_standard_all_LUDA.csv",
                 nmap={"CD8+ T cell": "CD8+ T cell", "cDC1": "cDC1", "cDC2": "cDC2_CD1C",
                       "Follicular B cells": "Follicular B cells",
                       "follicular helper T cell": "follicular helper T cell",
                       "IgG+ Plasma cells": "IgG+ Plasma B cells", "Macro_ISG15": "Macro_ISG15",
                       "cDC_LAMP3 (Mature DCs)": "cDC_LAMP3", "myCAF": "myCAF",
                       "Macro_PPARG": "Macro_PPARG", "Macro_SPP1": "Macro_SPP1", "Treg": "Treg"},
                 scs=["GSE127465", "GSE131907", "GSE148071"],
                 bulks=["GSE31210", "GSE3141", "GSE37745", "GSE68465", "GSE72094"]),
    "GC": dict(ct_col="annotation", tcol="Time", scol="Status",
               gold="gold_standard_all_GC.csv",
               nmap={"eCAF": "eCAF", "DC": "DC",
                     "exhausted CD8+ T cells": "CD8-LAYN-exhausted",
                     "iCAF": "iCAF", "FOXP3+ treg cells": "FOXP3-Treg",
                     "Treg": "FOXP3-Treg", "plasma cell": "Plasma.cells", "CAFs": "eCAF"},
               scs=["GSE183904"],
               bulks=["GSE15459", "GSE26253", "GSE26899", "GSE26901", "GSE28541",
                      "GSE29272", "GSE34942", "GSE57303", "GSE66229", "GSETCGA"]),
}

# 全部 (cancer, sc, bulk) 配对
PAIRS = [(cn, sc, bk)
         for cn, c in CANCERS.items()
         for sc in c["scs"] for bk in c["bulks"]]


def run_batch(algo, run_script, lang="R", out_suffix=".csv", timeout=None, base_dir=None):
    """
    通用批量调度器：逐配对独立起进程运行 run 脚本，跳过已有输出，
    单个失败不中断全局，最终写 _batch_log.txt。

    algo        算法名（输出文件前缀）
    run_script  run 脚本名（位于本算法文件夹内）
    lang        'R' 用 Rscript，'python' 用当前解释器
    out_suffix  输出文件后缀：A 类为 '.csv'，B 类为 '_proportions.csv'
    timeout     单配对超时秒数。None 时从环境变量 TIPHD_PAIR_TIMEOUT 读取，
                默认 7200（2小时）。GC（~137K细胞）的深度学习配对可能需要 30-120 分钟。
    base_dir    算法文件夹路径（通常由 batch_run.py 传入自身所在目录）
    """
    import subprocess
    import time

    if timeout is None:
        timeout = int(os.environ.get("TIPHD_PAIR_TIMEOUT", "7200"))

    if base_dir is None:
        base_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    script = os.path.join(base_dir, run_script)

    log_path = os.path.join(base_dir, "_batch_log.txt")
    logf = open(log_path, "w", encoding="utf-8")

    def log(msg):
        print(msg)
        logf.write(msg + "\n")
        logf.flush()

    total = len(PAIRS)
    failed, skipped, done = [], 0, 0
    log(f"=== {algo} batch run: {total} pairs ===")
    log(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    for i, (cancer, sc, bulk) in enumerate(PAIRS):
        label = f"[{i+1}/{total}] {cancer}|{sc}|{bulk}"
        out_csv = os.path.join(base_dir, cancer, f"{algo}_{bulk}_{sc}{out_suffix}")
        if os.path.exists(out_csv):
            log(f"{label} SKIP (exists)")
            skipped += 1
            continue

        log(f"{label} ...")
        cmd = ([RSCRIPT, script, cancer, sc, bulk] if lang == "R"
               else [sys.executable, script, cancer, sc, bulk])
        # Limit threads to prevent pthread_create errors on high-core machines
        run_env = os.environ.copy()
        run_env["OMP_NUM_THREADS"] = "1"
        run_env["OPENBLAS_NUM_THREADS"] = "1"
        run_env["MKL_NUM_THREADS"] = "1"
        run_env["NUMEXPR_NUM_THREADS"] = "1"
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=timeout,
                               env=run_env,
                               encoding="utf-8", errors="replace")
            tail = [l for l in r.stdout.strip().split("\n")[-5:] if l.strip()]
            for line in tail:
                log(f"  {line.strip()}")
            if r.returncode != 0:
                log(f"  FAIL rc={r.returncode}")
                for line in [l for l in r.stderr.strip().split("\n") if l.strip()][-5:]:
                    log(f"  ERR: {line.strip()}")
                failed.append(label)
            else:
                log("  DONE")
                done += 1
        except subprocess.TimeoutExpired:
            log(f"  TIMEOUT (>{timeout}s)")
            failed.append(label)
        except Exception as e:
            log(f"  EXCEPTION: {e}")
            failed.append(label)

    log("\n=== Summary ===")
    log(f"Total: {total}, Done: {done}, Skipped: {skipped}, Failed: {len(failed)}")
    if failed:
        log("Failed pairs:")
        for f in failed:
            log(f"  {f}")
    log(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logf.close()
