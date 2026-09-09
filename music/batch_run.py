"""MuSiC 批量调度：逐配对独立进程运行 run_MuSiC_pair.R（输出 proportions）"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="MuSiC", run_script="run_MuSiC_pair.R", lang="R",
                 out_suffix="_proportions.csv",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
