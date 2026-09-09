"""scAB 批量调度：逐配对独立进程运行 run_scAB_pair.R"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="scAB", run_script="run_scAB_pair.R", lang="R",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
