"""scSurv 批量调度：逐配对独立进程运行 run_scSurv_pair.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="scSurv", run_script="run_scSurv_pair.py", lang="python",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
