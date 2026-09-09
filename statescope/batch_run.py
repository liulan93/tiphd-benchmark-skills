"""Statescope 批量调度：逐配对独立进程运行 run_Statescope_pair.py（真·BLADE，输出 proportions）"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="Statescope", run_script="run_Statescope_pair.py", lang="python",
                 out_suffix="_proportions.csv",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
