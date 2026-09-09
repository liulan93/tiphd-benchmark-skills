"""scPER 批量调度：逐配对独立进程运行 run_scPER_pair.py（真·ADAE+xgboost，输出 proportions）"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="scPER", run_script="run_scPER_pair.py", lang="python",
                 out_suffix="_proportions.csv",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
