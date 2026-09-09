"""TiRank 批量调度：逐配对独立进程运行 run_TiRank_pair.py（真·基因对+神经网络，输出 Rank_Label）"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="TiRank", run_script="run_TiRank_pair.py", lang="python",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
