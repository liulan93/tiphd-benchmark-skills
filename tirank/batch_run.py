"""TiRank batch runner: invoke run_TiRank_pair.py in a separate subprocess per pair (gene-pair network + neural net; outputs Rank_Label)"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="TiRank", run_script="run_TiRank_pair.py", lang="python",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
