"""Statescope batch runner: invoke run_Statescope_pair.py in a separate subprocess per pair (BLADE; outputs proportions)"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="Statescope", run_script="run_Statescope_pair.py", lang="python",
                 out_suffix="_proportions.csv",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
