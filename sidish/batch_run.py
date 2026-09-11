"""SIDISH batch runner: invoke run_SIDISH_pair.py in a separate subprocess per pair"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_toolkit"))
import config

config.run_batch(algo="SIDISH", run_script="run_SIDISH_pair.py", lang="python",
                 base_dir=os.path.dirname(os.path.abspath(__file__)))
