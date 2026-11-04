"""Shared setup for the plot scripts (import it first).

DATA  - where the test outputs (transfer_*.csv, oh2d_*.csv, oh2d_shift/, oh2d_jitter/) are read from:
        $BC_DATA, default <package>/data (where the tests write their CSVs).
PLOTS - plot_scripts/plots; the script switches into it, so relative savefig names land there.
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("BC_DATA", os.path.join(_HERE, "..", "data"))
PLOTS = os.path.join(_HERE, "plots")
os.makedirs(PLOTS, exist_ok=True)
DATA = os.path.abspath(DATA)
os.chdir(PLOTS)
