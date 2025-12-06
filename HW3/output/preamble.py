from pathlib import Path

import cvxpy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

from PIL import Image as image
from scipy.fft import dst, idst
from tqdm.auto import tqdm, trange

# from numpy import fft
from pyfftw.interfaces import numpy_fft as fft
import pyfftw

# --- 
HW3_DIR = Path().resolve()
OUTPUT_DIR = Path("./output")
OUTPUT_DIR.mkdir(exist_ok=True)

mm = 1 / 25.4  # mm in inches
PAPER_WIDTH_IN = (210 - 15 - 25) * mm

EXPORT_DPI = 600

sns.set_theme(context="paper", style="whitegrid")
gray_norm = mpl.colors.Normalize(vmin=0, vmax=1)

# --- 
# Based on: https://stackoverflow.com/a/29261252/5036246
def sci_notation(number, sig_fig=2, latex_style=True, omit_coeff_if_unity=True):
    ret_string = "{0:.{1:d}e}".format(number, sig_fig)
    a, b = ret_string.split("e")

    # Calculate factor part
    factor = ""
    joiner = ""
    if omit_coeff_if_unity and float(a) == 1.0:
        factor = ""
        joiner = ""
    else:
        factor = a
        joiner = "×"

    # Calculate exponent part
    b = int(b)
    exponent = ""
    if latex_style:
        exponent = f"10^{{{str(b)}}}"
    else:
        exponent = f"10^{str(b)}"

    return factor + joiner + exponent

# --- 
pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(60.0)  # A minute
