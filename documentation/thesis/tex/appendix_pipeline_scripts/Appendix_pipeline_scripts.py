# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import os
import sys
import pandas as pd
import numpy as np

# #!pip install scipy
from scipy import stats


# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt

# import qvalue
from IPython.display import display, HTML, Latex, display_latex, display_html

import warnings

warnings.filterwarnings("ignore")

import subprocess

# #!pip install seaborn
import seaborn as sns
import pylab
import scipy.stats as stats
from scipy.stats import chi2, norm

# %matplotlib inline
# #!pip install venn
import venn

import re

# Also generate tables in latex format, with formatting opts
pd.set_option("display.latex.repr", True)
pd.set_option("display.latex.longtable", True)
pd.set_option("display.latex.multirow", True)
pd.set_option("display.latex.multicolumn", True)
pd.set_option("display.max_colwidth", 2000)

# Switch when converting to PDF ('LaTeX') or ODT ('HTML')
# formatting = 'HTML'
formatting = "LaTeX"

if formatting == "HTML":
    print(
        "FORMATTED FOR REVIEW PURPOSES ONLY. SEE CORRESPONDING PDF FOR CORRECT LAYOUT."
    )
# -

dummy = HTML(
    """<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>."""
)

if formatting == "LaTeX":
    display(Latex(r"\doublespacing"))
    display(Latex(r"\setlength{\parskip}{5mm}"))
    display(Latex(r"\newpage"))

# +
#gitspec = subprocess.check_output(
#    ["git", "describe", "--all", "--dirty", "--long"]
#).strip()
#print("This notebook was generated from git revision " + str(gitspec, "utf-8"))
# -


# # Appendix: Pipeline scripts and notebooks
#
# All scripts appearing in these diagrams can be found in the \hyperref[supp:sf5]{Supplementary File 5} zip file, Figure \ref{fig:chapter2-pipeline_2} scripts in the target-DBs directory and Figure \ref{fig:chapter3-pipeline} in the magma_analysis directory, once unzipped.

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[htpb]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Pipeline used to generate DUGGIE and STITCH. \normalfont{Green parallelograms represent data sets, blue diamonds analysis steps.}}"
        )
    )
    display(Latex(r"\includegraphics[width=160mm]{../../img/Chapter2-pipeline_2}"))
    display(Latex(r"\label{fig:chapter2-pipeline_2}"))
    display(Latex(r"\end{figure}"))

# <img src="../../img/Chapter2-pipeline_2.png">

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[htpb]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Pipeline used to generate drug-disease association results. \normalfont{Green parallelograms represent data sets, blue diamonds analysis steps.}}"
        )
    )
    display(Latex(r"\includegraphics[width=160mm]{../../img/Chapter3-pipeline}"))
    display(Latex(r"\label{fig:chapter3-pipeline}"))
    display(Latex(r"\end{figure}"))

# <img src="../../img/Chapter3-pipeline.png">
