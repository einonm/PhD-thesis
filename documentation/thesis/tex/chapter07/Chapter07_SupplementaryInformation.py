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


# # Supplementary Information
#
# **Supplementary File 1: SF1-dti_lit_search_full.zip** - Results from a search of the PubMed website, using the search string ‘(((drug) AND (gene)) AND (interaction)) AND (database)’, constrained to the human species and freely available publications only. The 646 publication results were limited to the past 10 years from March 2019.

display(Latex(r"\label{supp:sf1}"))

# **Supplementary File 2: SF2-dti_lit_search_filtered.zip** - Remaining 24 publications resulting from a manual curation of supplementary file 1, containing only publications reporting on available drug-gene interaction databases.

display(Latex(r"\label{supp:sf2}"))

# **Supplementary File 3: SF3-DUGGIE_DGI_DB.zip** - The DUGGIE (DrUG-Gene IntEractions) drug-gene interaction database, with 1,323 approved drugs having unique drug target lists of 5 or more targets, listed by ATC code. There are 5,600 unique genes given as targets of these drugs in 64,312 unique interactions.

display(Latex(r"\label{supp:sf3}"))

# **Supplementary File 4: SF4-STITCH_DGI_DB.zip** - The STITCH drug-gene interaction database, the largest contributing database to DUGGIE, formatted and quality controlled in an identical manner to DUGGIE fro comparison purposes.

display(Latex(r"\label{supp:sf4}"))

# **Supplementary File 5: SF5-code_archive.zip** - Archive of bash scripts and python 3 jupyter notebook code used to conduct this study.

display(Latex(r"\label{supp:sf5}"))

# **Supplementary File 6: SF6-appendix-B_archive.zip** - Archive of scripts and results supporting the mini analysis in Appendix B.

display(Latex(r"\label{supp:sf6}"))

# **Supplementary File 7: SF7-results-data_archive.tar.gz** - Archive of all result data, sufficient to recreate this thesis document using the original thesis jupyter notebooks found in \ref{supp:sf5}.

display(Latex(r"\label{supp:sf7}"))
