# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## List Duplicate ATC Codes ##
#
# We'd like to see what codes are duplicated, where multiple codes are used for the same chemical substance (as is common in the ATC, where there are several therapeutic effects of the drug). This only needs to be done for the final level, level 5.

import os
import pandas as pd
from IPython.display import display
import re

atc = pd.read_csv("../data/atc/level-5-atc.txt", sep="\ ")

atc
