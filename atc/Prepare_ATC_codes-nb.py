# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Prepare ATC code description table

import os
import pandas as pd
from IPython.display import display
import re

# The original list looks like this, a hierarchical tree, obtained from:
#
# <tt>http://www.genome.jp/kegg-bin/get_htext?br08303.keg</tt>

atc_path = "../data/atc/"
atc = pd.read_csv(atc_path + "br08303.keg", sep="\t", header=None, index_col=0)
atc.head(50)

# The original list/tree prefixes each level by a letter, A-E. Using this, we can create 5 files, one for each level. First, split them into 5 lists:

# +
atc_codes = [[] for i in range(6)]

code_labels = ["", "A", "B", "C", "D", "E"]

for index, row in atc.iterrows():
    for atc_level in range(1, 6):
        if index[0] == code_labels[atc_level]:
            atc_codes[atc_level].append(index[1:].strip())
# -

atc_codes

# Then write them to 5 files:

for i in range(1, 6):
    file = "level-" + str(i) + "-atc.txt"
    with open(atc_path + file, "w") as output:
        for line in atc_codes[i]:
            output.write(str(line) + "\n")
