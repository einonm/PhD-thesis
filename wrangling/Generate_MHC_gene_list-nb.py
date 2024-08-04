# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.7.0
# ---

# # Generate MHC gene list
#
# Create a list of all gene IDs that are in the MHC region - the Major Histocompatibility Complex block, which is on chr 6 between postions 26000000 and 33000000

import os
import pandas as pd
from IPython.display import display

# Use data/Ensembl.pos.txt to exclude those genes in the MHC region
df = pd.read_csv("/home/mpnme/data/Ensembl_pos.tsv", sep="\t")
df.head()

mhc_df = df[
    (df.chrom == "6")
    & (
        ((df.chromStart < 33000000) & (df.chromStart > 26000000))
        | ((df.chromEnd < 33000000) & (df.chromEnd > 26000000))
    )
]
display(mhc_df.head())
len(mhc_df)

mhc_df.to_csv("/home/mpnme/data/mhc_genes.txt", index=False, sep="\t")
