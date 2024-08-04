# ---
# jupyter:
#   jupytext:
#     notebook_metadata_filter: language_info
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
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
#     version: 3.8.5
# ---

# # Remove MHC genes from the proximity annotation geneset

import os
import pandas as pd
from IPython.display import display

# +
# Load the annotation file
# Note that this file must be changed to have tabs as the first two delimiters
# In bash: cat <file> | sed 's/\t/ /g' | sed 's/\([^\s]*\)\s/\1\t/' | sed 's/\s/\t/'/'
prox_annot = pd.read_csv(
    "../data/wrangling/magma_g1000eur_NCBI37.3.prox.genes.annot",
    sep="\t",
    comment="#",
    header=None,
    index_col=0,
)

print(len(prox_annot))
display(prox_annot.head())

# +
# Load the MHC genes and remove them from the annotation list
mhc_genes = pd.read_csv("../qtl/data/mhc_genes.txt", sep="\t")["Ensembl"]
mhc_genes = mhc_genes.unique()
prox_annot = prox_annot.loc[~prox_annot.index.isin(mhc_genes)]

print(len(prox_annot))
display(prox_annot.head())

prox_annot.to_csv(
    "../data/wrangling/magma_g1000eur_NCBI37.3.prox.genes.annot.no_mhc",
    sep="\t",
    header=False,
)
