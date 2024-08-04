# ---
# jupyter:
#   jupytext:
#     notebook_metadata_filter: language_info
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
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
#     version: 3.8.17
# ---

# # Generate MHC gene list
#
# Create a list of all gene IDs that are in the MHC region - the Major Histocompatibility Complex block, which is on chr 6 between postions 26000000 and 33000000

import os
import pandas as pd
from IPython.display import display

# Use data/Ensembl.pos.txt to exclude those genes in the MHC region
mhc_df = pd.read_csv("../qtl/data/mhc_genes.txt", sep="\t")
mhc_df.head()

len(mhc_df)

duggie_df = pd.read_csv("../data/target-dbs/all_dgi_targets_atc_ensembl.csv",sep="\t", header=None, index_col=0)

duggie_df.head()

pattern = '|'.join(mhc_df['Ensembl'])

duggie_df[2] = duggie_df[1].str.replace(pattern, '')


def count_genes(string):
    return len(string.split())


duggie_df['orig_count'] = duggie_df[1].apply(lambda x: count_genes(x))
duggie_df['no_mhc_count'] = duggie_df[2].apply(lambda x: count_genes(x))

duggie_df[duggie_df['orig_count'] != duggie_df['no_mhc_count']]

# There are only 4 drugs affected by having MHC targets, affecting 1% or less of the targets of those drugs - so we can ignore MHC considerations.
