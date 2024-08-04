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

# ## Parse HD PED & MAP
# Here we want to convert the SNP ID column in chr/pos/allele format to an rs number.
#
# The order of the map file must be kept

import os
import numpy as np
import pandas as pd
from IPython.display import display

map_path = os.path.join("/home/mpnme/GWA4_12345.map.no_rs")
map_pairs = pd.read_csv(
    map_path, sep="\t", header=None, names=["chr", "snp_id", "dist", "pos"]
)
display(map_pairs.head())

# Convert chr and pos to an rs number using Kaviar

kaviar = pd.read_csv("../kaviar/data/kaviar_parsed.tsv", sep="\t")

# Convert chr column to string, to allow merge below to complete

# kaviar['chr'] = kaviar['chr'].astype(str)
# kaviar.columns = ['chr', 'snp_id', 'dist', 'bp']
display(kaviar.head())

joined_df2 = pd.merge_ordered(map_pairs, kaviar, on=["chr", "pos"], how="left")
# joined_df.columns=['chr', 'snp_id_x', 'dist', 'bp', 'delete']
display(joined_df2.head())

missing_list = joined_df2[joined_df2.snp_id_y.notnull()]
missing_list.drop(["chr", "dist", "pos"], axis=1, inplace=True)
missing_list.drop_duplicates(subset="snp_id_x", inplace=True)
display(missing_list.head())
missing_list.to_csv(
    os.path.join("/home/mpnme/", "missing_rs_ids.txt"),
    sep="\t",
    index=False,
    header=None,
)
## Now use plink to remove these from the bed/bim/fam
