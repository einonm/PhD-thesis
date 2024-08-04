# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# ## Generate a functional and positional hybrid data set
#
# We want to use all functional gene sets, and add any whole gene sets from positional data that are not present in the functional data. This notebook can do this for all QTL datasets - GTeX _ methylation/CNC, eQTLgen and psychEncode.

import os
import sys
import pandas as pd
import numpy as np
from IPython.display import display

# +
# get list of qtls we want to combine.
# Currently in use, uncomment one set to run.

#combined_prefix = "all-eqtl-brain"
#combined_prefix = "all-brain"
#combined_prefix = "all-eqtl-blood"
#combined_prefix = "all-eqtl-tissues"
#combined_prefix = "all-tissues"

# Previously used data sets

combined_prefix = "combined"
#combined_prefix = "new_combined"
#combined_prefix = "everything_combined"

# load in data
func_file = "../data/magma_analysis/" + combined_prefix + ".gene.annot-QTLs"

func_df = pd.read_csv(
    func_file, sep="\t", header=None, index_col=0, names=["gene", "pos", "snps"]
).sort_index()
display(func_df.head())
# -

# Load in positional data
pos_file = "../data/wrangling/magma_g1000eur_NCBI37.3.prox.genes.annot.ensemble"
pos_df = pd.read_csv(
    pos_file,
    sep="\t",
    header=None,
    index_col=0,
    names=["gene", "pos", "snps"],
    comment="#",
).sort_index()
display(pos_df.head())

# ## 'Missing gene' hybrid strategy
#
# This approach only includes positional genes that are missing from the functional set

# select genes only in the positional list
only_pos_df = (
    pd.merge(
        pos_df, func_df, left_index=True, right_index=True, how="outer", indicator=True
    )
    .query('_merge == "left_only"')
    .drop(["_merge"], axis=1)
    .sort_index()
)
only_pos_df.drop(["pos_y", "snps_y"], axis=1, inplace=True)
only_pos_df.columns = ["pos", "snps"]
display(only_pos_df.head())

# add only_pos_df to func_df
hybrid_df = pd.concat([only_pos_df, func_df]).sort_index()
display(hybrid_df.head())

hybrid_df.to_csv(
    "../data/magma_analysis/" + combined_prefix + "_hybrid.gene.annot-QTLs", sep="\t", header=False
)

# ## Results QC

hybrid_df

hybrid_df[hybrid_df.index.duplicated() == True]

hybrid_df.isnull().any()

# Find average & std dev of snps in annotation data set
hybrid_df['counts'] = hybrid_df['snps'].str.split().str.len()
hybrid_df.describe()

# ## 'Combined gene' hybrid approach
#
# We want to take all functional and all positional data sets and amalgamate the two - so each gene has SNPs from both the functional and positional data sets.

display(pos_df.head())
display(len(pos_df))
display(func_df.head())  # GRCh37
display(len(func_df))

combined_df = pd.merge(
    pos_df,
    func_df,
    how="outer",
    left_index=True,
    right_index=True,
    suffixes=["_pos", "_qtl"],
)

# +
# positional pos_x is out by 10k as it includes window applied by MAGMA
# for index, row in combined_df.iterrows():
#    if row['pos_x'] != row['pos_y']:
#        display(row['pos_x'] + "  " + row['pos_y'])
# -

combined_df["pos_qtl"].fillna(combined_df["pos_pos"], inplace=True)

combined_df.fillna("", inplace=True)

combined_df["snps"] = combined_df["snps_pos"] + " " + combined_df["snps_qtl"]

combined_df["dedup_snps"] = [
    " ".join(sorted(set(x.split()))) for x in combined_df["snps"]
]

combined_df.drop(["pos_pos", "snps_pos", "snps_qtl", "snps"], axis=1, inplace=True)

combined_df.to_csv(
    "../data/magma_analysis/" + combined_prefix + "_hybridboth.gene.annot-QTLs",
    sep="\t",
    header=False,
)

# ## Results QC

combined_df

combined_df[combined_df.index.duplicated() == True]

combined_df.isnull().any()

# Find average & std dev of snps in annotation data set
combined_df['counts'] = combined_df['dedup_snps'].str.split().str.len()
combined_df.describe()


