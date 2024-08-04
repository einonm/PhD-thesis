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

# ## Examine functional - positional annotation differences ##
#
# We'd like to measure and visualise the differences between the functional and positional annotation of genes - each gene will be assigned a set of SNPs, which may be different between the annotations.

# +
import os
import sys
import numpy as np
import pandas as pd

# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns

from IPython.display import display

from pandas_profiling import ProfileReport

import warnings

warnings.filterwarnings("always")
# -

func_df = pd.read_csv(
    "../data/magma_analysis/combined.gene.annot-allQTLs", sep="\t", header=None
)

func_df.describe()

pos_df = pd.read_csv(
    "../data/wrangling/magma_g1000eur_NCBI37.3.prox.genes.annot.ensemble",
    sep="\t",
    header=None,
)

pos_df.describe()

intersect_df = pd.merge(pos_df, func_df, on=0, how="inner", suffixes=("_pos", "_func"))
intersect_df.head()

snp_lists_df = intersect_df.drop(["1_pos", "1_func"], axis=1)

snp_lists_df.head()

snp_lists_df["snps_pos"] = snp_lists_df["2_pos"].str.split()
snp_lists_df["snps_func"] = snp_lists_df["2_func"].str.split()
snp_lists_df.drop(["2_pos", "2_func"], axis=1, inplace=True)

# +
# So we can store the set as a dataframe item
snp_lists_df["intersect"] = snp_lists_df["snps_pos"].astype(object)

for index, row in snp_lists_df.iterrows():
    snp_lists_df.at[index, "snps_pos"] = set(snp_lists_df.at[index, "snps_pos"])
    snp_lists_df.at[index, "snps_func"] = set(snp_lists_df.at[index, "snps_func"])
    snp_lists_df.at[index, "len_pos"] = len(snp_lists_df.at[index, "snps_pos"])
    snp_lists_df.at[index, "len_func"] = len(snp_lists_df.at[index, "snps_func"])
    intersect = snp_lists_df.at[index, "snps_pos"].intersection(
        snp_lists_df.at[index, "snps_func"]
    )
    snp_lists_df.at[index, "intersect"] = set(intersect)
    snp_lists_df.at[index, "len_intersect"] = len(intersect)
    snp_lists_df.at[index, "len_pos_uniq"] = (
        snp_lists_df.at[index, "len_pos"] - snp_lists_df.at[index, "len_intersect"]
    )
    snp_lists_df.at[index, "len_func_uniq"] = (
        snp_lists_df.at[index, "len_func"] - snp_lists_df.at[index, "len_intersect"]
    )
snp_lists_df.head()
# -

snp_lists_df.describe()

# +
# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 1

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
# fig.subplots_adjust(hspace=0.5, wspace=0.3)
ax1 = plt.gca()
# plt.subplot(fig_height, fig_width, 1, sharex=ax1)
temp = plt.boxplot(
    [
        snp_lists_df["len_intersect"],
        snp_lists_df["len_pos"],
        snp_lists_df["len_func"],
        snp_lists_df["len_pos_uniq"],
        snp_lists_df["len_func_uniq"],
    ]
)

plt.legend()
plt.show()

# -
profile = ProfileReport(snp_lists_df)
profile.to_widgets()
