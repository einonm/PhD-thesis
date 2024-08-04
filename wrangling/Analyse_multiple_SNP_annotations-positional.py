# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Analyse multiple SNP annotations
#
# When a SNP is assigned to a gene in an annotation process, it may get assigned to more than one gene. This will reduce the power.
#
# How much of this is going on with the QTL and proximity datasets we are using?

# +
import os
import pandas as pd
from IPython.display import display

import numpy as np

# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns

# -

# Proximity dataset
prox_df = pd.read_csv(
    "../data/wrangling/magma_gene_atlas.prox.genes.annot.ensemble",
    sep="\t",
    names=["Ensembl", "pos", "snplist"],
)
prox_df.drop(["pos"], axis=1, inplace=True)
prox_df.set_index("Ensembl", inplace=True)

prox_df["snplist"] = prox_df["snplist"].str.split(" ")

# +
# pivot table to one of SNPs vs genelist, from genes vs snplist
snpgenes_df = pd.DataFrame(columns=["snp", "genelist"])
# snpgenes_df.set_index('snp', inplace=True)

for index, row in prox_df.iterrows():
    for snp in row["snplist"]:
        snpgenes_df = snpgenes_df.append(
            {"snp": snp, "genelist": index}, ignore_index=True
        )
# -

groups = snpgenes_df.groupby(["snp"])

sizes = []
for name, group in groups:
    sizes.append(len(group))

plt.hist(sizes)

d = {x: sizes.count(x) for x in sizes}

counts = {
    1: 3193857,
    2: 369281,
    3: 42714,
    5: 2561,
    4: 14188,
    7: 3793,
    21: 143,
    16: 1737,
    6: 4069,
    10: 1003,
    14: 649,
    15: 320,
    8: 4786,
    12: 235,
    20: 578,
    19: 330,
    35: 53,
    9: 874,
    11: 172,
    30: 220,
    23: 156,
    18: 416,
    22: 188,
    24: 412,
    13: 269,
    32: 55,
    29: 47,
    40: 34,
    28: 187,
    42: 23,
    38: 9,
    17: 33,
    26: 128,
    31: 46,
    39: 19,
    33: 25,
    25: 23,
    36: 3,
    37: 7,
}

df = pd.DataFrame(list(counts.items()), columns=["genes-per-snp", "snpcount"])
df.sort_values("genes-per-snp", inplace=True)
df.set_index(["genes-per-snp"], inplace=True)

# How many snps are there in total?
totalsnps = df["snpcount"].sum()
single_gene_snps = df.at[1, "snpcount"]
print(totalsnps)
multiple_snp_genes = totalsnps - single_gene_snps
print(multiple_snp_genes)

# percentage of SNPs assigned to > 1 gene?
(multiple_snp_genes / totalsnps) * 100

df["log10_snpcount"] = np.log10(df["snpcount"])

df

dummy = df["snpcount"].plot.bar(figsize=(12, 10))

dummy = df["log10_snpcount"].plot.bar(figsize=(12, 12))
