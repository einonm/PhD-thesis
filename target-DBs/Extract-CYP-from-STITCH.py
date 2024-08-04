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

import os
import numpy as np
import pandas as pd
from IPython.display import display, HTML, Latex, display_latex, display_html

stitch_df = pd.read_csv(
    "../data/target-dbs/stitch_dgi_targets_atc_ensembl.csv",
    sep="\t",
    header=None,
    names=['atc','genelist_str'])

stitch_df

cyp_list = pd.read_csv("../data/wrangling/CYP-list2.txt", header=None)

stitch_df['genelist_lst'] = stitch_df['genelist_str'].str.split(" ")

stitch_df.at[index, 'nocyp_lst'] = "str"
for index, row in stitch_df.iterrows():
    stitch_df.at[index, 'nocyp_lst'] = [gene for gene in row['genelist_lst'] if gene not in list(cyp_list[0])]

for index, row in stitch_df.iterrows():
    stitch_df.at[index, "nocyp_str"] = " ".join(
        [str(x) for x in sorted(list(set(stitch_df.at[index, 'nocyp_lst'])))]
    )

stitch_df.head()

stitch_df["old_count"] = [len(x) for x in stitch_df['genelist_lst']]
stitch_df["new_count"] = [len(x) for x in stitch_df['nocyp_lst']]

stitch_df['diff'] = stitch_df["old_count"] - stitch_df["new_count"]

len(stitch_df[stitch_df['diff'] != 0])

stitch_nocpy_df = stitch_df[['atc', 'nocyp_str']]

stitch_nocpy_df

stitch_nocpy_df.to_csv(
    "../data/target-dbs/stitch_nocyp.tsv",
    sep="\t",
    header=None,
    index=None,
)


