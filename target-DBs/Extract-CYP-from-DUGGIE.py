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

duggie_df = pd.read_csv(
    "../data/target-dbs/all_dgi_targets_atc_ensembl.csv",
    sep="\t",
    header=None,
    names=['atc','genelist_str'])

duggie_df

cyp_list = pd.read_csv("../data/wrangling/CYP-list2.txt", header=None)

duggie_df['genelist_lst'] = duggie_df['genelist_str'].str.split(" ")

duggie_df.at[index, 'nocyp_lst'] = "str"
for index, row in duggie_df.iterrows():
    duggie_df.at[index, 'nocyp_lst'] = [gene for gene in row['genelist_lst'] if gene not in list(cyp_list[0])]

for index, row in duggie_df.iterrows():
    duggie_df.at[index, "nocyp_str"] = " ".join(
        [str(x) for x in sorted(list(set(duggie_df.at[index, 'nocyp_lst'])))]
    )

duggie_df.head()

duggie_df["old_count"] = [len(x) for x in duggie_df['genelist_lst']]
duggie_df["new_count"] = [len(x) for x in duggie_df['nocyp_lst']]

duggie_df['diff'] = duggie_df["old_count"] - duggie_df["new_count"]

len(duggie_df[duggie_df['diff'] != 0])

duggie_nocpy_df = duggie_df[['atc', 'nocyp_str']]

duggie_nocpy_df

duggie_nocpy_df.to_csv(
    "../data/target-dbs/duggie_nocyp.tsv",
    sep="\t",
    header=None,
    index=None,
)
