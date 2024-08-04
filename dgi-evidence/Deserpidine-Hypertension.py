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

import os
import numpy as np
import pandas as pd
from IPython.display import display
import webbrowser

# We want to explore the drug-gene interaction evidence behind the Deserpidine-Hypertension association, ATC code C02AA05
#
# * ENSG00000159640 ACE 6.6907
# * ENSG00000197635 DPP4 4.4732
# * ENSG00000142002 DPP9 3.4908
# * ENSG00000165646 SLC18A2 1.6812
# * ENSG00000074603 DPP8 1.4693

targets = [
    "ENSG00000159640",
    "ENSG00000197635",
    "ENSG00000142002",
    "ENSG00000165646",
    "ENSG00000074603",
]
ATC = "C02AA05"
pubchemCID = 8550

# ### Evidence from DrugBank

# Check drugbank
db_df = pd.read_csv("../data/target-dbs/drugbank_evidence.csv", sep="\t", index_col=2)

db_df.head()

db_df[db_df["atc_codes"].str.contains(ATC)]

# +
for target in targets:
    if target in db_df[db_df["atc_codes"].str.contains(ATC)].index:
        results = db_df.loc[target]


bad_chars = ["[", "]", '"', "'", " "]
for index, row in results.iterrows():
    if ATC in row["atc_codes"]:
        string = row.refs
        for s in bad_chars:
            string = string.replace(s, "")

lst = string.split(",")

display(lst)

for i in lst:
    webbrowser.open("https://www.ncbi.nlm.nih.gov/pubmed/" + i)
# -

# ### Evidence from GtoPdb

# Check GtoPdb
gtopdb_df = pd.read_csv("../data/target-dbs/gtopdb_dgi_evidence.csv", sep="\t")
gtopdb_df.head()


gtopdb_df[gtopdb_df["pubchem_cid"] == pubchemCID]

# ### Evidence from T3DB

t3db_df = pd.read_csv("../data/target-dbs/T3DB/t3db_evidence.tsv", sep="\t")
t3db_df.head()

t3db_df[t3db_df["pubchem_cid"] == pubchemCID]
