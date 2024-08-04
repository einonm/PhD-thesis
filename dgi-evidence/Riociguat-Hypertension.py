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

# We want to explore the drug-gene interaction evidence behind the Riociguat-Hypertension association, ATC code C02KX05
#
# * ENSG00000164116    GUCY1A3  	9.152
# * ENSG00000061918    GUCY1B3    6.8972
# * ENSG00000152402    GUCY1A2    3.778
# * ENSG00000188153    COL4A5     nan

targets = ["ENSG00000164116", "ENSG00000061918", "ENSG00000152402", "ENSG00000188153"]
ATC = "C02KX05"
pubchemCID = 11304743

# ### Evidence from DrugBank

db_df = pd.read_csv("../data/target-dbs/drugbank_evidence.csv", sep="\t", index_col=2)
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

gtopdb_df = pd.read_csv("../data/target-dbs/gtopdb_dgi_evidence.csv", sep="\t")
gtopdb_df[gtopdb_df["pubchem_cid"] == pubchemCID]


# ### Evidence from T3DB

t3db_df = pd.read_csv("../data/target-dbs/T3DB/t3db_evidence.tsv", sep="\t")
t3db_df[t3db_df["pubchem_cid"] == pubchemCID]
