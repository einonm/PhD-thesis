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

import os
import numpy as np
import pandas as pd
from IPython.display import display

# The DsigDB database is available from http://dsigdb.tanlab.org/Downloads/DSigDB_All_detailed.txt (last update: May 2015, release 1). It's split into 4 parts - D1-D4. We want D1 and D2.

dsigdb_df = pd.read_csv("../data/target-dbs/DsigDB/DSigDB_All_detailed.txt", sep="\t")

dsigdb_df["Source"].unique()

# +
# We want D1 - 'D1 PubChem', 'D1 ChEMBL' and D2 - 'Roche', 'GSK', 'FDA', 'Kinome Scan', 'RBC', 'MRC', 'LINCS'.
accepted_sources = [
    "D1 PubChem",
    "D1 ChEMBL",
    "Roche",
    "GSK",
    "FDA",
    "Kinome Scan",
    "RBC",
    "MRC",
    "LINCS",
]

dsigdb_sources_df = dsigdb_df[dsigdb_df["Source"].isin(accepted_sources)]
dsigdb_sources_df["Drug"] = dsigdb_sources_df["Drug"].str.lower()
# -

# We're only interested in drugs and genes, so remove type and source then get rid of duplicates
dsigdb_filtered_df = dsigdb_sources_df.drop(["Type", "Source"], axis=1)
dsigdb_filtered_df = dsigdb_filtered_df.drop_duplicates()

single_atc_df = pd.read_csv("../data/atc/single_drugs_atc5.tsv", sep="\t")

# merge if possible
dsigdb_atc_df = pd.merge(
    dsigdb_filtered_df,
    single_atc_df,
    left_on="Drug",
    right_on="description",
    how="left",
)
dsigdb_atc_df

# How many drugs do not have ATC codes?
missing_atcs = dsigdb_atc_df[dsigdb_atc_df["atc_code"].isnull()]["Drug"].unique()
len(missing_atcs)

# How many drugs?
len(dsigdb_atc_df["atc_code"].unique())

# This should be less than the measurements above!
missing_atcs = dsigdb_atc_df[dsigdb_atc_df["atc_code"].isnull()]["Drug"].unique()
display(missing_atcs)
# How many drugs do not have ATC codes?
len(missing_atcs)

dsigdb_atc_df = dsigdb_atc_df.dropna()

dsigdb_atc_df

# How many drugs have we got now?
len(dsigdb_atc_df["atc_code"].unique())

# ### Convert HGNC to entrez using the Biomart conversion table

conversion_df = pd.read_csv(
    "/home/mpnme/data/phd/wrangling/Ensembl_prot_entrez.tsv", sep="\t"
)
conversion_df.head()

conv_df = pd.merge(
    dsigdb_atc_df, conversion_df, left_on="Gene", right_on="HGNC symbol", how="left"
)
conv_df.drop(["Ensembl", "Ensembl_prot"], axis=1, inplace=True)
conv_df = conv_df.drop_duplicates()
conv_df.dropna(subset=["Entrez ID"], inplace=True)
conv_df["Entrez ID"] = conv_df["Entrez ID"].astype(int).astype(str)
conv_df.drop(["description", "Drug", "Gene", "HGNC symbol"], axis=1, inplace=True)
display(conv_df.tail())
display(conv_df.shape)

conv_df = conv_df.dropna()
display(conv_df.shape)

# collect all entrez per ATC code and combine in one row
combined_df = (
    conv_df.groupby("atc_code")
    .apply(lambda x: " ".join(x["Entrez ID"]))
    .to_frame(name="Entrez ID")
)

combined_df

combined_df.to_csv("../data/target-dbs/DsigDB/dsigdb_targets.tsv", sep="\t")

combined_df.at["L01XE22", "Entrez ID"]
