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

# We want to transform the list of drug/target interactions into a list of pubchem CIDs vs Entrez IDs, with ATC codes if possible. Also include the evidence presented for the interaction (e.g. PMID).
#
# Taken from http://drugcentral.org/download

import numpy as np
import pandas as pd

df = pd.read_csv("../data/target-dbs/DrugCentral/drug.target.interaction.tsv", sep="\t")

display(df.shape)
display(len(df.columns))
display(df.columns)
df.head(10)

df2 = df[df["ORGANISM"] == "Homo sapiens"]
df3 = df2.drop(
    [
        "STRUCT_ID",
        "TARGET_NAME",
        "TARGET_CLASS",
        "ACCESSION",
        "SWISSPROT",
        "ACT_VALUE",
        "ACT_SOURCE",
        "ACT_TYPE",
        "ACT_UNIT",
        "ACT_COMMENT",
        "ACT_SOURCE",
        "RELATION",
        "MOA",
        "MOA_SOURCE",
        "ACT_SOURCE_URL",
        "MOA_SOURCE_URL",
        "ACTION_TYPE",
        "TDL",
        "ORGANISM",
    ],
    axis=1,
)
df3

single_atc_df = pd.read_csv("../data/atc/single_drugs_atc5.tsv", sep="\t")

# merge if possible
dc_atc_df = pd.merge(
    df3, single_atc_df, left_on="DRUG_NAME", right_on="description", how="left"
)
dc_atc_df

# How many drugs do not have ATC codes?
missing_atcs = dc_atc_df[dc_atc_df["atc_code"].isnull()]["DRUG_NAME"].unique()
len(missing_atcs)

# How many drugs?
len(dc_atc_df["atc_code"].unique())

result_df = dc_atc_df[~dc_atc_df["atc_code"].isnull()]

result_df.drop(["DRUG_NAME", "description"], axis=1, inplace=True)

result_df

# ### Convert HGNC to entrez using the Biomart conversion table

conversion_df = pd.read_csv(
    "/home/mpnme/data/phd/wrangling/Ensembl_prot_entrez.tsv", sep="\t"
)
conversion_df.head()

conv_df = pd.merge(
    dc_atc_df, conversion_df, left_on="GENE", right_on="HGNC symbol", how="left"
)
conv_df.drop(["Ensembl", "Ensembl_prot"], axis=1, inplace=True)
conv_df = conv_df.drop_duplicates()
conv_df.dropna(subset=["Entrez ID"], inplace=True)
conv_df["Entrez ID"] = conv_df["Entrez ID"].astype(int).astype(str)
conv_df.drop(["DRUG_NAME", "description", "GENE", "HGNC symbol"], axis=1, inplace=True)
display(conv_df.tail())
display(conv_df.shape)

conv_df = conv_df.dropna()
conv_df

# collect all entrez per ATC code and combine in one row
combined_df = (
    conv_df.groupby("atc_code")
    .apply(lambda x: " ".join(x["Entrez ID"]))
    .to_frame(name="Entrez ID")
)

combined_df.to_csv("../data/target-dbs/DrugCentral/drugcentral_targets.tsv", "\t")

combined_df
