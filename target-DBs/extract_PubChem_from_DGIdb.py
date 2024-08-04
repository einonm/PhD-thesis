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

import numpy as np
import pandas as pd

# We want to transform the list of drug/target interactions into a list of pubchem CIDs vs Entrez IDs, with ATC codes. Data downloaded from http://www.dgidb.org/downloads (Last updated 2018-01-25).

df = pd.read_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/dgidb/dgidb_interactions.tsv",
    sep="\t",
)

display(df.shape)
df.head(100)

df["interaction_claim_source"].unique()

df = df[
    [
        "gene_name",
        "entrez_id",
        "drug_name",
        "drug_chembl_id",
        "interaction_claim_source",
        "PMIDs",
    ]
]
df.dropna(subset=["entrez_id", "drug_chembl_id"], inplace=True)

display(len(df))

len(df["drug_chembl_id"].unique())

# +
drug_df = pd.DataFrame(columns=["gene_list", "gene_count"])
drug_df.index.name = "drug_chembl_id"

for index, row in df.iterrows():
    drug = row["drug_chembl_id"]
    if drug not in drug_df.index.values:
        drug_df.at[drug, "gene_list"] = [int(row["entrez_id"])]
    else:
        if row["entrez_id"] not in drug_df.at[drug, "gene_list"]:
            drug_df.at[drug, "gene_list"].append((int(row["entrez_id"])))

drug_df["gene_list"] = [sorted(l) for l in drug_df["gene_list"]]

display(drug_df.head())
# -

drug_df["gene_count"] = drug_df["gene_list"].str.len()
drug_df["gene_list"] = [" ".join(str(l) for l in x) for x in drug_df["gene_list"]]

display(len(drug_df.index.unique()))
display(len(drug_df[drug_df["gene_count"] >= 5]))
display(drug_df.head())

chembl_only = drug_df.index
chembl_series = chembl_only.to_series()
chembl_series.to_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/dgidb/chembl.list",
    header=None,
    index=None,
)

# **Submit list as a job to https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi. Convert from Chembl to pubchem CID.**

pubchem_cid = pd.read_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/dgidb/2777745265049466990.txt",
    sep="\t",
    header=None,
    names=["drug_chembl_id", "pubchem_cid"],
)

display(len(pubchem_cid))
pubchem_cid.dropna(inplace=True)
display(len(pubchem_cid))
pubchem_cid.head()

# merge IDs back


result_df = pd.merge(drug_df, pubchem_cid, how="left", on="drug_chembl_id")


display(len(result_df))
result_df.dropna(subset=["pubchem_cid"], inplace=True)
display(len(result_df))
result_df["pubchem_cid"] = result_df["pubchem_cid"].astype(int)
result_df.head()
display(len(result_df[result_df["gene_count"] >= 5]))

result_df.drop(["drug_chembl_id", "gene_count"], axis=1)
result2_df = result_df[["pubchem_cid", "gene_list"]]

result2_df.head()

result2_df.to_csv("../data/target-dbs/dgidb_targets.csv", index=None, sep="\t")
