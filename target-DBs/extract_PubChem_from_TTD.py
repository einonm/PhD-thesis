# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import os
import numpy as np
import pandas as pd
from IPython.display import display

# The Therapeutic Target Database (TTD) is available from  http://db.idrblab.net/ttd/full-data-download (last update: Jul 14th, 2019). We need to filter this data for targets of approved drugs, which the TTD marks as 'successful'.
#
# Read in the downloaded full TTD dataset:

df = pd.read_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/TTD/P1-01-TTD_download.txt",
    sep="\t",
    header=None,
    names=["tid", "item", "value"],
    index_col=0,
)

# Take a look at the format of the data
df.head()

# Examine all the items for a random example target
df[df.index == "T08123"]

# Data is listed per type, so look at all the types, pivot the table and remove unwanted ones.
tdd_cols = df["item"].unique()
display(len(tdd_cols))
display(tdd_cols)

# How many unique targets are there?
tdd_ids = df.index.unique()
display(len(tdd_ids))

# +
# Now pivot the table of targets, so each item type is a column of multiple values and each row a target
new_df = pd.DataFrame(columns=tdd_cols)

for index, row in df.iterrows():
    if index not in new_df.index.values or pd.isnull(new_df.loc[index][row["item"]]):
        new_val = row["value"]
    else:
        new_val = str(new_df.loc[index, row["item"]]) + ";" + row["value"]

    new_df.at[index, row["item"]] = new_val
# -

# Only look at successful targets (targeted by at least one approved drug), how many are there?
new_df = new_df[new_df["Type of target"] == "Successful target"]
display("There are " + str(len(new_df)) + " drugs with approved targets")

# For these successful targets, only keep the target ID (UniProt) and list of drugs targeting it
success_df = new_df.drop(
    [
        "Inhibitor (gating inhibitor)",
        "Regulator",
        "Inducer",
        "Immunomodulator (Immunostimulant)",
        "Blocker (channel blocker)",
        "Regulator (upregulator)",
        "Modulator (allosteric modulator)",
        "Cofactor",
        "Opener",
        "Immunomodulator",
        "Antagonist",
        "Blocker",
        "PDB Structure",
        "Breaker",
        "MetaCyc/BioCyc Pathway",
        "Suppressor",
        "Stablizer",
        "Binder (minor groove binder)",
        "Intercalator",
        "QSAR Model",
        "Inhibitor",
        "Enhancer",
        "Stimulator",
        "Wiki Pathway",
        "PathWhiz Pathway",
        "Activator",
        "Agonist",
        "Binder",
        "Reactome Pathway",
        "Modulator",
        "Pathway Interaction Database",
        "Function",
        "PANTHER Pathway",
        "NetPath Pathway",
        "KEGG Pathway",
        "Name",
        "Type of target",
        "Synonyms",
        "BioChemical Class",
        "EC Number",
        "Target Validation",
        "Disease",
    ],
    axis=1,
)
success_df.head()

# How many targets are there?
success_df.shape

# Drop those targets that do not have a UniProt ID listed. How many are left?
success_df.dropna(subset=["UniProt ID"], inplace=True)
success_df.shape

# +
# Take this list of drug sets per uniprot ID, and change it to a list of drugs against uniprot IDs
# Also tally how many proteins each drug has
drug_df = pd.DataFrame(columns=["uniprot_ids", "prot_count"])
drug_df.index.name = "drug"

for index, row in success_df.iterrows():
    drugs = row["Drug(s)"].split(";")
    prot = row["UniProt ID"]
    for drug in drugs:
        if drug not in drug_df.index.values:
            new_val = prot
            count = 1
        else:
            new_val = str(drug_df.loc[drug].uniprot_ids) + ";" + prot
            count = drug_df.loc[drug]["prot_count"] + 1

        drug_df.at[drug, "uniprot_ids"] = new_val
        drug_df.at[drug, "prot_count"] = count
# -

# What does the table of drugs vs protein targets look like?
drug_df.sort_values("prot_count", ascending=False).head()

# ### Convert UniProt IDs to EntrezGene

drug_df["uniprot_ids"] = [l.split(";") for l in drug_df["uniprot_ids"]]
len(drug_df)

# +
# Use idmapping table from UniProt found at http://www.uniprot.org/downloads/
idmap_full = pd.read_csv(
    "../data/drugbank/HUMAN_9606_idmapping_selected.tab",
    index_col=0,
    sep="\t",
    names=[
        "UniProtKB-AC",
        "UniProtKB-ID",
        "(EntrezGene)",
        "RefSeq",
        "GI",
        "PDB",
        "GO",
        "UniRef100",
        "UniRef90",
        "UniRef50",
        "UniParc",
        "PIR",
        "NCBI-taxon",
        "MIM",
        "UniGene",
        "PubMed",
        "EMBL",
        "EMBL-CDS",
        "Ensembl",
        "Ensembl_TRS",
        "Ensembl_PRO",
        "Add PubMed",
    ],
)

display(idmap_full.head())
idmap = idmap_full[["(EntrezGene)"]]

# Remove multiple gene IDs if they exist
# TODO - why are there > 1 genes listed for some proteins?
for index, row in idmap.iterrows():
    if str(row["(EntrezGene)"]).find(";") != -1:
        idmap.loc[index]["(EntrezGene)"] = row["(EntrezGene)"][
            : row["(EntrezGene)"].index(";")
        ]

# Remove MHC region genes?
# mhc_genes = pd.read_csv('../qtl/data/mhc_genes.txt', sep='\t')['(EntrezGene)']
# mhc_genes = mhc_genes.unique()
# idmap = idmap.loc[~idmap['(EntrezGene)'].isin(mhc_genes)]

idmap["(EntrezGene)"] = idmap["(EntrezGene)"].astype(str)
display(idmap.head())

# +
drug_genes_df = drug_df.copy()

drug_genes_df["genes"] = [[] for i in range(len(drug_genes_df))]
drug_genes_df["gene_count"] = 0

# +
for index, row in drug_genes_df.iterrows():
    for prot in row["uniprot_ids"]:
        if prot in idmap.index:
            gene = idmap.at[prot, "(EntrezGene)"]
            if gene not in row["genes"]:
                drug_genes_df.at[index, "genes"] = drug_genes_df.at[index, "genes"] + [
                    gene
                ]
                drug_genes_df.at[index, "gene_count"] = (
                    drug_genes_df.at[index, "gene_count"] + 1
                )


# del drug_genes_df['uniprot_ids']
drug_genes_df["diff"] = drug_genes_df["prot_count"] - drug_genes_df["gene_count"]
# -

len(
    drug_genes_df[drug_genes_df["gene_count"] > 1].sort_values(
        "gene_count", ascending=False
    )
)

drug_genes_df["genes"] = [sorted(l) for l in drug_genes_df["genes"]]
drug_genes_df["targets"] = [" ".join(g) for g in drug_genes_df["genes"]]
del drug_genes_df["genes"]

drug_genes_df.reset_index(inplace=True)
display(drug_genes_df.head())
display(len(drug_genes_df))

drug_genes_df["drug"] = [d.replace("'", "prime") for d in drug_genes_df["drug"]]

# +
del drug_genes_df["uniprot_ids"]
del drug_genes_df["prot_count"]
del drug_genes_df["gene_count"]
del drug_genes_df["diff"]

drug_genes_df.dropna(subset=["targets"], inplace=True)

display(drug_genes_df.head())
display(len(drug_genes_df))
# -


# Write this DF to a file
drug_genes_df.to_csv("../data/target-dbs/ttd_drug_genes.csv", sep="\t", index=None)

len(drug_genes_df.dropna())

# ### Convert drug names to PubChem CIDs with ATC codes.

# Convert drug names to IDs
ddf = pd.read_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/TTD/P1-02-TTD_crossmatching.txt",
    sep="\t",
    header=None,
    names=["did", "type", "item", "value"],
    index_col=0,
)

# Data is listed per type, so look at all the types, pivot the table and remove unwanted ones.
d_cols = ddf["item"].unique()
display(len(d_cols))
display(d_cols)
display(ddf.head())

d_ids = ddf.index.unique()
display(len(d_ids))
display(d_ids)

# +
new_ddf = pd.DataFrame(columns=d_cols)

for index, row in ddf.iterrows():
    if index not in new_ddf.index.values or pd.isnull(new_ddf.loc[index][row["item"]]):
        new_val = row["value"]
    else:
        new_val = str(new_ddf.loc[index, row["item"]]) + ";" + row["value"]

    new_ddf.at[index, row["item"]] = new_val
# -

display(new_ddf.head())
display(len(new_ddf))

# stats for all drugs in TTD
print("out of:- " + str(len(new_ddf)))
print("CAS Number: " + str(len(new_ddf[~(new_ddf["CAS Number"].isna())])))
print("Formular: " + str(len(new_ddf[~(new_ddf["Formular"].isna())])))
print("PubChem Compound: " + str(len(new_ddf[~(new_ddf["PubChem CID"].isna())])))
print("PubChem Substance: " + str(len(new_ddf[~(new_ddf["PubChem SID"].isna())])))
print("ChEBI: " + str(len(new_ddf[~(new_ddf["ChEBI ID"].isna())])))
print("SuperDrug ATC: " + str(len(new_ddf[~(new_ddf["SuperDrug ATC"].isna())])))

# Look at stats for drugs with approved drug targets
result = pd.merge(
    drug_genes_df, new_ddf, how="right", left_on="drug", right_on="DrugName"
)

display(result.head())
display(len(result))

# stats for all drugs in TTD
print("out of:- " + str(len(result)))
print("CAS Number: " + str(len(result[~(result["CAS Number"].isna())])))
print("Formular: " + str(len(result[~(result["Formular"].isna())])))
print("PubChem Compound: " + str(len(result[~(result["PubChem CID"].isna())])))
print("PubChem Substance: " + str(len(result[~(result["PubChem SID"].isna())])))
print("ChEBI: " + str(len(result[~(result["ChEBI ID"].isna())])))
print("SuperDrug ATC: " + str(len(result[~(result["SuperDrug ATC"].isna())])))

# Choose PubChem Compound ID, and tidy up
result.drop(
    [
        "drug",
        "DrugName",
        "CAS Number",
        "Formular",
        "PubChem SID",
        "ChEBI ID",
        "SuperDrug CAS",
    ],
    axis=1,
    inplace=True,
)

result.dropna(subset=["PubChem CID"], inplace=True)
result["PubChem CID"] = [x.strip().replace("CID ", "") for x in result["PubChem CID"]]
display(result.head())
display(len(result))

# Which of these are approved? Let's just use ATC codes for now...
result.dropna(subset=["SuperDrug ATC"], inplace=True)
result = result[["PubChem CID", "SuperDrug ATC", "targets"]]
result.rename(
    columns={"SuperDrug ATC": "ATC", "targets": "targets_entrez"}, inplace=True
)
result["ATC"] = [x.strip().replace(";", " ") for x in result["ATC"]]
display(result.head())
display(len(result))

result = result[result["targets_entrez"] != ""]
result["targets_entrez"] = result["targets_entrez"].astype(str)

len(result)

# Remove all duplicates, by joining results for ATC and entrez
result = (
    result.groupby("PubChem CID")
    .agg({"ATC": " ".join, "targets_entrez": " ".join})
    .reset_index()
)
display(result[result.duplicated(subset=["PubChem CID"], keep=False)])

display(result[result["PubChem CID"] == "5359476"])
display(result.shape)
display(len(result["PubChem CID"].unique()))

# de-dup ATC codes and entrez
result["ATC"] = [set(x) for x in result["ATC"].str.split()]
result["ATC"] = [x.remove("nan") if "nan" in x else x for x in result["ATC"]]
result["ATC"] = [" ".join(list(x)) if x != None else "nan" for x in result["ATC"]]
result.head()

result["targets_entrez"] = [set(x) for x in result["targets_entrez"].str.split()]
result["targets_entrez"] = [
    x.remove("nan") if "nan" in x else x for x in result["targets_entrez"]
]
result["targets_entrez"] = [
    " ".join(list(x)) if x != None else "nan" for x in result["targets_entrez"]
]
result.head()

result["targets_entrez"] = result["targets_entrez"].astype(str)
result = result[~result["PubChem CID"].str.contains(";")]
display(result.shape)
display(result.head())

# Write this DF to a file
result.to_csv("../data/target-dbs/ttd_ATC_targets.csv", sep="\t", index=None)

# ### Evidence of relationship checks (from individual drug results)
# To see what evidence this database provides for a drug-target interaction appearing in summary analysis results

# CACNB2 & Nifedipine?
# Nifedipine - C08CA05 C07FB03 C08GA01 C08CA55
result[result["PubChem CID"].str.contains("4485")]

result[result["ATC"].str.contains("C01EB10")]
result[result["ATC"].str.contains("Telapre")]
