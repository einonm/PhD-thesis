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

# # Extract ChEMBL from T3DB
#
# This notebook takes the full drugbank database and extracts each drug that has a protein target, the protein target and the drug ChEMBL code

import os
import csv
import gzip
import collections
import re
import io
import json
import xml.etree.ElementTree as ET

# #!pip install requests
import requests
import pandas as pd
from IPython.display import display

# Load in the T3DB XML (All Toxin Records (with Toxin-Target Mechanisms of Action and References) from here http://www.t3db.ca/downloads)

# +
data_path = "../data/target-dbs/T3DB/"

tree = ET.parse(data_path + "toxins.xml")
root = tree.getroot()
# -

rows = list()
for i, drug in enumerate(root):
    row = collections.OrderedDict()
    row["pubchem_cid"] = drug.findtext("pubchem_compound_id")
    row["uniprot"] = drug.findtext("targets/target/uniprot_id")
    row["references"] = [
        target.text
        for target in drug.findall(
            "targets/target[uniprot_id='"
            + str(row["uniprot"])
            + "']/references/reference/pubmed_id"
        )
    ]
    for index, ref in enumerate(row["references"]):
        if ref == None:
            row["references"][index] = "No pubmed, see T3DB"

    if row["references"] == []:
        row["references"] = ["No pubmed, see T3DB"]
    rows.append(row)

columns = ["pubchem_cid", "uniprot", "references"]
t3db_df = pd.DataFrame.from_dict(rows)[columns]

t3db_df  # .head()

# Load in a mapping of UniProt to gene ID

# +
# convert uniprot IDs to Ensembl
data_path = "../data/drugbank/"
idmap_path = data_path + "HUMAN_9606_idmapping_selected.tab"
idmap_full = pd.read_csv(
    idmap_path,
    index_col=0,
    sep="\t",
    names=[
        "UniProtKB-AC",
        "UniProtKB-ID",
        "EntrezGene",
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
        "Additional PubMed",
    ],
)

idmap = idmap_full[["EntrezGene", "Ensembl"]]

# Remove multiple gene IDs if they exist
# TODO - why are there > 1 genes listed for some proteins?
for index, row in idmap.iterrows():
    if str(row["EntrezGene"]).find(";") != -1:
        idmap.at[index, "EntrezGene"] = row["EntrezGene"][
            : row["EntrezGene"].index(";")
        ]
    if str(row["Ensembl"]).find(";") != -1:
        idmap.at[index, "Ensembl"] = row["Ensembl"][: row["Ensembl"].index(";")]

# Remove MHC region genes?
# mhc_genes = pd.read_csv('../qtl/data/mhc_genes.txt', sep='\t')['Ensembl']
# mhc_genes = mhc_genes.unique()
# idmap = idmap.loc[~idmap['Ensembl'].isin(mhc_genes)]

display(idmap.head())
# -

result_df = pd.merge(t3db_df, idmap, left_on="uniprot", right_index=True)
result_df.drop(["uniprot"], axis=1, inplace=True)
result_df.head()

result_df.to_csv(
    "../data/target-dbs/T3DB/t3db_evidence.tsv",
    sep="\t",
    index=None,
    columns=["pubchem_cid", "Ensembl", "references"],
)

# Now transform the result to one row per drug, listing it's targets


# +
result_df.dropna(subset=["EntrezGene"], inplace=True)

drug_df = pd.DataFrame(columns=["gene_list"])
drug_df.index.name = "pubchem_cid"

for index, row in result_df.iterrows():
    drug = row["pubchem_cid"]
    if drug not in drug_df.index.values:
        drug_df.at[drug, "gene_list"] = [int(row["EntrezGene"])]
    else:
        if row["EntrezGene"] not in drug_df.at[drug, "gene_list"]:
            drug_df.at[drug, "gene_list"].append((int(row["EntrezGene"])))

drug_df["gene_list"] = [sorted(set(l)) for l in drug_df["gene_list"]]

# -

drug_df["gene_list"] = [" ".join(str(l) for l in x) for x in drug_df["gene_list"]]
drug_df.dropna(inplace=True)
drug_df = drug_df[drug_df.index != ""]

drug_df.head(n=100)

drug_df.to_csv("../data/target-dbs/T3DB/t3db_targets.tsv", sep="\t")
