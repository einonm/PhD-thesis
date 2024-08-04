# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
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

# + [markdown] {"toc": true}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Extract-ChEMBL-from-DrugBank" data-toc-modified-id="Extract-ChEMBL-from-DrugBank-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Extract ChEMBL from DrugBank</a></span></li></ul></div>
# -

# # Extract ChEMBL from DrugBank
#
# This notebook takes the full drugbank database and extracts each drug that has a protein target, the protein target and the drug ChEMBL code.

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

# Load in the drugbank XML (version 5.1.6, obtained from https://www.drugbank.ca/releases/latest, last updated 2020-04-22):

# +
data_path = "../data/drugbank/"

tree = ET.parse(data_path + "full-drugbank-database-5.1.6.xml")
root = tree.getroot()
# -

# Now extract a list of table rows, with drugbank ID, name, ATC codes and protein target information:

# +
ns = "{http://www.drugbank.ca}"

rows = list()
for i, drug in enumerate(root):
    row = collections.OrderedDict()
    assert drug.tag == ns + "drug"
    row["drugbank_id"] = drug.findtext(ns + "drugbank-id[@primary='true']")
    row["name"] = drug.findtext(ns + "name")
    row["type"] = drug.attrib["type"]
    row["PubChem Compound"] = drug.findtext(
        "{ns}external-identifiers/{ns}external-identifier[{ns}resource='PubChem Compound']/{ns}identifier".format(
            ns=ns
        )
    )
    row["approved"] = drug.findtext("{ns}groups/{ns}group".format(ns=ns))
    row["atc_codes"] = [
        code.get("code")
        for code in drug.findall("{ns}atc-codes/{ns}atc-code".format(ns=ns))
    ]
    row["targets"] = [
        target.text
        for target in drug.findall(
            (
                "{ns}targets/{ns}target/{ns}polypeptide/{ns}external-identifiers/"
                "{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier"
            ).format(ns=ns)
        )
    ]

    rows.append(row)
# +
columns = ["targets", "atc_codes", "approved", "type", "PubChem Compound"]

drugbank_df = pd.DataFrame.from_dict(rows)[columns]

display(len(drugbank_df[drugbank_df["approved"] == "approved"]))
drugbank_df = drugbank_df[(drugbank_df["approved"] == "approved")]
# drugbank_df = drugbank_df[(drugbank_df['approved'] == 'approved') & ((drugbank_df['type'] == 'small molecule'))]
display(len(drugbank_df))
drugbank_df.drop(["approved", "type"], axis=1, inplace=True)
display(drugbank_df.head())

# +
drugbank_size = drugbank_df.shape[0]

# cleanup rows with no target
drugbank_df = drugbank_df[drugbank_df["targets"].str.len() != 0]
drugbank_notarget_size = drugbank_df.shape[0]

drugbank_df.dropna(subset=["PubChem Compound"], inplace=True)
drugbank_nochem_size = drugbank_df.shape[0]

drugbank_df = drugbank_df[drugbank_df["atc_codes"].str.len() != 0]

print(
    "cleaning drugs with no targets removed "
    + str(drugbank_size - drugbank_notarget_size)
    + " drugs from "
    + str(drugbank_size)
    + ", leaving "
    + str(drugbank_notarget_size)
)
print(
    "cleaning drugs with no PubChem ID removed "
    + str(drugbank_notarget_size - drugbank_nochem_size)
    + ", leaving "
    + str(drugbank_nochem_size)
)
print(
    "cleaning drugs with no ATC removed "
    + str(drugbank_nochem_size - drugbank_df.shape[0])
    + ", leaving "
    + str(drugbank_df.shape[0])
)

display(drugbank_df.head())
# -

print("out of:- " + str(drugbank_notarget_size))
print(
    "PubChem Compound: "
    + str(len(drugbank_df[~(drugbank_df["PubChem Compound"].isna())]))
)
print("ATC codes: " + str(len(drugbank_df[~(drugbank_df["atc_codes"].isna())])))

# Now we need to convert the protein targets to genes, using the ID mapping table from http://www.uniprot.org/downloads

# +
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
display(idmap_full.head())


# +
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


# Tidy up lists of multiple proteins:

# +
atc_level_genes = drugbank_df.copy()

atc_level_genes["genes"] = [[] for i in range(len(atc_level_genes))]

for index, row in atc_level_genes.iterrows():
    for prot in row["targets"]:
        if prot in idmap.index:
            gene = idmap.at[prot, "EntrezGene"]
            if gene not in atc_level_genes.loc[index]["genes"]:
                atc_level_genes.at[index, "genes"] = atc_level_genes.at[
                    index, "genes"
                ] + [str(gene)]

display(atc_level_genes.head())

# +
# We perform this cutoff once all DGI DB results are integrated, not here - so set to 1.
geneset_cutoff = 1
print("geneset list length before: " + str(len(atc_level_genes)))
atc_level_genes_genefilter = atc_level_genes[
    atc_level_genes["genes"].str.len() >= geneset_cutoff
]
print("geneset list length after: " + str(len(atc_level_genes_genefilter)))
print(
    "genesets >= 5 genes: "
    + str(len(atc_level_genes[atc_level_genes["genes"].str.len() >= 5]))
)

display(atc_level_genes_genefilter.head())
# -

atc_level_genes_genefilter["genes"] = [
    sorted(l) for l in atc_level_genes_genefilter["genes"]
]
atc_level_genes_genefilter["genes"] = atc_level_genes_genefilter["genes"].apply(
    " ".join
)
atc_level_genes_genefilter["atc_codes"] = atc_level_genes_genefilter["atc_codes"].apply(
    " ".join
)

display(atc_level_genes_genefilter.shape)
display(atc_level_genes_genefilter.head())

atc_level_genes_genefilter = atc_level_genes_genefilter[
    atc_level_genes_genefilter["genes"] != "nan"
]
display(atc_level_genes_genefilter.shape)

# save the atc level 5 group of drugs against their names to compare with other drug target databases
compare_save_df = atc_level_genes_genefilter[["PubChem Compound", "atc_codes", "genes"]]
####compare_save_df.drop_duplicates('PubChem Compound', inplace=True)
compare_save_df = compare_save_df[compare_save_df["atc_codes"] != ""]
display(compare_save_df.head())
display(len(compare_save_df))
compare_save_df.to_csv(
    "../data/target-dbs/drugbank_atc_targets.csv", index=False, sep="\t"
)

compare_save_df[compare_save_df["atc_codes"].str.contains("L01XE4")]

compare_save_df[compare_save_df["PubChem Compound"] == 5359476]

# Also, create a table of PMID references as evidence for each drug-target interaction:

evidence = list()
for i, drug in enumerate(root):
    assert drug.tag == ns + "drug"
    PubChem_Compound = drug.findtext(
        "{ns}external-identifiers/{ns}external-identifier[{ns}resource='PubChem Compound']/{ns}identifier".format(
            ns=ns
        )
    )
    approved = drug.findtext("{ns}groups/{ns}group".format(ns=ns))
    atc_codes = [
        code.get("code")
        for code in drug.findall("{ns}atc-codes/{ns}atc-code".format(ns=ns))
    ]

    for target in drug.findall(
        (
            "{ns}targets/{ns}target/{ns}polypeptide/{ns}external-identifiers/"
            "{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier"
        ).format(ns=ns)
    ):
        row = collections.OrderedDict()
        row["PubChemCompound"] = PubChem_Compound
        row["approved"] = approved
        row["atc_codes"] = atc_codes
        row["target"] = target.text
        if target.text in idmap.index:
            row["target"] = idmap.at[target.text, "Ensembl"]
        row["refs"] = [
            ref.text
            for ref in drug.findall(
                (
                    "{ns}targets/{ns}target/{ns}references/{ns}articles/"
                    "{ns}article/{ns}pubmed-id"
                ).format(ns=ns)
            )
        ]
        evidence.append(row)

columns = ["PubChemCompound", "atc_codes", "approved", "target", "refs"]
evidence_df = pd.DataFrame.from_dict(evidence)[columns]
display(evidence_df.head())
# evidence_df = evidence_df[(evidence_df['approved'] == 'approved') & ((evidence_df['type'] == 'small molecule'))]
evidence_df = evidence_df[(evidence_df["approved"] == "approved")]
evidence_df.drop(["approved"], axis=1, inplace=True)
evidence_df["target"] = evidence_df["target"].astype(str)
display(len(evidence_df))
display(len(evidence_df[~evidence_df["target"].str.startswith("ENS")]))
evidence_df.to_csv("../data/target-dbs/drugbank_evidence.csv", index=False, sep="\t")
