# ---
# jupyter:
#   anaconda-cloud: {}
#   jupytext:
#     notebook_metadata_filter: anaconda-cloud,language_info
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.7.6
# ---

# # Extract ATC code gene sets from DrugBank
#
# This notebook takes the full drugbank database and extracts each drug that has a protein target, the protein target and the drug ATC code.

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

# Load in the drugbank XML (5.1.0 here, obtained from https://www.drugbank.ca/releases/latest):

# +
data_path = "../data/drugbank/"

tree = ET.parse(data_path + "full-drugbank-database-5.1.0.xml")
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
    row["atc_codes"] = [
        code.get("code")
        for code in drug.findall("{ns}atc-codes/{ns}atc-code".format(ns=ns))
    ]
    row["target"] = [
        target.text
        for target in drug.findall(
            "{ns}targets/{ns}target/{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier".format(
                ns=ns
            )
        )
    ]

    rows.append(row)
# -

# Tidy up lists of multiple proteins:

# +
columns = ["drugbank_id", "name", "atc_codes", "target"]
drugbank_df = pd.DataFrame.from_dict(rows)[columns]
drugbank_df.set_index("drugbank_id", inplace=True)

drugbank_df.head()
# -

# ### Only consider drugs whose target gene sets have more than 5 members
#
# (also see cutoff following conversion to ensembl below)

# * We exclude drugs that have one (or perhaps two) targets, as they are uninteresting for enrichment analysis - you can't find an enrichment of just one gene, it has to be a group.
#
# * The set of gene targets for each ATC group (at levels 1,2,3,4 and 5) is compared to the set of gene targets from the GWAS/QTL analysis. Strong associations between them are what we are after (as found by magma).

# +
drugbank_size = drugbank_df.shape[0]

# Add other manually-curated ATC codes of drugbank
# drugbank_df_noatc = drugbank_df_genefilter[drugbank_df_genefilter.atc_codes.str.len() == 0]
# drugbank_df_noatc.to_csv("data/no_atc_code.csv", sep='\t', header=None)

# cleanup rows with no ATC code
drugbank_df_atc = drugbank_df[drugbank_df.atc_codes.str.len() != 0]

# cleanup rows with no target
drugbank_df_atc = drugbank_df_atc[drugbank_df_atc["target"].str.len() != 0]

print(
    "cleaning drugs with no ATC code or targets removed "
    + str(drugbank_size - drugbank_df_atc.shape[0])
    + " drugs from "
    + str(drugbank_size)
    + ", leaving "
    + str(drugbank_df_atc.shape[0])
)

display(drugbank_df_atc.head())
# -

# Set the ATC level of the gene set to be generated (1-5) in atc_level below.

# +
atc_level = 5

atc_level_chars = [
    0,
    1,
    3,
    4,
    5,
    7,
]  # number of characters in the ATC code for each level
atc_level_total = [0, 14, 93, 267, 885, 4823]  # total number of entries at each level
# -

# Now, let's get a list of all proteins at the given ATC level.

# +
atc_level_str = "atc_level_" + str(atc_level)
drugbank_df_atc[atc_level_str] = ""
atc_level_full = pd.DataFrame(columns=(atc_level_str, "drugbank_ids", "targets"))
atc_level_full.set_index(atc_level_str, inplace=True)

for index, row in drugbank_df_atc.iterrows():
    for item in row["atc_codes"]:
        atc_nth = item[: atc_level_chars[atc_level]]
        if atc_nth in atc_level_full.index:
            atc_level_full.at[atc_nth, "drugbank_ids"] = (
                atc_level_full.at[atc_nth, "drugbank_ids"] + "    " + index
            )
            atc_level_full.at[atc_nth, "targets"] = (
                atc_level_full.at[atc_nth, "targets"] + row["target"]
            )
        else:
            atc_level_full.at[atc_nth, "drugbank_ids"] = index
            atc_level_full.at[atc_nth, "targets"] = row["target"]


print(
    "There are "
    + str(len(atc_level_full))
    + " level "
    + str(atc_level)
    + " atc drug/target groups (out of "
    + str(atc_level_total[atc_level])
    + " available groups)"
)
display(atc_level_full.head(20))
# -

# Now we need to convert the protein targets to genes, using the ID mapping table from http://www.uniprot.org/downloads

# +
idmap_path = data_path + "idmapping-ENS.tab"
idmap_full = pd.read_csv(
    idmap_path,
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
        "EMBL",
        "EMBL-CDS",
        "Ensembl",
        "Ensembl_TRS",
        "Ensembl_PRO",
        "PubMed",
    ],
)

idmap = idmap_full[["Ensembl"]]

# Remove multiple gene IDs if they exist
# TODO - why are there > 1 genes listed for some proteins?
for index, row in idmap.iterrows():
    if str(row["Ensembl"]).find(";") != -1:
        idmap.loc[index]["Ensembl"] = row["Ensembl"][: row["Ensembl"].index(";")]

# Remove MHC region genes?
mhc_genes = pd.read_csv("../qtl/data/mhc_genes.txt", sep="\t")["Ensembl"]
mhc_genes = mhc_genes.unique()
idmap = idmap.loc[~idmap["Ensembl"].isin(mhc_genes)]

display(idmap.head())


# +
atc_level_genes = atc_level_full.copy()

atc_level_genes["genes"] = [[] for i in range(len(atc_level_genes))]

for index, row in atc_level_genes.iterrows():
    for prot in row["targets"]:
        if prot in idmap.index:
            gene = idmap.loc[prot]["Ensembl"]
            if gene not in atc_level_genes.loc[index]["genes"]:
                atc_level_genes.loc[index]["genes"] = atc_level_genes.at[
                    index, "genes"
                ] + [gene]

del atc_level_genes["drugbank_ids"]

display(atc_level_genes.head())

# +
geneset_cutoff = 5
print("geneset list length before: " + str(len(atc_level_genes)))
atc_level_genes_genefilter = atc_level_genes[
    atc_level_genes["genes"].str.len() >= geneset_cutoff
]
print("geneset list length after: " + str(len(atc_level_genes_genefilter)))

display(atc_level_genes_genefilter.head())

# +
atc_level_genes_genefilter["targets"] = atc_level_genes_genefilter["genes"].apply(
    " ".join
)
del atc_level_genes_genefilter["genes"]

display(atc_level_genes_genefilter.head())
# -

# ## Add some canary drugs
#
# Add the SNCA and MAPT genes, which are known to be associated with PD.
#
# We do this by creating a 'wonderdrug' with ATC code V07AA01 (Basically, a type of plaster!) to avoid future confusion

# +
# Disabled for now, as the outliers it introduces can skew correlation results

# new_genes = 'ENSG00000145335 ENSG00000186868' # SNCA and MAPT
# new_atc = 'V07AA01'

# if atc_level == 1:
#    atc_level_genes_genefilter.loc['V']['targets'] += new_genes
# else:
#    atc_level_genes_genefilter.loc[new_atc[:atc_level_chars[atc_level]]] = new_genes

# print(new_atc[:atc_level_chars[atc_level]] + " value is " +
#      atc_level_genes.loc[new_atc[:atc_level_chars[atc_level]]]['genes'])
# -

# ## Don't include ATC entries with matching target lists
#
# This removes entries with duplicate gene lists.

# +
dedup_df = atc_level_genes_genefilter.copy()

dedup_df.drop_duplicates(inplace=True)
dedup_df.reset_index(inplace=True)

print("geneset list length before dedup: " + str(len(atc_level_genes_genefilter)))
print("geneset list length after dedup: " + str(len(dedup_df)))

dedup_df = dedup_df.groupby(by="targets")[atc_level_str, "targets"].head()

# TODO - using the amalgamated name breaks the results presentation.
# Use a seperate file for amalgamated ATC codes?
dedup_df.set_index(atc_level_str, inplace=True)

dedup_df.head()
# -

# Sort for ease of comparison and write these results to a file:

# +
dedup_df.sort_index(inplace=True)

dedup_df.to_csv(
    data_path + "dbank_gene_set-atc" + str(atc_level) + "no_mhc.txt",
    sep="\t",
    header=None,
)
# -
