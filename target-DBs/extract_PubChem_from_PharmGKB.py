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

# # Parse PharmGKB
#
# We want to obtain a list of drug ATC codes vs. Genes
#
# https://www.pharmgkb.org/downloads/
#
# ### Drugs / Gene file

import os
import pandas as pd
from IPython.display import display

# +
pgkb_genes = pd.read_csv("../data/pharmgkb/genes/genes.tsv", sep="\t")
display(pgkb_genes.head())
display(pgkb_genes.columns)

pgkb_genes = pgkb_genes[
    [
        "PharmGKB Accession Id",
        "Ensembl Id",
        "HGNC ID",
        "Chromosome",
        "Chromosomal Start - GRCh37.p13",
        "Chromosomal Stop - GRCh37.p13",
    ]
]

display(pgkb_genes.head())
print("PharmGKB has " + str(len(pgkb_genes)) + " drugs")
# -

# Remove drugs that have no Ensembl ID
pgkb_genes = pgkb_genes[pgkb_genes["Ensembl Id"].notnull()]
display(pgkb_genes.head())
print(
    "PharmGKB has "
    + str(len(pgkb_genes))
    + " drugs, after removing drugs with no Ensembl ID"
)

# ### Drugs / ATC file

# +
pgkb_atc = pd.read_csv("../data/pharmgkb/drugs/drugs.tsv", sep="\t")

pgkb_atc = pgkb_atc[["PharmGKB Accession Id", "External Vocabulary"]]

display(pgkb_atc.head())

print("PharmGKB has " + str(len(pgkb_atc)) + " drugs")
# -

# Remove drugs that have no ATC code
pgkb_atc = pgkb_atc[pgkb_atc["External Vocabulary"].notnull()]
display(pgkb_atc.head())
print(
    "PharmGKB has "
    + str(len(pgkb_atc))
    + " drugs, after removing drugs with no ATC code"
)
