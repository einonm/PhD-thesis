# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Parse eQTLgen data
#
# Inputs:
#  1. eQTLgen cis-eQTLs with less than FDR0.05 (from https://www.eqtlgen.org/cis-eqtls.html)
#
# Outputs:
#  1. SNP to gene mapping file
#
# Methods:
# Data wrangling with pandas

# +
import os
import pandas as pd
import numpy as np
import sys
from IPython.display import display

import re

# -

# Don't do this lengthy step, the file eqtlgen_out.tsv probably exists - goto #LOAD EQTLGEN_OUT
eqtlgen_df = pd.read_csv(
    "../data/wrangling/eQTLgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
    sep="\t",
    usecols=["SNP", "Gene", "GeneChr", "GenePos"],
)

# Format gene position - not no end position is available! Not including as it's not used
eqtlgen_df["pos"] = (
    eqtlgen_df["GeneChr"].astype(str) + ":" + eqtlgen_df["GenePos"].astype(str)
)
eqtlgen_df.drop(["GeneChr", "GenePos"], axis=1, inplace=True)

# +
eqtlgen_df.set_index("SNP", inplace=True)
display(eqtlgen_df)

uniq_genes = len(pd.value_counts(eqtlgen_df["Gene"]))
print("There are " + str(uniq_genes) + " unique genes.")
# -

# ## Create MAGMA annotation file
#
# We'd now like to create an annotation file (mapping SNPs to genes) to use with MAGMA. By default, MAGMA uses genome proximity of SNPs to genes to do this. We are using sQTL data, so we have to manually create the annotation file.
#
# An annotation file is a set of lines, each of the format:
#
# <tt>gene_name   chr:start_pos:stop_pos   snp1 snp2 snp3 snp4...etc</tt>
#
# Where we've chosen Ensembl IDs to identify genes.
#
# So let's list, for each unique Ensembl ID, the list of snps associated with it:

# +
annot = pd.DataFrame(columns=("Gene", "pos", "snplist"))
annot.set_index("Gene", inplace=True)

for index, row in eqtlgen_df.iterrows():
    if row["Gene"] in annot.index:
        annot.loc[row["Gene"]]["snplist"] = (
            annot.at[row["Gene"], "snplist"] + " " + index
        )
    else:
        annot.at[row["Gene"], "snplist"] = index
        annot.at[row["Gene"], "pos"] = row["pos"]

display(annot.head())
# -

# Save intermediate file after previous lengthy processing
annot.to_csv("../data/wrangling/eqtlgen_out.tsv", sep="\t", header=False)

# LOAD EQTLGEN_OUT
# The eqtlgen_out file takes some time to generate, so re-load here
annot2 = pd.read_csv(
    "../data/wrangling/eqtlgen_out.tsv",
    sep="\t",
    header=None,
    names=["Gene", "snplist"],
    index_col=0,
)
annot2

# +
gene_pos = pd.read_csv("../data/wrangling/Ensembl_pos.tsv", sep="\t", index_col=0)
gene_pos = gene_pos[["chrom", "chromStart", "chromEnd"]]

display(gene_pos.head())
# -

annot_mrg = pd.merge(gene_pos, annot2, how="right", left_index=True, right_index=True)
annot_mrg[annot_mrg.index.duplicated(keep=False)]

# +
# Ensembl_pos.tsv may have duplicated lines for alternate Entrez IDs, but with matching positions - remove dupes.
annot_mrg = annot_mrg.drop_duplicates()

# Not all genes have positions available, remove these too
annot_mrg = annot_mrg[pd.isnull(annot_mrg["chrom"]) == False]
display(annot_mrg.head())
display(len(annot_mrg))
# -

annot_mrg["pos"] = (
    annot_mrg["chrom"].astype(str)
    + ":"
    + annot_mrg["chromStart"].astype(int).astype(str)
    + ":"
    + annot_mrg["chromEnd"].astype(int).astype(str)
)
annot_mrg.drop(["chrom", "chromStart", "chromEnd", "Gene"], axis=1, inplace=True)
annot_mrg = annot_mrg[["pos", "snplist"]]

annot_mrg

# Convert snplist string into a list and sort it
annot_mrg["snplist"] = annot_mrg["snplist"].str.split()
for snplist in annot_mrg["snplist"]:
    snplist.sort()

# Remove any snps without an rs number
for index, row in annot_mrg.iterrows():
    # Make a copy, as we can't remove items from a list being iterated over
    snplist = row["snplist"].copy()
    for rsnum in row["snplist"]:
        if ":" in rsnum:
            snplist.remove(rsnum)
        else:
            # list is alphabetical, so all 'chr:xxx' items are gone - try the next list
            break
    annot_mrg.at[index, "snplist"] = snplist

# Check if any snplist contains duplicate SNPs
annot_mrg["duplicate"] = [
    0 if len(set(x)) == len(x) else 1 for x in annot_mrg["snplist"]
]
annot_mrg["duplicate"].any()

annot_mrg["snpstring"] = [" ".join(x) for x in annot_mrg["snplist"]]

annot_mrg.drop(["snplist", "duplicate"], axis=1, inplace=True)

annot_mrg

annot_mrg.to_csv("../data/wrangling/eqtlgen_final.tsv", header=None, sep="\t")
