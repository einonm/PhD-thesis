# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Parse PsychENCODE data
#
# Inputs:
#  1. PsychENCODE eQTL set with FDR<0.05 and a filter requiring genes to have an expression > 0.1 FPKM in at least 10 samples (DER-08a_hg19_eQTL.significant.txt from http://resource.psychencode.org/)
#  2. SNP information for all QTLs considered, including rsIDs (if available), location, and reference and alternate alleles (http://resource.psychencode.org/Datasets/Derived/QTLs/SNP_Information_Table_with_Alleles.txt)
#
# Outputs:
#  1. SNP to gene mapping file
#
# Methods:
# Data wrangling with pandas

# +
import os
import pandas as pd
import sys
from IPython.display import display

import re

# -

# Don't do this lengthy step, the file psychencode_out.tsv probably exists - goto #LOAD PSYCHENCODE_OUT
psychencode_df = pd.read_csv(
    "../data/wrangling/PsychENCODE/DER-08a_hg19_eQTL.significant.txt",
    sep="\t",
    usecols=["gene_id", "SNP_id", "gene_chr", "gene_start", "gene_end"],
)

psychencode_df

alleles_df = pd.read_csv(
    "../data/wrangling/PsychENCODE/SNP_Information_Table_with_Alleles.txt",
    sep="\t",
    usecols=["PEC_id", "Rsid"],
)

alleles_df

merged_df = pd.merge(
    psychencode_df, alleles_df, how="left", left_on="SNP_id", right_on="PEC_id"
)

# Remove those entries not having an rs number
merged_df = merged_df[merged_df["Rsid"].str.startswith("rs")]
merged_df.drop(["SNP_id", "PEC_id"], axis=1, inplace=True)

# Remove trailing '.1' numbers from gene IDs
merged_df["gene_id"] = [x[0:15] for x in merged_df["gene_id"]]

# +
#### eqtlgen_df.set_index('SNP', inplace=True)
display(merged_df.head())

uniq_genes = len(pd.value_counts(merged_df["gene_id"]))
print("There are " + str(uniq_genes) + " unique genes.")
# -

merged_df.set_index("Rsid", inplace=True)

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

for index, row in merged_df.iterrows():
    if row["gene_id"] in annot.index:
        annot.loc[row["gene_id"]]["snplist"] = (
            annot.at[row["gene_id"], "snplist"] + " " + index
        )
    else:
        annot.at[row["gene_id"], "snplist"] = index

display(annot.head())
# -

annot

# Save intermediate file after previous lengthy processing
annot.to_csv("../data/wrangling/psychencode_out.tsv", sep="\t", header=False)

# LOAD PSYCHENCODE_OUT
# The psychencode_out file takes some time to generate, so re-load here
annot2 = pd.read_csv(
    "../data/wrangling/psychencode_out.tsv",
    sep="\t",
    header=None,
    names=["Gene", "snplist"],
    index_col=0,
)

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

annot_mrg.to_csv("../data/wrangling/psychencode_final.tsv", header=None, sep="\t")

# # QC

annot_mrg

annot_mrg[annot_mrg.index.duplicated() == True]

annot_mrg.isnull().any()

# Find average & std dev of snps in annotation data set
annot_mrg['counts'] = annot_mrg['snpstring'].str.split().str.len()
annot_mrg.describe()


