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

# # MAGMA/Biomart gene position table comparison
#
# Using hg19, which table - MAGMA NCBI37 or Biomart download, is better at having more relevant genes?
#
# Testing using various QTL databases - PsychENCODE/eQTLgen

# +
import os
import pandas as pd
import numpy as np
import sys
from IPython.display import display

import re

# -

# Read eQTLgen dataset
eqtlgen_df = pd.read_csv(
    "../data/wrangling/eQTLgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
    sep="\t",
    usecols=["SNP", "Gene", "GeneChr", "GenePos"],
)

unique_eg_genes = pd.value_counts(eqtlgen_df["Gene"])

# Read MAGMA gene-pos and Biomart gene-pos tables
magma_gene_pos = pd.read_csv(
    "../data/wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    index_col=0,
    names=["Ensembl", "Entrez", "chrom", "chromStart", "chromEnd", "HGNC"],
)
biomart_gene_pos = pd.read_csv(
    "../data/wrangling/Ensembl_pos.tsv", sep="\t", index_col=0
)

# +
# How many eQTLgen genes are in the MAGMA table?
annot_mrg = pd.merge(
    magma_gene_pos, unique_eg_genes, how="right", left_index=True, right_index=True
)

# Matches are False, misses are True
pd.isnull(annot_mrg["chrom"]).value_counts()

# +
# How many eQTLgen genes are in the Biomart table?
annot_mrg = pd.merge(
    biomart_gene_pos, unique_eg_genes, how="right", left_index=True, right_index=True
)

# Matches are False, misses are True
pd.isnull(annot_mrg["chrom"]).value_counts()
# -

# Read PsychENCODE dataset
psychencode_df = pd.read_csv(
    "../data/wrangling/PsychENCODE/DER-08a_hg19_eQTL.significant.txt",
    sep="\t",
    usecols=["gene_id", "SNP_id", "gene_chr", "gene_start", "gene_end"],
)

# +
unique_pe_genes = pd.value_counts(psychencode_df["gene_id"])
annot_mrg = pd.merge(
    magma_gene_pos, unique_pe_genes, how="right", left_index=True, right_index=True
)

pd.isnull(annot_mrg["chrom"]).value_counts()

# +
biomart_gene_pos = pd.read_csv(
    "../data/wrangling/Ensembl_pos.tsv", sep="\t", index_col=0
)
annot_mrg = pd.merge(
    biomart_gene_pos, unique_pe_genes, how="right", left_index=True, right_index=True
)

pd.isnull(annot_mrg["chrom"]).value_counts()
# -
# From these two tests, Biomart outperforms the MAGMA table in matching genes in QTL datasets
