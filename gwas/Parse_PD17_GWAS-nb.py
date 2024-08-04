# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Parse HD GWAS
#
# Here we want to convert the SNP ID column in chr/pos/allele format to an rs number.

import os
import pandas as pd
from IPython.display import display

signif_pairs_path = os.path.join(
    "/home/mpnme/data/PD-23am/METAANALYSIS_PdGeneAnd23AndMe1.tbl"
)
signif_pairs = pd.read_csv(signif_pairs_path, sep="\t")
display(signif_pairs.head())

# Extract the MarkerName to chr / bp columns

tester_df = signif_pairs.head()
signif_pairs["chr"] = signif_pairs["MarkerName"].str.split(":").str[0].str[3:]
signif_pairs["bp"] = signif_pairs["MarkerName"].str.split(":").str[1]
signif_pairs.head()

signif_pairs["chr"] = pd.to_numeric(signif_pairs["chr"])
signif_pairs["bp"] = pd.to_numeric(signif_pairs["bp"])
display(type(signif_pairs["chr"][0]))
display(type(signif_pairs["bp"][0]))

signif_pairs.count()

# Convert chr and pos to an rs number using Kaviar

gene_pos = pd.read_csv(
    "../data/wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    index_col=0,
    names=["Ensembl", "Entrez", "chrom", "chromStart", "chromEnd", "HGNC"],
)

# Convert chr column to string, to allow merge below to complete

# kaviar['chr'] = kaviar['chr'].astype(str)
kaviar.columns = ["chr", "bp", "snp_id"]
display(kaviar.head())

display(type(kaviar["chr"][0]))
display(type(kaviar["bp"][0]))

joined_df = pd.merge(kaviar, signif_pairs, on=["chr", "bp"])
display(joined_df.head())

joined_df.isnull().any(axis=1)

joined_df.to_csv(
    os.path.join("/home/mpnme/data/PD-23am/", "METAANALYSIS_PdGeneAnd23AndMe1-rs.tbl2"),
    sep="\t",
    index=False,
)

joined_df[joined_df.duplicated(subset="chr") == True]
