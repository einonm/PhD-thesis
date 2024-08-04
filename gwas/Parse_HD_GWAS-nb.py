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

# ## Parse HD GWAS
#
# Here we want to convert the SNP ID column in chr/pos/allele format to an rs number.

import os
import pandas as pd
from IPython.display import display

signif_pairs_path = os.path.join("../data/gwas/", "GWA12345_summary_rs_snps.txt")
signif_pairs = pd.read_csv(signif_pairs_path, sep="\t")
display(signif_pairs.head())

# Convert chr and pos to an rs number using Kaviar
data_path = "../data/"
kaviar = pd.read_csv(data_path + "../data/wrangling/kaviar_parsed.tsv", sep="\t")

# Convert chr column to string, to allow merge below to complete

# kaviar['chr'] = kaviar['chr'].astype(str)
kaviar.columns = ["chr", "bp", "snp_id"]
display(kaviar.head())

joined_df = pd.merge(kaviar, signif_pairs, on=["chr", "bp"])
display(joined_df.head())

joined_df.to_csv(
    os.path.join("../data/gwas/", "GWA12345_summary_rs_snps.txt"), sep="\t", index=False
)
