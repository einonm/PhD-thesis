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

# ## Parse COVID-19 suscepibility GWAS
#
# Here we want to convert the SNP ID column in chr/pos/allele format to an rs number.
#
# Taken from https://www.covid19hg.org/results/

import os
import pandas as pd
from IPython.display import display

signif_pairs_path = os.path.join(
    "/scratch/c.mpnme/data/gnt-data/gwas/COVID19_HGI_ANA5_20200513.txt"
)
signif_pairs = pd.read_csv(signif_pairs_path, sep="\t")
display(signif_pairs.head())

len(signif_pairs)

# Extract the MarkerName to chr / bp columns

# Convert chr and pos to an rs number using Kaviar

kaviar = pd.read_csv(
    "../data/wrangling/kaviar_parsed.tsv",
    sep="\t",
    index_col=0,
    names=["#CHR", "POS", "RS"],
)

# Convert chr column to string, to allow merge below to complete

kaviar.head()

joined_df = pd.merge(kaviar, signif_pairs, on=["#CHR", "POS"])
display(joined_df.head())

joined_df = joined_df[["#CHR", "POS", "RS", "all_inv_var_meta_p", "all_inv_var_het_p"]]
joined_df.rename(columns={"#CHR": "CHR"}, inplace=True)
joined_df.head()

joined_df.to_csv(
    os.path.join(
        "/scratch/c.mpnme/data/gnt-data/gwas/", "COVID19_HGI_ANA5_20200513_rs.txt"
    ),
    sep="\t",
    index=False,
    header=True,
)

joined_df[joined_df.duplicated(subset="CHR") == True]
