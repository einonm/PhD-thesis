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

# # Exclude MHC block from GWAS
#
# (Major Histocompatibility Complex block)
#
# For each association file (product of the GWAS), remove any SNPs within the MHC block, defined to be on Chr 6 between position 26000000 and 33000000. save it with the .no_mhc extension.

import os
import pandas as pd
from IPython.display import display

# +
prefix = "/home/mpnme/data/GWAS/"
gwas_list = [
    "Cluster1.assoc",
    "Cluster2.assoc",
    "Cluster3.assoc",
    "Cluster4.assoc",
    "daner_PGC_SCZ52_0513a.txt",
    "Log.McGill.assoc",
    "PD_2014_META.txt",
]

for gwas_file in gwas_list:
    # Import some columns as text, as conversion leads to rounding errors
    gwas = pd.read_csv(prefix + gwas_file, delim_whitespace=True, dtype=object)
    chr_col = ""
    bp_col = ""
    for col in list(gwas):
        if col.casefold() == "chr":
            chr_col = col
        if col.casefold() == "bp":
            bp_col = col

    discard_rows = gwas.loc[
        (gwas[chr_col] == "6")
        & (pd.to_numeric(gwas[bp_col]) >= 26000000)
        & (pd.to_numeric(gwas[bp_col]) <= 33000000)
    ]
    display(discard_rows.head())

    filtered_gwas = (
        pd.merge(gwas, discard_rows, how="outer", indicator=True)
        .query('_merge == "left_only"')
        .drop(["_merge"], axis=1)
    )
    filtered_gwas.to_csv(prefix + gwas_file + ".no_mhc", index=False, sep=" ")
    print(
        "orig: "
        + str(len(gwas))
        + "  discarded: "
        + str(len(discard_rows))
        + "  filtered: "
        + str(len(filtered_gwas))
    )
# -
