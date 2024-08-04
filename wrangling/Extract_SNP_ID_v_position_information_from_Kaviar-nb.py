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

import os
import pandas as pd

# +
# Read in table using an iterator, 10000 lines at a time.
iter_csv = pd.read_csv(
    "/home/mpnme/data/phd/wrangling/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19-trim.vcf",
    sep="\t",
    comment="#",
    usecols=[0, 1, 2, 7],
    header=None,
    names=["chr", "pos", "snp_id", "info"],
    iterator=True,
    chunksize=10000,
)

# Remove rows with no SNP ID
kaviar = pd.concat(chunk[chunk["snp_id"] != "."] for chunk in iter_csv)

kaviar.head()
# -

kaviar.shape

# remove x and y chroms, to prevent the 'chr' type being a string (and making table merging difficult)
kaviar_snps = kaviar[kaviar["chr"] != "X"]

print(kaviar_snps.shape)

kaviar_snps = kaviar_snps[kaviar_snps["chr"] != "Y"]
print(kaviar_snps.shape)

kaviar_snps = kaviar_snps[kaviar_snps["chr"] != "M"]
print(kaviar_snps.shape)

kaviar_snps.to_csv(
    "~/data/gnt-data/wrangling/kaviar_parsed.tsv",
    sep="\t",
    columns=["chr", "pos", "snp_id"],
    index=False,
)
kaviar_snps.head()
