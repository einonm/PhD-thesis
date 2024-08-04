# ---
# jupyter:
#   jupytext:
#     notebook_metadata_filter: language_info
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
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
#     version: 3.7.11
# ---

# ## Gather all gene / snp files into one list
#
# This collates gene-snp list files from two directories into one combined file & doubles check for duplicates in snp list once collated.

# imports
import os
import numpy as np
import pandas as pd
from IPython.display import display
import re

# +
# get list of qtls we want to combine.
# Commented data sets are currently in use, uncomment one set to run.

combined_prefix = "all-eqtl-brain"
file_prefix1 = "../data/qtl/all-eqtl-brain/"
file_prefix2 = "" # No second dir 

#combined_prefix = "all-brain"
#file_prefix1 = "../data/qtl/all-eqtl-brain/"
#file_prefix2 = "../data/qtl/mqtl"

#combined_prefix = "all-eqtl-blood"
#file_prefix1 = "../data/qtl/all-eqtl-blood/"
#file_prefix2 = "" # No second dir 

#combined_prefix = "all-eqtl-tissues"
#file_prefix1 = "../data/qtl/all-eqtl-tissues/"
#file_prefix2 = "" # No second dir 

#combined_prefix = "all-tissues"
#file_prefix1 = "../data/qtl/all-eqtl-tissues/"
#file_prefix2 = "../data/qtl/mqtl"

# Previously used data sets

#combined_prefix = "combined"
#file_prefix1 = "../data/qtl/qtls-all/"
#file_prefix2 = "" # No second dir

#combined_prefix = "new_combined"
#file_prefix1 = "../data/qtl/new-qtl/"
#file_prefix2 = "" # No second dir

#combined_prefix = "everything_combined"
#file_prefix1 = "../data/qtl/new-qtl/"
#file_prefix2 = "../data/qtl/qtl-all/"

# +
qtl_files1 = []
qtl_files2 = []
for root, dirs, files in os.walk(file_prefix1):
    qtl_files1 = files

for root, dirs, files in os.walk(file_prefix2):
    qtl_files2 = files
# -

combined_df = pd.read_csv(
    os.path.join(file_prefix1, qtl_files1.pop()),
    sep="\t",
    header=None,
    index_col=0,
    names=["Ensembl", "position", "snplist"],
)

# Append qtls for each
def combine_qtls(combined_df, new_df):
    combined_df = pd.merge(combined_df, new_df, how="outer", on=["Ensembl"])
    combined_df.fillna("", inplace=True)
    # Sort out position column, create one column with a non-empty string
    combined_df["position"] = np.where(
        combined_df["position_x"] == "",
        combined_df["position_y"],
        combined_df["position_x"],
    )
    combined_df.drop(["position_x", "position_y"], axis=1, inplace=True)
    # Sort out snplist, create one long string of all snps
    combined_df["snplist"] = combined_df["snplist_x"] + " " + combined_df["snplist_y"]
    combined_df.drop(["snplist_x", "snplist_y"], axis=1, inplace=True)
    return combined_df


# +
for new_file in qtl_files1:
    new_df = pd.read_csv(
        os.path.join(file_prefix1, new_file),
        sep="\t",
        header=None,
        index_col=0,
        names=["Ensembl", "position", "snplist"],
    )
    combined_df = combine_qtls(combined_df, new_df)

for new_file in qtl_files2:
    new_df = pd.read_csv(
        os.path.join(file_prefix2, new_file),
        sep="\t",
        header=None,
        index_col=0,
        names=["Ensembl", "position", "snplist"],
    )
    combined_df = combine_qtls(combined_df, new_df)
# -

combined_df.shape

# remove duplicate SNPs from SNP lists
combined_df["snplist"] = combined_df.apply(
    lambda x: " ".join(sorted(set(x["snplist"].split()))), axis=1
)
combined_df.shape

# Remove any X/Y chromosomes
combined_df = combined_df[~combined_df["position"].str.startswith("X")]
combined_df = combined_df[~combined_df["position"].str.startswith("Y")]
combined_df.shape

# Remove duplicate SNP lists from genes
new_combined_df = combined_df.drop_duplicates()
new_combined_df.shape

new_combined_df.to_csv("../data/magma_analysis/" + combined_prefix + ".gene.annot-QTLs",
                   sep = "\t", header=False)

# ## Results QC

new_combined_df

# There should be no duplicates
new_combined_df[new_combined_df.index.duplicated() == True]

# No values should be null, so this should be false
new_combined_df.isnull().any()

# Find average & std dev of snps in annotation data set
counts_df = pd.DataFrame()
counts_df['counts'] = new_combined_df['snplist'].str.split().str.len()
counts_df.describe()
