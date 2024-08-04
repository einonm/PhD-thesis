# ---
# jupyter:
#   jupytext:
#     notebook_metadata_filter: language_info
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Python 3
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
#     version: 3.8.5
# ---

# # Convert proximity annotation file to Ensembl
#
# In order to perform a gene analysis using the default MAGMA proximity annotation method, along with ATC gene sets, we'll need to convert the generated proximity annotation file from entrez to ensembl.

import os
import sys
import pandas as pd
from IPython.display import display
import warnings

warnings.filterwarnings("always")

# Use the script `create_magma_prox_annot.sh` to generate a gene annotation from `magma_auxiliary_files/g1000_eur/NCBI37.3.gene.loc`
#
# As the generated gene annot file has a variable number of columns (due to the snp list), outside of the notebook, replace the first two occurrences of '\t' with ';' e.g. by twice running:
#
#     :%s/^I/;/
#
# in vim. Save as 'file-comsep'.

# +
annot_loc = pd.read_csv(
    "/c8000xd3/big-mpnme/data/magma_g1000eur_NCBI37.3.prox.genes.annot-comsep",
    comment="#",
    sep=";",
    header=None,
    names=["Entrez ID", "chr:start:end", "snplist"],
)

display(annot_loc.head())
print("number of genes: " + str(len(annot_loc)))

# +
# remove \t's from snplist
annot_loc["snplist"] = annot_loc["snplist"].str.replace("\t", " ")

display(annot_loc.head())

# +
data_path = "/home/mpnme/data/phd/wrangling/"
idmap = pd.read_csv(data_path + "Ensembl_pos.tsv", sep="\t")

#  remove trailing '.0' from EntrezID
idmap["Entrez ID"] = idmap["Entrez ID"].astype("int").astype("str")

display(idmap.head())
print("number of genes: " + str(len(idmap)))
# -

# Look at this! multiple Entrez IDs for one HGNC/Ensembl code - wtf??
# Let's not use the Entrez ID from this table as it can't be trusted / is not understood - let's remove it
# See https://www.biostars.org/p/16505/
display(idmap[idmap["Ensembl"] == "ENSG00000235249"])

# idmap.drop('Entrez ID', axis=1, inplace=True)
display(idmap.head())
idmap.drop_duplicates(inplace=True)
print("number of genes: " + str(len(idmap)))

idmap["Entrez ID"] = pd.to_numeric(idmap["Entrez ID"])
# idmap.set_index('Ensembl', inplace=True)
display(idmap.head())
display([idmap.index.duplicated() == True])

annot_ensembl = pd.merge(
    annot_loc, idmap, left_on="Entrez ID", right_on="Entrez ID", how="inner"
)  # [['Ensembl','Entrez ID', 'chr', 'start', 'end', 'gene']]
annot_ensembl.set_index("Ensembl", inplace=True)
annot_ensembl.sort_index(inplace=True)
display(annot_ensembl.head())
print("number of genes: " + str(len(annot_ensembl)))

# +
# check for NaN in index, and remove
annot_ensembl = annot_ensembl[annot_ensembl.index.notnull()]

# check for duplicated Ensembl IDs, following conversion from Entrez
annot_ensembl = annot_ensembl[~annot_ensembl.index.duplicated()]

display(annot_ensembl.shape)
# -

# remove columns not required
annot_ensembl_final = annot_ensembl.drop(
    ["Entrez ID", "HGNC symbol", "chrom", "chromStart", "chromEnd", "version"], axis=1
)

# +
display(annot_ensembl_final.shape)
annot_ensembl_final.head()

## NOTE - why is this at 19,863 greater than the original annot_loc list at 18,155?
# -

annot_ensembl_final.to_csv(
    data_path + "../../" + "magma_g1000eur_NCBI37.3.prox.genes.annot.ensemble",
    sep="\t",
    header=False,
)
