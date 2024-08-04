# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Convert magma gene loc file to Ensembl
#
# In order to perform a gene analysis using the default MAGMA proximity annotation method, along with ATC gene sets, we'll need to convert the magma gene location file from entrez to ensembl.

import os
import sys
import pandas as pd
from IPython.display import display
import warnings

warnings.filterwarnings("always")

# +
annot_loc = pd.read_csv(
    "../data/wrangling/NCBI37.3.gene.loc",
    comment='#',
    delim_whitespace=True,
    header=None,
    names=['Entrez ID', 'chr', 'start', 'end', 'sense', 'gene'],
)

annot_loc["Entrez ID"] = annot_loc["Entrez ID"].astype("int").astype("str")

display(annot_loc.head())
print("number of genes: " + str(annot_loc.shape[0]))

# +
data_path = '../data//wrangling/'
idmap = pd.read_csv(
    data_path + 'Ensembl_pos.tsv', sep='\t'
)  #[['Entrez ID', 'Ensembl']] #index_col='Entrez ID',

#  remove trailing '.0' from EntrezID
idmap["Entrez ID"] = idmap["Entrez ID"].astype("int").astype("str")

display(idmap.head())
print("number of genes: " + str(len(idmap)))

# +
# Look at this! multiple Entrez IDs for one HGNC/Ensembl code - wtf??
# Let's not use the Entrez ID from this table as it can't be trusted / is not understood - let's remove it
# See https://www.biostars.org/p/16505/
display(idmap[idmap["Ensembl"] == "ENSG00000235249"])

idmap.drop("Entrez ID", axis=1, inplace=True)
display(idmap.head())
idmap.drop_duplicates(inplace=True)
print("number of genes: " + str(len(idmap)))
# -

annot_ensembl = pd.merge(
    annot_loc, idmap, left_on="gene", right_on="HGNC symbol", how="left"
)[["Ensembl", "Entrez ID", "chr", "start", "end", "gene"]]
annot_ensembl.set_index("Ensembl", inplace=True)
display(annot_ensembl.head())
print("number of genes: " + str(len(annot_ensembl)))

# +
# annot_ensembl.to_csv(data_path + 'NCBI37.3.ensembl.gene.loc', sep = "\t", header=False)
# -
