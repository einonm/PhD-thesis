# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Generate An Ensembl protein vs entrez table
#
# Using tables obtained from online databases, create a table of Ensemble / Gene name / Entrez ID  for use in data conversions

import os
import pandas as pd
from IPython.display import display
import re

# #!pip install biomart
from biomart import BiomartServer

# +
server = BiomartServer("http://grch37.ensembl.org/biomart")
data_path = "/home/mpnme/data/phd/wrangling/"
response_file = data_path + "biomart_response.tsv"

# server.show_databases() #ENSEMBL_MART_ENSEMBL
# server.show_datasets() #hsapiens_gene_ensembl

# use the 'hsapiens_gene_ensembl' dataset
hs_gene_ensembl = server.datasets["hsapiens_gene_ensembl"]

# hs_gene_ensembl.show_filters() # none
# hs_gene_ensembl.show_attributes()  #as below

response = hs_gene_ensembl.search(
    {
        "attributes": [
            "ensembl_gene_id",
            "ensembl_peptide_id",
            "entrezgene_id",
            "hgnc_symbol",
        ]
    }
)

biomart_response = open(response_file, "w")

for line in response.iter_lines():
    biomart_response.write("\t".join(line.decode("utf-8").split("\t")) + "\n")

biomart_response.close()


# +
idmap = pd.read_csv(
    response_file,
    sep="\t",
    names=["Ensembl", "Ensembl_prot", "Entrez ID", "HGNC symbol"],
)

display(idmap.head())
print("Found " + str(len(idmap)) + " genes.")
# -

idmap.dropna(subset=["Ensembl_prot"], inplace=True)
idmap["Entrez ID"] = idmap["Entrez ID"].astype(int)
display(idmap.head())
print("remaining " + str(len(idmap)) + " genes.")

idmap.to_csv(os.path.join(data_path, "Ensembl_prot_entrez.tsv"), index=False, sep="\t")
