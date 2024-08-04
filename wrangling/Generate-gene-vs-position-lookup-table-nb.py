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

# ## Generate gene vs position table
#
# Using tables obtained from online databases, create a table or Ensemble / Gene name / Entrez ID vs chromosome and position for use in data conversions

import os
import pandas as pd
from IPython.display import display
import re

# #!pip install biomart
from biomart import BiomartServer

# +
server = BiomartServer("http://grch37.ensembl.org/biomart")
data_path = "../data/wrangling/"
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
            "entrezgene_id",
            "hgnc_symbol",
            "chromosome_name",
            "start_position",
            "end_position",
            "version",
        ]
    }
)

biomart_response = open(response_file, "w")

for line in response.iter_lines():
    biomart_response.write("\t".join(line.decode("utf-8").split("\t")) + "\n")

biomart_response.close()


# +
idmap_full = pd.read_csv(
    response_file,
    sep="\t",
    names=[
        "Ensembl",
        "Entrez ID",
        "HGNC symbol",
        "chrom",
        "chromStart",
        "chromEnd",
        "version",
    ],
)

idmap = idmap_full[pd.notnull(idmap_full["Entrez ID"])]

display(idmap.head())
print("Found " + str(len(idmap)) + " genes.")

# +
# TODO - strategy to drop dupes here may not be most efficient in using genes in magma list?
# idmap.sort_values(by=['Entrez ID', 'version'])
idmap_dedup = idmap  # .drop_duplicates(subset='Entrez ID', keep='first')
# display(idmap_dedup.head())
# print(len(idmap_dedup))

# Look at this! multiple Entrez IDs for one HGNC/Ensembl code - wtf??
# display(idmap_dedup[idmap_dedup['HGNC symbol'] == 'OR4F29'])
# -

idmap_dedup.to_csv(os.path.join(data_path, "Ensembl_pos.tsv"), index=False, sep="\t")

# +
# look at how many have been discarded
merged = idmap.merge(idmap_dedup, indicator=True, how="left")
display(merged.head())

len(merged[merged["_merge"] != "both"])
