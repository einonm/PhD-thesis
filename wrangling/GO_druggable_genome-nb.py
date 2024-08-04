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

# # Generate a gene ontology of only druggable genome pathways
#
# We'd like to run a gene set analysis between:
#
# * functional/hybrid mapped disease GWAS genes (with p-values)
#
# and
# * pathways that have at least one gene in the druggable genome
#
# So below we create the set of GO pathways that involve the druggable genome.
#

import os
import sys
import pandas as pd
import numpy as np
from IPython.display import display, HTML

HTML(
    """<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>."""
)

# * Load in the druggable genome, taken from the supplimentary information of https://doi.org/10.1126%2Fscitranslmed.aag1166

# +
# load in druggable genome data
dg_file = "~/data/druggable_genome.csv"

dg_df = pd.read_csv(dg_file, sep="\t", index_col=0)
dg_df.drop(
    [
        "druggability_tier",
        "strand",
        "no_of_gwas_regions",
        "small_mol_druggable",
        "bio_druggable",
        "adme_gene",
    ],
    axis=1,
    inplace=True,
)
display(dg_df.head())
print("Total genes in druggable genome = " + str(dg_df.shape[0]))
# -

# * Load in a GO derived set of pathways (from Janet's / Andrew's annotation pipeline)

# +
go_file = "/home/mpnme/data/GO_ALL_PC_20-2000_MAGMA.LH.Magma.Format.Final.txt"

go_df = pd.read_csv(go_file, sep=" ", header=None)

display(go_df.head())
print("Size of GO table = " + str(go_df.shape[0]))
# -

# Convert the Entrez IDs to Ensembl, and remove any duplicates (as more than one Ensembl ID maps to an Entrez ID)

# +
# Convert Entrez IDs to Ensembl
rosetta_file = "/home/mpnme/data/phd/wrangling/NCBI37.3.ensembl.gene.loc"
rosetta_df = pd.read_csv(
    rosetta_file,
    sep="\t",
    header=None,
    names=["ensembl", "entrez", "chr", "start", "end", "gene"],
)
display(rosetta_df.head())
print("Size of rosetta table = " + str(rosetta_df.shape[0]))

# merged table is bigger than original GO table as more than one Ensembl ID maps to aan Entrez ID.
merge_df = pd.merge(go_df, rosetta_df, how="inner", right_on="entrez", left_on=1)
print("Size of merged GO table = " + str(merge_df.shape[0]))

merge_df.drop(["entrez", "chr", "start", "end", "gene"], axis=1, inplace=True)
merge_df.drop_duplicates(subset=[0, 1], inplace=True)
merge_df.drop([1], axis=1, inplace=True)
display(merge_df.head())
print("Size of merged GO table after duplicates dropped = " + str(merge_df.shape[0]))
# -

# 'Transpose' the table to give one entry per GO term, each listing the genes involved in the pathway.

# +
result = merge_df.groupby(by=0)
counts = result.count()

final = result.apply(lambda x: " ".join(x["ensembl"]))
final_df = final.to_frame()
final_df.rename
final_df.set_index(0)

go_long_df = pd.merge(counts, final_df, left_index=True, right_index=True, how="inner")
go_long_df.index.names = ["go"]
go_long_df.columns = ["gene_count", "gene_list"]
display(go_long_df.head())

# -

# Now remove any table entries which do not have at least one druggable gene.

# +
for index, row in go_long_df.iterrows():
    # look at ensemble gene list, make into list
    # Also order list, then write back to df
    genes = row["gene_list"].split(" ")
    genes.sort()
    go_long_df.set_value(index, "gene_list", " ".join(genes))
    for gene in genes:
        # is gene in druggable genome?
        if gene in dg_df.index:
            # yes - go to next index
            go_long_df.set_value(index, "in_dg", 1)
            break
        else:
            # no - are we the last list item?
            if gene == genes[-1]:
                # yes - remove from list
                go_long_df.set_value(index, "in_dg", 0)

display(go_long_df.head())
print(
    "Size of GO table after non-druggable-genome pathways dropped = "
    + str(go_long_df[go_long_df["in_dg"] == 1].shape[0])
)
print(
    "Number of non-druggable-genome pathways dropped = "
    + str(go_long_df[go_long_df["in_dg"] == 0].shape[0])
)

# +
# drop rows without druggable genome genes
druggable_df = go_long_df[go_long_df["in_dg"] == 1]
display(druggable_df.shape)

# strip gene count and in_dg from df
stripped_df = druggable_df.drop(["gene_count", "in_dg"], axis=1)
stripped_df.shape

# show duplicates
# stripped_df[stripped_df.duplicated(subset='gene_list') == True]
# stripped_df.drop_duplicates(inplace=True)
# new_size = stripped.shape
# -

stripped_df.to_csv(
    "/home/mpnme/data/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt",
    header=False,
    sep="\t",
)
