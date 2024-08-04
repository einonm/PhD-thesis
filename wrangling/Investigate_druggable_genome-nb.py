# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.1
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
#     version: 3.7.0
# ---

# # Investigate the overlap of ATC drug targets and the druggable genome

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
The raw code for this jupyter notebook is by default hidden for easier reading.
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

# * Load in ATC level 5 drug lists then 'transpose' so that each row is a gene, listing drugs that target it.

# +
# load in atc5 gene data
atc_file = "../data/drugbank/dbank_gene_set-atc5.txt"

atc_df = pd.read_csv(atc_file, sep="\t", header=None, index_col=0)

atc_genes_df = pd.DataFrame(columns=("gene", "atc_drugs"))
atc_genes_df.set_index("gene", inplace=True)

print("Total ATC target drugs = " + str(atc_df.shape[0]))

for index, row in atc_df.iterrows():
    for gene in row[1].split(" "):
        if gene in atc_genes_df.index:
            atc_genes_df.loc[gene]["atc_drugs"] = (
                atc_genes_df.at[gene, "atc_drugs"] + " " + index
            )
        else:
            atc_genes_df.set_value(gene, "atc_drugs", index)

display(atc_genes_df.head())
print("Total ATC target genes = " + str(atc_genes_df.shape[0]))
# -

# * Merge the ATC gene list with the druggable genome, in order to find the overlap

merged_df = pd.merge(
    atc_genes_df, dg_df, how="outer", left_index=True, right_index=True
)
display(merged_df.head())
print("Total genes, combining ATC + druggable genome = " + str(merged_df.shape[0]))

# * Look at genes that are ATC drug targets and are in the druggable genome (those genes with a non-empty atc_drugs list and also a druggable genome entry)

atc_in_dg = merged_df[
    ~pd.isnull(merged_df["atc_drugs"]) & ~pd.isnull(merged_df["chr_b37"])
]
display(atc_in_dg.head())
print("ATC genes found in druggable genome = " + str(atc_in_dg.shape[0]))

# * Conversely look at those genes that are ATC drug targets, but not in the druggable genome (Entries with a non-empty atc_drugs list but no druggable genome entry)

atc_not_in_dg = merged_df[
    ~pd.isnull(merged_df["atc_drugs"]) & pd.isnull(merged_df["chr_b37"])
]
display(atc_not_in_dg.head())
print("ATC genes not found in druggable genome = " + str(atc_not_in_dg.shape[0]))

# So how many drugs do we throw away if we only look at the druggable genome (i.e. how many drugs have targets that are all missing from the druggable genome?).

# +
missing = atc_not_in_dg.drop(
    ["hgnc_names", "chr_b37", "start_b37", "end_b37", "description"], axis=1
)

missing_drugs_df = pd.DataFrame(columns=("atc_drug", "genes"))
missing_drugs_df.set_index("atc_drug", inplace=True)

for index, row in missing.iterrows():
    for drug in row[0].split(" "):
        if drug in missing_drugs_df.index:
            missing_drugs_df.loc[drug]["genes"] = (
                missing_drugs_df.at[drug, "genes"] + " " + index
            )
        else:
            missing_drugs_df.set_value(drug, "genes", index)

display(missing_drugs_df.head())
print(
    "Potentially there are "
    + str(missing_drugs_df.shape[0])
    + " drugs that could be discarded."
)
# -

# Now list drugs that do not have druggable genome targets:

# +
missing_merged_df = pd.merge(
    atc_df, missing_drugs_df, how="inner", left_index=True, right_index=True
)

# For each drug, sort each gene list and compare (just check lengths of gene lists)
for index, row in missing_merged_df.iterrows():
    if len(row[1]) == len(row["genes"]):
        display(index + " : " + row[1])

# -

# ### So only one drug will be discarded from the analysis, as well as 273 genes.
