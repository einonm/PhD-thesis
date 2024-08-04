# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# We want to transform the list of drug/target interactions into a list of pubchem CIDs vs Entrez IDs, with ATC codes if possible.

import numpy as np
import pandas as pd

df = pd.read_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/stitch/9606.protein_chemical.links.detailed.v5.0.tsv",
    sep="\t",
)

display(df.shape)
df.head()

# Remove any predicted drug targets, leaving experimental and database scores only.
noexp_df = df[(df["experimental"] > 0) | (df["database"] > 0)]
noexp_df.shape

noexp_df[
    (noexp_df["database"] > 0)
    & (noexp_df["experimental"] > 0)
    & (noexp_df["combined_score"] < 700)
].head(100)

# Set an acceptable threshold to scores. This combined score also
# takes into account textmining and prediction , so use experimental and database scores below instead
# thresh_df = noexp_df[noexp_df['combined_score'] > 700]
thresh_df = noexp_df[(noexp_df["experimental"] > 700) | (noexp_df["database"] > 700)]
thresh_df.shape

# +
thresh_df.drop(
    ["experimental", "prediction", "database", "textmining", "combined_score"],
    axis=1,
    inplace=True,
)
display(len(thresh_df["chemical"].unique()))

thresh_df.protein = thresh_df.protein.str[5:]
display(thresh_df.head())
# -

thresh_df[thresh_df['chemical'] == 'CIDm00002767'] # (maps to 2770, which is cisplatin?)

# ### Convert entrez to ensembl using the Biomart conversion table

conversion_df = pd.read_csv(
    "/home/mpnme/data/phd/wrangling/Ensembl_prot_entrez.tsv", sep="\t"
)
conversion_df.head()

# +
 ## START HERE!! 8370 may introduce a conversion error?
# ENSP00000289352, maps to ENSG00000158406 / HIST1H4H is shared

conversion_df[conversion_df['Entrez ID'] == 8370.0]
# -

conv_df = pd.merge(
    thresh_df, conversion_df, left_on="protein", right_on="Ensembl_prot", how="left"
)
conv_df.drop(
    ["Ensembl", "protein", "Ensembl_prot", "HGNC symbol"], axis=1, inplace=True
)
display(conv_df.shape)
conv_df.dropna(subset=["Entrez ID"], inplace=True)
conv_df["chemical"] = conv_df.chemical.str[4:]
conv_df["Entrez ID"] = conv_df["Entrez ID"].astype(int)
display(conv_df.head())
display(conv_df.shape)

drug_df = pd.DataFrame(sorted(conv_df["chemical"].unique()))
drug_df["gene_list"] = np.empty((len(drug_df), 0)).tolist()
drug_df["gene_count"] = 0
drug_df.set_index(0, inplace=True)
drug_df.index.name = "PubChem_CID"

# +
for index, row in conv_df.iterrows():
    drug_df.at[row["chemical"], "gene_list"].append((int(row["Entrez ID"])))

drug_df["gene_list"] = [list(set(x)) for x in drug_df["gene_list"]]
drug_df["gene_list"] = [sorted(l) for l in drug_df["gene_list"]]
# -

drug_df.head()

drug_df["gene_count"] = drug_df["gene_list"].str.len()
drug_df["gene_list"] = [" ".join(str(l) for l in x) for x in drug_df["gene_list"]]

drug_df.gene_count.sum()

display(len(drug_df.index))
display(len(drug_df[drug_df["gene_count"] >= 5]))
display(drug_df.head())

list(drug_df[drug_df.index == '00002767']['gene_list'])

drug_df = drug_df.drop(["gene_count"], axis=1)
drug_df.to_csv("../data/target-dbs/stitch_targets2.csv", sep="\t")
