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

# We want to transform the list of drug/target interactions into a list of pubchem CIDs vs Entrez IDs, with ATC codes if possible. Also include the evidence presented for the interaction (e.g. PMID).

import numpy as np
import pandas as pd

# +
# df = pd.read_csv("/home/mpnme/source/genomic-novel-targets/data/target-dbs/GtoPdb/gtopdb_ligands.tsv", sep='\t')

# +
# display(df.shape)
# display(df.columns)
# df.head(100)
# -

df2 = pd.read_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/GtoPdb/gtopdb_interactions.tsv",
    sep="\t",
)

display(df2.shape)
display(len(df2.columns))
display(df2.columns)
display(df2["ligand_species"].unique())
display(len(df2[df2["target_species"] == "Human"]))
df2.head()

display(df2.shape)
df2 = df2[df2["target_species"] == "Human"]
df3 = df2.drop(
    [
        "target",
        "target_id",
        "target_gene_symbol",
        "target_ligand",
        "target_ligand_id",
        "target_ensembl_gene_id",
        "target_ligand_gene_symbol",
        "target_ligand_ensembl_gene_id",
        "target_ligand_uniprot",
        "target_ligand_pubchem_sid",
        "ligand_species",
        "ligand_gene_symbol",
        "affinity_low",
        "affinity_median",
        "affinity_high",
        "affinity_units",
        "original_affinity_units",
        "concentration_range",
        "ligand_context",
        "endogenous",
        "original_affinity_low_nm",
        "original_affinity_median_nm",
        "original_affinity_high_nm",
        "original_affinity_relation",
        "action_comment",
        "selectivity",
        "action",
        "receptor_site",
        "assay_description",
        "type",
        "primary_target",
        "ligand",
        "ligand_id",
        "target_species",
    ],
    axis=1,
)
display(df3.shape)
df3 = df3.dropna()
display(df3.shape)
df3["ligand_pubchem_sid"] = df3["ligand_pubchem_sid"].astype(int)
display(df3.head())

# +
# Some gene ids are a collection of gene IDs with '\' between, remove these for now
df4 = df3[df3["target_uniprot"].apply(lambda x: x.find("|") == -1)]
display(df4.shape)

# Look at evidence for the drug amlodipine
df4[df4["ligand_pubchem_sid"] == 178103560]
# -

pubchem_sid_series = pd.DataFrame(df4["ligand_pubchem_sid"].unique())
pubchem_sid_series.to_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/GtoPdb/pubchem_sid.list",
    header=None,
    index=None,
)

# **Submit list as a job to https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi. Convert from pubchem SID to pubchem CID.**

pubchem_cid = pd.read_csv(
    "/home/mpnme/source/genomic-novel-targets/data/target-dbs/GtoPdb/293344786362584733.txt",
    sep="\t",
    header=None,
    names=["ligand_pubchem_sid", "pubchem_cid"],
)

display(len(pubchem_cid))
pubchem_cid.dropna(inplace=True)
pubchem_cid["pubchem_cid"] = pubchem_cid["pubchem_cid"].astype(int)
display(len(pubchem_cid))
pubchem_cid.head()

result_df = pd.merge(df4, pubchem_cid, how="left", on="ligand_pubchem_sid")


display(len(result_df))
result_df.dropna(subset=["pubchem_cid"], inplace=True)
display(len(result_df))
result_df["pubchem_cid"] = result_df["pubchem_cid"].astype(int)
result_df.drop(["ligand_pubchem_sid"], axis=1, inplace=True)
display(result_df.head())

# +
data_path = "../data/drugbank/"
idmap_path = data_path + "HUMAN_9606_idmapping_selected.tab"
idmap_full = pd.read_csv(
    idmap_path,
    index_col=0,
    sep="\t",
    names=[
        "UniProtKB-AC",
        "UniProtKB-ID",
        "(EntrezGene)",
        "RefSeq",
        "GI",
        "PDB",
        "GO",
        "UniRef100",
        "UniRef90",
        "UniRef50",
        "UniParc",
        "PIR",
        "NCBI-taxon",
        "MIM",
        "UniGene",
        "PubMed",
        "EMBL",
        "EMBL-CDS",
        "Ensembl",
        "Ensembl_TRS",
        "Ensembl_PRO",
        "Additional PubMed",
    ],
)

idmap = idmap_full[["Ensembl", "(EntrezGene)"]]
idmap.rename(columns={"(EntrezGene)": "EntrezGene"}, inplace=True)

# Remove multiple gene IDs if they exist
# TODO - why are there > 1 genes listed for some proteins?
for index, row in idmap.iterrows():
    if str(row["EntrezGene"]).find(";") != -1:
        idmap.at[index, "EntrezGene"] = row["EntrezGene"][
            : row["EntrezGene"].index(";")
        ]
    if str(row["Ensembl"]).find(";") != -1:
        idmap.at[index, "Ensembl"] = row["Ensembl"][: row["Ensembl"].index(";")]

display(idmap.head())
# -

result2_df = pd.merge(
    result_df, idmap, how="left", left_on="target_uniprot", right_index=True
)
display(result2_df.shape)
result2_df.drop(["target_uniprot"], axis=1, inplace=True)
result2_df.dropna(inplace=True)
display(result2_df.shape)
result2_df.head()

# save the result2_df table in a seperate file, used to track evidence for a DGI
result2_df.to_csv(
    "../data/target-dbs/gtopdb_dgi_evidence.csv",
    sep="\t",
    index=None,
    columns=["pubmed_id", "pubchem_cid", "Ensembl"],
)

# +
drug_df = pd.DataFrame(columns=["gene_list", "gene_count"])
drug_df.index.name = "drug_pubchem_cid"

for index, row in result2_df.iterrows():
    drug = row["pubchem_cid"]
    if drug not in drug_df.index.values:
        drug_df.at[drug, "gene_list"] = [row["EntrezGene"]]
    else:
        if row["EntrezGene"] not in drug_df.at[drug, "gene_list"]:
            drug_df.at[drug, "gene_list"].append(row["EntrezGene"])

drug_df["gene_list"] = [sorted(l) for l in drug_df["gene_list"]]

display(drug_df.head(100))
# -

drug_df["gene_count"] = drug_df["gene_list"].str.len()
drug_df["gene_list"] = [" ".join(str(l) for l in x) for x in drug_df["gene_list"]]

display(len(drug_df.index.unique()))
display(len(drug_df[drug_df["gene_count"] >= 5]))
display(drug_df.head(100))

# merge IDs back


drug_df.drop(["gene_count"], axis=1, inplace=True)
drug_df.head(50)

drug_df.to_csv("../data/target-dbs/gtopdb_targets.csv", sep="\t")

drug_df
