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

# ## Parse AD (2022) GWAS
#
# Here we want to convert the SNP ID column in chr/pos/allele format to an rs number.

# +
import os
import pandas as pd
from IPython.display import display

from pyliftover import LiftOver
# -

signif_pairs_path = os.path.join(
    "/home/c.mpnme/scratch-data/gwas/EADB_proxy2_meta.txt"
)
signif_pairs = pd.read_csv(signif_pairs_path, delim_whitespace=True)
display(signif_pairs)

signif_pairs.BEGIN = signif_pairs.BEGIN.astype(int)
signif_pairs.rename(columns = {"BEGIN":"pos"}, inplace=True)
signif_pairs = signif_pairs[["PVALUE", "chr", "pos"]]
display(signif_pairs)

# Generate file to pass to LiftOver to get hg19/GRc37 
lo = LiftOver('../data/gwas/hg38ToHg19.over.chain.gz')

signif_pairs['pos_hg19'] = 0
for index, row in signif_pairs.iterrows():
    result = lo.convert_coordinate('chr' + str(row['chr'].astype(int)), row['pos'].astype(int))
    if result:
        signif_pairs.at[index, "pos_hg19"] = result[0][1]

display(signif_pairs)

signif_pairs[signif_pairs['pos_hg19'] == 0]

# remove unconverted entries
signif_pairs = signif_pairs[~(signif_pairs['pos_hg19'] == 0)]
signif_pairs.drop(['pos'], axis=1, inplace=True)
signif_pairs.rename(columns = {"pos_hg19":"pos"}, inplace=True)

signif_pairs

# Convert chr and pos to an rs number using Kaviar

kaviar_df = pd.read_csv(
    "../data/wrangling/kaviar_parsed.tsv",
    sep="\t",
)

# Convert chr column to string, to allow merge below to complete

# kaviar['chr'] = kaviar['chr'].astype(str)
#kaviar.columns = ["chr", "bp", "snp_id"]
display(kaviar_df.head())

display(type(kaviar_df["pos"][0]))
display(type(signif_pairs["pos"][0]))

joined_df = pd.merge(kaviar_df, signif_pairs, on=["chr", "pos"], how="right")
display(joined_df.head())

len(joined_df[~joined_df['snp_id'].isnull()])

joined_df = joined_df[~joined_df['snp_id'].isnull()]

len(joined_df)

joined_df.isnull().any(axis=1)

joined_df.drop(['chr', 'pos'], axis=1, inplace=True)

joined_df

joined_df.to_csv(
    os.path.join("../data/gwas/", "AD2022.tsv"),
    sep="\t",
    index=False,
)

joined_df[joined_df.duplicated(subset="chr") == True]
