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

# ## Parse geneatlas hypertension (I10)
#
# Here we want to convert the SNP ID column in chr/pos/allele format to an rs number.

import os
import pandas as pd
from IPython.display import display

signif_pairs_path = os.path.join('../data/gwas/PD-23am/nallsEtAl2019_excluding23andMe_allVariants.tab')
signif_pairs = pd.read_csv(signif_pairs_path, sep="\t")
display(signif_pairs.head())

# Extract the MarkerName to chr / bp columns

tester_df = signif_pairs.head()
signif_pairs['chr'] = signif_pairs['SNP'].str.split(':').str[0].str[3:]
signif_pairs['bp'] = signif_pairs['SNP'].str.split(':').str[1]
signif_pairs.head()

signif_pairs['chr'] = pd.to_numeric(signif_pairs['chr'])
signif_pairs['bp'] = pd.to_numeric(signif_pairs['bp'])
display(type(signif_pairs['chr'][0]))
display(type(signif_pairs['bp'][0]))

signif_pairs.count()

# Convert chr and pos to an rs number using Kaviar

kaviar = pd.read_csv('../data/wrangling/kaviar_parsed.tsv', sep='\t')

kaviar

# Convert chr column to string, to allow merge below to complete

#kaviar['chr'] = kaviar['chr'].astype(str)
kaviar.columns = ['chr','bp', 'snp_id']
display(kaviar.head())

display(type(kaviar['chr'][0]))
display(type(kaviar['bp'][0]))

joined_df = pd.merge(kaviar, signif_pairs, on=['chr', 'bp'])
display(joined_df.head())

joined_df.isnull().any(axis=1)

joined_df.to_csv(os.path.join('../data/gwas/PD-23am/', 'METAANALYSIS_PdGene-rs.tbl'), sep = "\t", index=False)

joined_df[joined_df.duplicated(subset='chr') == True]



