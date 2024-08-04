# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Parse eQTL data
#
# Inputs: 
#
#  1. GTEx V7 significant pairs (https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz)
#  2. GTEx variant ID to dbSNP rs ID lookup table (https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz)
#
# Outputs:
#
# Methods:
#
# Here we want to get a list of expression-QTL SNP and gene/protein associations, using the GTEx analysis results as a source.
#
# **The Parse_eQTL_data.sh script automates this method**

import os
import pandas as pd
import sys
from IPython.display import display

# Let's look at the QTL data - b37 used here.
# * N.B. variant ID is chr-pos-ref-alt-build.
# * **WARNING - slope gives the effect slopt w.r.t. ALT allele, not REF!!**

# pval_nominal in the table below is the p value we need to consider (slope and slope se for effect size)

# +
#ip = get_ipython()

#if ip.has_trait('kernel'):
#    signif_pairs_path  = os.path.join('/neurocluster/databank/QTLs/GTEx/v7/',
#                                      'Brain_Anterior_cingulate_cortex_BA24.signifpairs.txt')
#    outfile = 'data/qtl/qtls-all/Brain_Anterior_cingulate_cortex_BA24.gene.annot'
#else:
signif_pairs_path = sys.argv[1]
outfile = sys.argv[2]

signif_pairs = pd.read_csv(signif_pairs_path, sep="\t")

# Remove X chromosome SNPs
signif_pairs = signif_pairs[~signif_pairs['variant_id'].astype(str).str.startswith('X')]

display(signif_pairs.head())

uniq_genes = len(pd.value_counts(signif_pairs['gene_id']))
print("There are " + str(uniq_genes) + " unique genes.")
# -

# Check if any beta values are too low to be be considered (should be 0)

print(len(signif_pairs[abs(signif_pairs['slope']) < 0.10]))

# From the significant SNP/gene pairs, map to rs numbers. GTEx provide a mapping table of rs IDs to GTEx variant_ids.

var_rs = pd.read_csv('/neurocluster/databank/QTLs/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt', sep='\t')
var_rs.head()

# Remove entries that have no rs number
var_rs = var_rs[~var_rs['rs_id_dbSNP147_GRCh37p13'].astype(str).str.startswith('.')]
display(var_rs.describe())

# +
joined_df = pd.merge(var_rs, signif_pairs, on='variant_id')
display(joined_df.head())

print(str(len(pd.value_counts(joined_df['gene_id']))))
# -

joined_df.drop(['chr','variant_pos'], axis=1, inplace=True)
joined_df.set_index(['rs_id_dbSNP147_GRCh37p13'], inplace=True)
joined_df.head()

# ## Create MAGMA annotation file
#
# We'd now like to create an annotation file (mapping SNPs to genes) to use with MAGMA. By default, MAGMA uses genome proximity of SNPs to genes to do this. We are using sQTL data, so we have to manually create the annotation file.
#
# An annotation file is a set of lines, each of the format:
#
# <tt>gene_name   chr:start_pos:stop_pos   snp1 snp2 snp3 snp4...etc</tt>
#
# Where we've chosen Ensembl IDs to identify genes.
#
# So let's list, for each unique Ensembl ID, the list of snps associated with it:

# +
annot = pd.DataFrame(columns=('Ensembl', 'position', 'snplist'))
annot.set_index('Ensembl', inplace=True)

for index, row in joined_df.iterrows():
    if row['gene_id'] in annot.index:
        annot.loc[row['gene_id']]['snplist'] = annot.at[row['gene_id'], 'snplist'] + ' ' + index
    else:
        annot.at[row['gene_id'], 'snplist'] = index

display(annot.head())
annot_pos = annot.copy()

# +
# Remove version numbers from Ensembl Id
annot_pos.index = annot_pos.index.to_series().str.replace(r'\..*','')
display(annot_pos.head())

print(len(annot_pos))

# +
gene_pos = pd.read_csv('../data/wrangling/NCBI37.3.ensembl.gene.loc', sep='\t', index_col=0,
                       names=['Ensembl','Entrez', 'chrom', 'chromStart', 'chromEnd', 'HGNC'])

display(gene_pos.head())
# We're only interested in Ensembl / chr / pos
gene_pos = gene_pos[['chrom', 'chromStart', 'chromEnd']]

display(gene_pos.head())
# -

annot_mrg = pd.merge(gene_pos, annot_pos, how='right', left_index=True, right_index=True)
annot_mrg.dropna(subset=['chrom'], inplace=True)
annot_mrg['chromStart'] = annot_mrg['chromStart'].astype(int)
annot_mrg['chromEnd'] = annot_mrg['chromEnd'].astype(int)
display(annot_mrg.head())
print("After merge, there are " + str(len(annot_mrg)) + " genes.")

# +
for index, row in annot_mrg.iterrows():
        annot_mrg.at[index, 'position'] = str(row['chrom']) + ":" +
                                           str(row['chromStart']) + ":" +
                                           str(row['chromEnd']))

annot_final = annot_mrg[['position', 'snplist']]
annot_final.drop_duplicates(subset='snplist', inplace=True)

display(annot_final.head())
print("There are " + str(len(annot_final)) + " genes for analysis, from " + str(uniq_genes) + " initially")
# -

annot_final.to_csv(outfile, sep = "\t", header=False)


