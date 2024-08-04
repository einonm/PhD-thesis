# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Parse CMC eQTL data
#
# Here we want to get a list of expression-QTL SNP and gene/protein associations, using the Common Mind Consortium analysis results as a source.

import os
import pandas as pd
from IPython.display import display

# pval_nominal in the table below is the p value we need to consider (slope and slope se for effect size)

# +
cmc_path = os.path.join('/neurocluster/databank/QTLs/common_mind',
                                 'CMC_MSSM-Penn-Pitt_DLPFC_mRNA_eQTL-adjustedSVA-binned.formated.txt')
cmc = pd.read_csv(cmc_path, sep=" ", index_col=0)[['Gene']]
display(cmc.head())

uniq_genes = len(pd.value_counts(cmc['Gene']))
print("There are " + str(uniq_genes) + " unique genes.")
# -

# Remove X chromosome snps
cmc = cmc[cmc.index.str.startswith('r')]

# ## Create MAGMA annotation file
#
# We'd now like to create an annotation file (mapping SNPs to genes) to use with MAGMA. By default, MAGMA uses genome proximity of SNPs to genes to do this.
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

for index, row in cmc.iterrows():
    if row['Gene'] in annot.index:
        annot.loc[row['Gene']]['snplist'] = annot.at[row['Gene'], 'snplist'] + ' ' + index
    else:   
        annot.at[row['Gene'], 'snplist'] = index
    
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

# remove X chromosome genes
annot_mrg = annot_mrg[~(annot_mrg['chrom'] == 'X')]

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

annot_final.to_csv(os.path.join('../data/qtl/qtls-all/', 'cmc_eqtl.gene.annot'), sep = "\t", header=False)


