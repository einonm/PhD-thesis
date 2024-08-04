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

# # Parse mQTL data

import os
import pandas as pd
from IPython.display import display

# Read an mQTL file, and use just the columns we're interested in - Gene, SNP and FDR value.
#
# Currently, there are 3 mQTL files used - CP, FC, SN, for 3 parts of the brain. 

# +
# FC, CP, SN
mqtl_type = 'SN'

# Read just the columns we're interested in
mqtl_path = os.path.join('../data/qtl/exonQTLs_mQTLs/', 
                                 'Annotated.' + mqtl_type + '.Cis.mQTL.FDR.April.2017.txt')
# column 8 is 'CPG_GENE', set it as the index
mqtl = pd.read_csv(mqtl_path, sep=" ", index_col=8)[['SNP', 'FDR', 'BETA']]
display(mqtl.head())
print("There are " + str(len(mqtl)) + " SNPs given.")

# -

# Just look at snps that have BETA > 0.10
# REMOVED TO REPLICATE LEON'S RESULTS
signif_mqtl = mqtl.copy()
#signif_mqtl = mqtl[mqtl['BETA'] > 0.10]
#display(signif_mqtl.head())
uniq_genes = len(pd.value_counts(signif_mqtl.index))
#print("After filtering BETA > 0.10, there are " + str(len(signif_mqtl)) + 
#      " SNPs remaining, linked to " + str(uniq_genes) + " genes")

# Out of interest, how many snps with BETA > 0.10 have an FDR > 1E-8?
print(len(signif_mqtl[signif_mqtl['FDR'] > 1E-8]))
signif_mqtl[signif_mqtl['FDR'] > 1E-8].head()

#Look at each SNP in the splice table and convert the official gene name to the corresponding Ensembl gene ID 
gene_pos = pd.read_csv('../data/wrangling/NCBI37.3.ensembl.gene.loc', sep='\t', 
                       index_col=5, header=None, names=['Ensembl','1','chrom','chromStart','chromEnd', 'HGNC'])[['Ensembl','chrom','chromStart','chromEnd']]
display(gene_pos.head())

# +
mqtl_conv = pd.merge(gene_pos, signif_mqtl, how='right', left_index=True, right_index=True)

mqtl_conv.dropna(inplace=True)
mqtl_conv.set_index(['Ensembl'], inplace=True)

mqtl_conv['chromStart'] = mqtl_conv['chromStart'].astype(int)
mqtl_conv['chromEnd'] = mqtl_conv['chromEnd'].astype(int)
display(mqtl_conv.head())
print("After conversion to Ensembl, there are " + str(len(pd.value_counts(mqtl_conv.index))) + 
      " genes remaining from " + str(uniq_genes))

# +
# Load the list of MHC genes, in order to remove them from the gene set
# REMOVED TO REPLICATE LEON'S RESULTS

#mhc_genes = pd.read_csv('data/mhc_genes.txt', sep='\t')['Ensembl']
#mhc_genes = mhc_genes.unique()
#mqtl_conv = mqtl_conv.loc[~mqtl_conv.index.isin(mhc_genes)]
#display(mqtl_conv.head())
#print("After removing MHC genes, there are " + str(len(pd.value_counts(mqtl_conv.index))) + 
#      " genes remaining from " + str(uniq_genes))

# +
# Create MAGMA annotation file
annot = pd.DataFrame(columns=('Ensembl', 'position', 'snplist'))
annot.set_index('Ensembl', inplace=True)

for index, row in mqtl_conv.iterrows():
    if index in annot.index:
        annot.loc[index]['snplist'] = annot.at[index, 'snplist'] + ' ' + row['SNP']
    else:  
        annot.at[index, 'position'] = str(row['chrom']) + ":" +
                                           str(row['chromStart']) + ":" +
                                           str(row['chromEnd'])
        annot.at[index, 'snplist'] = row['SNP']
    
display(annot.head())

print("After creating the annotation list, there are " + str(len(pd.value_counts(annot.index))) + 
      " genes remaining from " + str(uniq_genes))
# -

annot.to_csv(os.path.join('../data/qtl/qtls-all/', mqtl_type + '_mqtl.gene.annot'), 
             sep = "\t", header=False)


