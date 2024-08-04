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

# # Parse sQTL data
#
# Here we want to get a list of splice-QTL SNP and gene associations, in the form of a MAGMA .gene.annot annotaion file.

import os
import pandas as pd
from IPython.display import display
#!pip install scipy
import scipy.stats as stats
import warnings
warnings.filterwarnings('always')

# ## sQTL data files
# Let's look at the splice QTL data (from a paper, DOI 10.1038/ncomms14519) - this has already been filtered to remove all Double corrected P vals > 10^-8, and is mapped to hg19/GRCh37.

# original QTL file
orig_splice_path = os.path.join('data', 'Splice.qtls.tsv')
orig_splice = pd.read_csv(orig_splice_path, sep="\t")
orig_splice.head()

# Extract colums we're interested in - rs SNP ID, SNP position, Gene. 

# +
splice = orig_splice[['sQTL_SNP_ID', 'SNP_position_hg19', 'Gene', 'Double_corrected_P']]

splice.set_index(['Gene'], inplace=True)

display(splice.head(10))

uniq_genes = len(pd.value_counts(splice.index))

print("There are " + str(uniq_genes) + " unique genes with associated sQTL data")
# -

# Look at each SNP in the splice table and convert the official gene name to tthe corresponding Ensembl gene ID 

#Look at each SNP in the splice table and convert the official gene name to tthe corresponding Ensembl gene ID 
gene_pos = pd.read_csv(os.path.join('data', 'Ensembl_pos.tsv'), sep='\t', index_col=2)
display(gene_pos.head())

# +
splice_conv = pd.merge(gene_pos, splice, how='inner', left_index=True, right_index=True)

splice_conv.set_index(['Ensembl'], inplace=True)
splice_conv.drop(['Entrez ID'], inplace=True)

display(splice_conv.head())
print("After conversion, " + str(len(pd.value_counts(splice_conv.index))) + 
      " genes remaining from " + str(uniq_genes))
# -

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

for index, row in splice_conv.iterrows():
    if index in annot.index:
        annot.loc[index]['snplist'] = annot.at[index, 'snplist'] + ' ' + row['sQTL_SNP_ID']
    else:  
        annot.at[index, 'position'] = str(row['chrom']) + ":" +
                                           str(row['chromStart']) + ":" +
                                           str(row['chromEnd'])
        annot.at[index, 'snplist'] = row['sQTL_SNP_ID']
    
display(annot.head())

print("After conversion, " + str(len(pd.value_counts(annot.index))) + 
      " SNPS remaining from " + str(uniq_genes))
# -

annot.to_csv(os.path.join(os.path.pardir, 'magma_analysis', 'qtls', 'sqtl.gene.annot'), sep = "\t", header=False)


