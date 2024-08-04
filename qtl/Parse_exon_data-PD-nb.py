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

# ## Parse exon data for PD GWAS
#
# Here we want to get a list of exon and gene/protein associations
#
# **TODO - make this generic for any GWAS (don't limit exon data to a specific set of GWAS SNPs)**

import os
import pandas as pd
from IPython.display import display

# Load the GWAS summary results
gwas_path = os.path.join('/home', 'mpnme', 'data', 'GWAS', 
                            'PD_2014_META.txt')
gwas = pd.read_csv(gwas_path, sep="\t")[['MarkerName', 'Chr', 'Bp']]
#display(gwas.head())

# Load the position & gene of all exon table
exon_data_path = os.path.join('/home', 'mpnme', 'data', 'QTLs', 
                            'EXON_Transcript_Data_SUSHI.txt')
exon_data = pd.read_csv(exon_data_path, delim_whitespace=True, 
                        header=None, names=['chr', 'start_bp', 'end_bp', 'gene', 'strand', 'exon'])
#display(exon_data.head())

# +
# remove 'chr' from all chromosome number entries
exon_data['chr'] = pd.to_numeric(exon_data['chr'].str.replace(r'chr', ''), 
                             errors='coerce')

# look at exon data for chr1 to chr22 only
# This removes all entries where exon_data[0] == NaN
exon_data = exon_data[exon_data['chr'] == exon_data['chr']]
exon_data['chr'] = pd.to_numeric(exon_data['chr'], downcast='integer')

#display(exon_data.head())

# +
def get_annotations(chrm):
    import pandas as pd
    annot = pd.DataFrame(columns=('snplist', 'Gene'))
    annot.set_index('Gene', inplace=True)
    exon_data_chr = exon_data[exon_data['chr'] == chrm].sort_values('start_bp')
    exon_data_chr.set_index(['start_bp'], inplace=True)
    gwas_chr = gwas[gwas['Chr'] == chrm]
    gwas_chr.set_index(['Bp'], inplace=True)

    for index, row in gwas_chr.iterrows():
        matches = exon_data_chr.query('@index >= start_bp & @index <= end_bp')
        # In the case of multiple matches, assume they are all the same gene
        # TODO - check this?
        if not matches.empty:
            annot = annot.append({'snplist' : row['MarkerName'], 'Gene' : matches.iloc[0]['gene']}, ignore_index=True)

    return annot
        

#annot.to_csv(os.path.join(os.path.pardir, 'magma_analysis', 'data', 'PD_exon.gene.annot'), sep = "\t", header=False)

# +
import ipyparallel

clients = ipyparallel.Client()
dview = clients[:]

# Put a copy of our source data frames on all the clients
dview.push({'exon_data':exon_data})
dview.push({'gwas':gwas})

results = dview.map_sync(get_annotations, range(1,23))
# -

# Combine all results into one table
exons = pd.concat(results)
exons.set_index('Gene', inplace=True)
exons.head()

#Look at each SNP in the splice table and convert the official gene name to the corresponding Ensembl gene ID 
gene_pos = pd.read_csv(os.path.join('data', 'Ensembl_pos.tsv'), sep='\t', index_col=2)
display(gene_pos.head())

# +
exon_conv = pd.merge(gene_pos, exons, how='inner', left_index=True, right_index=True)

exon_conv.set_index(['snplist'], inplace=True)

display(exon_conv.head())
#print("After conversion to Ensembl, there are " + str(len(pd.value_counts(mqtl_conv.index))) + 
#      " genes remaining from " + str(uniq_genes))
# -

# From the significant SNP/gene pairs, split out the variant ID to chr/position and map to rs numbers. The [Kaviar database](http://db.systemsbiology.net/kaviar/) has been used to get a SNP rs number / position mapping (hg19).

# +
# For exon snps with chr/pos positions
exon_conv_chr = exon_conv.copy()
# For exon snps with rs positions
exon_conv_rs = exon_conv.copy()

exon_conv_chr['chr'] = ''
exon_conv_chr['pos'] = ''

for index, row in exon_conv_chr.iterrows():
    if index.startswith('chr'):
        exon_conv_chr.at[index, 'chr'] = index.split(':')[0][3:]
        exon_conv_chr.at[index, 'pos'] = index.split(':')[1]
        exon_conv_rs.drop(index, inplace=True)
    else:
        exon_conv_chr.drop(index, inplace=True)
        
display(exon_conv_chr.head())
display(exon_conv_rs.head())
# -

# Convert chr and pos to an rs number using Kaviar
kaviar = pd.read_csv('../kaviar/data/kaviar_parsed.tsv', sep='\t', dtype=object)

# +
# convert chr column to string, to allow merge below to complete
kaviar['chr'] = kaviar['chr'].astype(str)

display(kaviar.head())

# +
joined_df = pd.merge(kaviar, exon_conv_chr)

display(joined_df.head())
# -

joined_df.drop(['chr','pos'], axis=1, inplace=True)
joined_df.set_index(['snp_id'], inplace=True)
joined_df.head()

# Now join the converted exon_chr rows with the exon_rs rows
# Look out for duplicate rs nums?
all_df = pd.concat([joined_df, exon_conv_rs])
dup = all_df.index.duplicated()
all_df.head()

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
    if row['Ensembl'] in annot.index:
        annot.loc[row['Ensembl']]['snplist'] = annot.at[row['Ensembl'], 'snplist'] + ' ' + index
    else:   
        annot.at[row['Ensembl'], 'snplist'] = index

display(annot)

# +
gene_pos = pd.read_csv(os.path.join('data', 'Ensembl_pos.tsv'), sep='\t', index_col=0)

# We're only interested in Ensembl / chr / pos
gene_pos = gene_pos[['chrom', 'chromStart', 'chromEnd']]

gene_pos.head()
# -

annot_mrg = pd.merge(gene_pos, annot, how='inner', left_index=True, right_index=True)
display(annot_mrg.head())

# +
for index, row in annot_mrg.iterrows():
        annot_mrg.at[index, 'position'] = str(row['chrom']) + ":" +
                                           str(row['chromStart']) + ":" +
                                           str(row['chromEnd'])
        
annot_final = annot_mrg[['position', 'snplist']]
annot_final.drop_duplicates(inplace=True)
        
display(annot_final.head())
#print("There are " + str(len(annot_final)) + " genes for analysis, from " + str(uniq_genes) + " initially")
# -

annot_final.to_csv(os.path.join(os.path.pardir, 'magma_analysis', 'qtls', 'exon.gene.annot'), sep = "\t", header=False)


