# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Gene Annotation of SNPs to Ensembl
#
# This method takes a positional list of SNPs from either the 1000 genomes project or from a GWAS summary results set (assumed aligned to build GRCh37), filters out those with |MAF| <= 0.01 and uses it to annotate a list of GRCh37 genes (magma standard set). These gene IDs are then converted from Entrez to Ensembl.
#
# The resultant file can be used to perform a MAGMA gene analysis, associating genes with a disease.

import pandas as pd
import sys
from IPython.display import display
import subprocess as sproc

# snploc_file is the MAF filtered list of SNPS vs position
snploc_file = "../data/gwas/gene_atlas_snps/gene_atlas.loc"
# gene_annot_file is the list of genes annotated by snps, positionally
gene_annot_prefix = "../data/gwas/gene_atlas_snps/magma_gene_atlas.prox"
gene_annot_file = gene_annot_prefix + ".genes.annot"

# Magma takes a snp-loc file with 3 columns - SNP ID, chromosome, position.
# We first have to filter out low MAF values
unfiltered_gene_atlas_snps = pd.read_csv(
    "../data/gwas/gene_atlas_snps/snps.pos.maf.tsv", sep="\t"
)

display(len(unfiltered_gene_atlas_snps))
display(unfiltered_gene_atlas_snps.columns)

gene_atlas_snps = unfiltered_gene_atlas_snps[unfiltered_gene_atlas_snps["MAF"] >= 0.01]
gene_atlas_snps = gene_atlas_snps[gene_atlas_snps["SNP"].str.startswith("rs")]

display(gene_atlas_snps.describe())
display(len(gene_atlas_snps))

# Drop MAF values to get the SNP-loc files suitable for a magma annotation
gene_atlas_loc = gene_atlas_snps.drop(["MAF"], axis=1)

gene_atlas_loc.to_csv(snploc_file, header=None, sep="\t", index=None)

result = sproc.check_output(
    ["./create_magma_prox_annot.sh", snploc_file, gene_annot_prefix]
)

print(str(result, "ascii"))

# As the generated gene annot file has a variable number of columns (due to the snp list), replace the first two occurrences of '\t' with ';'

with open(gene_annot_file, "r") as rfp:
    with open(gene_annot_file + ".sep", "w") as wfp:
        line = "start"
        while line:
            line = rfp.readline()
            wfp.write(line.replace("\t", ";", 2))
rfp.close()
wfp.close()

gene_annot = pd.read_csv(
    gene_annot_file + ".sep",
    sep=";",
    comment="#",
    names=["Entrez ID", "chr:start:end", "snplist"],
)

display(len(gene_annot))
gene_annot.head()
