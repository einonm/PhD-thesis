# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
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

# # Gene set analysis results summary
#
# ## Drugs & pathways found
#
# This notebook gives summary analysis data showing:
#
# * which ATC code gene sets (i.e. genes identified from the protein targets of the drugs within that set) are most associated with genes implicated in a disease GWAS by a hybrid mapping of those genes to disease SNPs.
# * which GO pathways are associated with the same.
#
# The following SNP-gene annotation methods are analysed:
# * Positional - SNPs are assigned to genes based on their proximity to the gene (within a 10kb window)
# *
# hybrid, hybrid-both, QTL and positional SNP-gene annotation
#
# The analysis is run for each SNP-gene annotation using the following drug target interaction sets:
#
# * DUGGIE   - a set of ATC coded drugs derived from the 8 different drug-gene interaction sources.
# * STITCH   - the largest drug-gene interaction source contributing to DUGGIE
#
# *Uncorrected P-values used throughout, except for Mann-Whitney and comparison of DUGGIE and STITCH analysis performance, where Q-values are used*
#
# The analysis attempts to determine which DTI set and annotation method is most effective at identifying known treatment drugs for a disease (and also not falsely identifying other drugs not used to treat).

# +
import os
import sys
import pandas as pd
import numpy as np

import math

sys.path.append("../lib/")

# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm
import scipy.stats as stats

from IPython.display import display, HTML

import subprocess

import results_display as rdisp

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
# -

HTML(
    """<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>."""
)

# +
# ** Change these values specific to the pipleine runs of interest **
#gwas = "PGC3_SZ"
#run_id = "3qtl"

#gwas = "MDD"
#run_id = "ch5"

#gwas = "BP"
#run_id = "ch5"

#gwas = "AD2022"
#run_id = "32annot"

gwas = "PD2_no23am"
run_id = "ch5"

#gwas = "HD"
#run_id = "32annot"

#gwas = "PD2_no23am"
#run_id = "no_mapt"

dtis = ["duggie", "stitch"]

annots = [
    "prox",
    "alleqtlbrain",
    "alleqtlbrainhybrid",
    "alleqtlbrainhybridboth",
    "allbrain",
    "allbrainhybrid",
    "allbrainhybridboth",
    "alleqtlblood",
    "alleqtlbloodhybrid",
    "alleqtlbloodhybridboth",
    "alleqtltissues",
    "alleqtltissueshybrid",
    "alleqtltissueshybridboth",
    "alltissues",
    "alltissueshybrid",
    "alltissueshybridboth",
]

# +
pd.set_option("display.max_colwidth", -1)
pd.options.display.max_rows = 1000
pd.set_option('display.precision', 3)

data_path = "../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"

N_GENE_THRESH = 4
SIGNIF_THRESH = 0.05
MAGMA_VER = 109
# -

gitspec = subprocess.check_output(
    ["git", "describe", "--all", "--dirty", "--long"]
).strip()
print("This notebook was run on git revision " + str(gitspec, "utf-8"))

# # Analysis of DUGGIE vs STITCH

# +
# All ATC genesets
atc_genesets_df = pd.read_csv(
    data_path + "target-dbs/all_dgi_targets_atc_ensembl.csv",
    header=None,
    sep="\t",
    index_col=0,
)

atc_nostitch_genesets_df = pd.read_csv(
    data_path + "target-dbs/nostitch_dgi_targets_atc_ensembl.csv",
    header=None,
    sep="\t",
    index_col=0,
)

stitch_genesets_df = pd.read_csv(
    data_path + "target-dbs/stitch_dgi_targets_atc_ensembl.csv",
    header=None,
    sep="\t",
    index_col=0,
)

# All GO genesets
go_genesets_df = pd.read_csv(
    "../data/magma_analysis/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt",
    header=None,
    sep="\t",
    index_col=0,
)
# All GO + ATC genesets
genesets_df = pd.concat([atc_genesets_df, go_genesets_df])

nostitch_genesets_df = pd.concat([atc_nostitch_genesets_df, go_genesets_df])

# Ensembl to HGNC translation table
gene_hgnc_df = pd.read_csv(
    "../data/wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    header=None,
    names=["Ensembl", "1", "chrom", "chromStart", "chromEnd", "HGNC"],
)[["Ensembl", "HGNC", "chrom"]]

duggie_drug_dfs = []
stitch_drug_dfs = []
nostitch_drug_dfs = []

dti_df_sets = [duggie_drug_dfs, stitch_drug_dfs, nostitch_drug_dfs]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)


# -


# Set up some plotting variables
fig_width = 18
fig_height = 5
fig_height2 = 34


# +
def generate_plot_list(count):
    plot_list = []
    i=1

    # Create just enough plots for the annotation data sets we have
    while (i <= count):
        plot_list.append([i, i+1])
        i = i + 2

    return plot_list

plot_list = generate_plot_list(len(annots))
# -

# # Comparison of duggie/stitch for all drugs
# Duggie covers more drugs than stitch.

# +
stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        df = dti_df_sets[dtis.index(dti)][annots.index(annot)]
        df, group_labels, colours = rdisp.group_drugs(df, gwas)

    stats_df = stats_df.append(
        rdisp.show_treatment_drug_stats(dti_df_sets[dtis.index(dti)], 3, SIGNIF_THRESH)
    )

display(stats_df)
# -

# # Analysis of annotation methods

# +
annot_dfs = []

for annot in annots:
    annot_dfs.append(rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id))

# +
ATC_LEVELS = rdisp.get_atc_levels(data_path)

for dti in dtis:
    for annot in annots:
        rdisp.summarise_drug_results(
            gwas,
            annot,
            ATC_LEVELS,
            magma_results_dir,
            summary_results_dir,
            dti,
            run_id,
            MAGMA_VER,
            N_GENE_THRESH,
            SIGNIF_THRESH,
        )

for annot in annots:
    rdisp.summarise_gopath_results(
        gwas,
        annot,
        ATC_LEVELS,
        magma_results_dir,
        summary_results_dir,
        run_id,
        MAGMA_VER,
        N_GENE_THRESH,
        SIGNIF_THRESH,
    )
# -

# ## DUGGIE drug results
# ### Significant Q-values
#
# Only up to the first 100 table entries are listed for each drugset, and the first 10 for tables expanded with geneset members. Only the first 50 geneset members are displayed in any table.

for annot in annots:
#for annot in ['alleqtlbrainhybrid']:
    annot_index = annots.index(annot)
    print("")
    print("Significant " + annot + " Q-values, DUGGIE")
    rdisp.display_tables(
        gwas,
        run_id,
        annot,
        "duggie_qvals.tsv",
        "Q " + annot + " duggie",
        annot_dfs[annot_index],
        summary_results_dir,
        genesets_df,
        gene_hgnc_df,
        False,
    )

# ## Results with STITCH only - ATC drugs

# ### Significant Q-values
#
# Only up to the first 100 table entries are listed for each drugset, and the first 10 for tables expanded with geneset members. Only the first 50 geneset members are displayed in any table.

for annot in annots:
    annot_index = annots.index(annot)
    print("")
    print("Significant " + annot + " Q-values, STITCH")
    rdisp.display_tables(
        gwas,
        run_id,
        annot,
        "-stitch_qvals.tsv",
        "Q " + annot + " stitch",
        annot_dfs[annot_index],
        summary_results_dir,
        stitch_genesets_df,
        gene_hgnc_df,
    )

# ## GO path results

# ### Significant Q-values
#
# Only up to the first 100 table entries are listed for each drugset, and the first 10 for tables expanded with geneset members. Only the first 50 geneset members are displayed in any table.

for annot in annots:
    annot_index = annots.index(annot)
    print("")
    print("Significant " + annot + " GO paths, DUGGIE")
    rdisp.display_tables(
        gwas,
        run_id,
        annot,
        "go_qvals.tsv",
        "Q " + annot + " gopaths",
        annot_dfs[annot_index],
        summary_results_dir,
        genesets_df,
        gene_hgnc_df,
    )


