# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import os
import sys
import pandas as pd
import numpy as np

# #!pip install scipy
from scipy import stats

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score

# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt

# import qvalue
from IPython.display import display, HTML, Latex, display_latex, display_html

import warnings

warnings.filterwarnings("ignore")

import subprocess

# #!pip install seaborn
import seaborn as sns
import pylab
import scipy.stats as stats
from scipy.stats import chi2, norm

# %matplotlib inline
# #!pip install venn
import venn

import re

import upsetplot as usp

sys.path.append("../../../../lib/")

import results_display as rdisp

# Also generate tables in latex format, with formatting opts
pd.set_option("display.max_colwidth", 2000)
pd.set_option('display.precision', 3)

# Switch when converting to PDF ('LaTeX') or ODT ('HTML')
#formatting = 'HTML'
formatting = "LaTeX"


if formatting == "LaTeX":
    pd.set_option("display.latex.repr", True)
    pd.set_option("display.latex.longtable", True)
    pd.set_option("display.latex.multirow", True)
    pd.set_option("display.latex.multicolumn", True)

if formatting == "HTML":
    print(
        "FORMATTED FOR REVIEW PURPOSES ONLY. SEE CORRESPONDING PDF FOR CORRECT LAYOUT."
    )

# +
data_path = "../../../../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"

MAGMA_VER = 109
N_GENE_THRESH = 4
SIGNIF_THRESH = 0.05
# -

dummy = HTML(
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

if formatting == "LaTeX":
    display(Latex(r"\doublespacing"))
    display(Latex(r"\setlength{\parskip}{5mm}"))
    display(Latex(r"\newpage"))

# +
#gitspec = subprocess.check_output(
#    ["git", "describe", "--all", "--dirty", "--long"]
#).strip()
#print("This notebook was generated from git revision " + str(gitspec, "utf-8"))
# -


# # Appendix: full neuropsychiatric and neurodegenerative disease results
#
# This appendix lists the full set of results obtained from the 32 gene set analysis pipeline executions performed for each of the six diseases studied in chapter \ref{application-of-the-analysis-pipeline-to-neuropsychiatric-and-neurodegenerative-diseases}.
#
# ## Full results - Schizophrenia

gwas = "PGC3_SZ"
run_id = "3qtl"

# +
# All ATC genesets
atc_genesets_df = pd.read_csv(
    data_path + "target-dbs/all_dgi_targets_atc_ensembl.csv",
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
    data_path + "magma_analysis/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt",
    header=None,
    sep="\t",
    index_col=0,
)
# All GO + ATC genesets
genesets_df = pd.concat([atc_genesets_df, go_genesets_df])

# Emsembl to HGNC translation table
gene_hgnc_df = pd.read_csv(
    data_path + "wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    header=None,
    names=["Ensembl", "1", "chrom", "chromStart", "chromEnd", "HGNC"],
)[["Ensembl", "HGNC"]]

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
display_annots = [
    "Proximity",
    "Brain",
    "Brain-miss.",
    "Brain-comb.",
    "Brain-mQTL",
    "Brain-mQTL-miss.",
    "Brain-mQTL-comb.",
    "Blood",
    "Blood-miss.",
    "Blood-comb.",
    "AlleQTLs",
    "AlleQTLs-miss.",
    "AlleQTLs-comb.",
    "AlleQTLs-mQTL",
    "AlleQTLs-mQTL-miss.",
    "AlleQTLs-mQTL-comb.",
]

all_annots = [annot + " DUGGIE" for annot in display_annots] + [annot + " STITCH" for annot in display_annots]

dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)
        
ATC_LEVELS = rdisp.get_atc_levels(data_path)

db_indications_df = pd.DataFrame.from_dict(rdisp.indications, orient='index')

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 2

# +
stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            dti_df_sets[dtis.index(dti)][annots.index(annot)], gwas
        )
        stats_df = stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, 0.05
            )
        )

# Some formatting
stats_df["Sens."] = stats_df["Sens."].map("{:,.3f}".format)
stats_df["Spec."] = stats_df["Spec."].map("{:,.3f}".format)
stats_df["FDR"] = stats_df["FDR"].map("{:,.3f}".format)
stats_df.drop(["DTI DB"], axis=1, inplace=True)
stats_df['Annotation'] = all_annots
cmatrix_stats_df = stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            cmatrix_stats_df.to_latex(
                column_format="p{4.5cm}p{1cm}R{1.5cm}R{1.5cm}p{1.5cm}R{1.5cm}R{1.5cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Confusion matrix results for schizophrenia, from q-value results at a 5\% significance level",
                label="tab:sz_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
deriv_stats_df = stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

#if formatting == "HTML":
#    display(HTML(deriv_stats_df.to_html(index=False)))
#if formatting == "LaTeX":
#    display(Latex(r"\scriptsize"))
#    display(
#        Latex(
#            deriv_stats_df.to_latex(
#                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
#                multirow=True,
#                multicolumn=True,
#                index=False,
#                caption="Derived confusion matrix results for schizophrenia, from q-value results at a 5\% significance level",
#                label="tab:sz_deriv_drugs",
#                position="htbp",
#                escape=False,
#            )
#        )
#    )
#    display(Latex(r"\normalsize"))
# -

# \nopagebreak

# ### Overlap of significant drugs found per annotation for schizophrenia

sz_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, 
                                           gwas, run_id, SIGNIF_THRESH)

# +
usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(sz_drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=100)

usp_figure.set_size_inches(14,21)
# -

# ### Drug enrichment across annotation results for schizophrenia
#
# *Note that Doxepin appears twice, as it has two ATC codes, one sourced via STITCH, the other DUGGIE.

sz_top_drugs_df = rdisp.get_top_drugs_ch5(sz_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.sz_novel_drugs)

if formatting == "HTML":
    display(HTML(sz_top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            sz_top_drugs_df.to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5.5cm}p{0.8cm}R{6cm}",
                caption=r"Drugs ordered by frequency of appearance as significantly "
                        r"associated with schizophrenia across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated and if the association is novel according to the criteria "
                        r"listed in the aims for chapter 5, with drug name and drug indications taken from "
                        r"DrugBank.",
                label="tab:sz_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# \newpage
# ## Full results - Major depressive disorder

gwas = "MDD"
run_id = "ch5"

# +
dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 2

# +
stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            dti_df_sets[dtis.index(dti)][annots.index(annot)], gwas
        )
        stats_df = stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, 0.05
            )
        )

# Some formatting
stats_df["Sens."] = stats_df["Sens."].map("{:,.3f}".format)
stats_df["Spec."] = stats_df["Spec."].map("{:,.3f}".format)
stats_df["FDR"] = stats_df["FDR"].map("{:,.3f}".format)
stats_df.drop(["DTI DB"], axis=1, inplace=True)
stats_df['Annotation'] = all_annots

cmatrix_stats_df = stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            cmatrix_stats_df.to_latex(
                column_format="p{4.5cm}p{1cm}R{1.5cm}R{1.5cm}p{1.5cm}R{1.5cm}R{1.5cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Confusion matrix for major depressive disorder, from q-value results at a 5\% significance level",
                label="tab:mdd_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
deriv_stats_df = stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

#if formatting == "HTML":
#    display(HTML(deriv_stats_df.to_html(index=False)))
#if formatting == "LaTeX":
#    display(Latex(r"\scriptsize"))
#    display(
#        Latex(
#            deriv_stats_df.to_latex(
#                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
#                multirow=True,
#                multicolumn=True,
#                index=False,
#                caption="Derived confusion matrix results for schizophrenia, from q-value results at a 5\% significance level",
#                label="tab:sz_deriv_drugs",
#                position="htbp",
#                escape=False,
#            )
#        )
#    )
#    display(Latex(r"\normalsize"))
# -

# ### Overlap of significant drugs found per annotation for major depressive disorder

mdd_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, 
                                            gwas, run_id, SIGNIF_THRESH)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(mdd_drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=100)

# \newpage
# ### Drug enrichment across annotation results for major depressive disorder

mdd_top_drugs_df = rdisp.get_top_drugs_ch5(mdd_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.mdd_novel_drugs)

if formatting == "HTML":
    display(HTML(mdd_top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            mdd_top_drugs_df.to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5.5cm}p{0.8cm}R{6cm}",                                                                                                                              
                caption=r"Drugs ordered by frequency of appearance as significantly "
                        r"associated with major depressive disorder across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated and if the association is novel according to the criteria "
                        r"listed in the aims for chapter 5, with drug name and drug indications taken from "
                        r"DrugBank.",
                label="tab:mdd_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# \newpage
# ## Full results - Bipolar disorder

gwas = "BP"
run_id = "ch5"

# +
dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 2

# +
stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            dti_df_sets[dtis.index(dti)][annots.index(annot)], gwas
        )
        stats_df = stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, 0.05
            )
        )

# Some formatting
stats_df["Sens."] = stats_df["Sens."].map("{:,.3f}".format)
stats_df["Spec."] = stats_df["Spec."].map("{:,.3f}".format)
stats_df["FDR"] = stats_df["FDR"].map("{:,.3f}".format)
stats_df.drop(["DTI DB"], axis=1, inplace=True)
stats_df['Annotation'] = all_annots

cmatrix_stats_df = stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            cmatrix_stats_df.to_latex(
                column_format="p{4.5cm}p{1cm}R{1.5cm}R{1.5cm}p{1.5cm}R{1.5cm}R{1.5cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Confusion matrix for bipolar disorder, from q-value results at a 5\% significance level",
                label="tab:bp_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
deriv_stats_df = stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

#if formatting == "HTML":
#    display(HTML(deriv_stats_df.to_html(index=False)))
#if formatting == "LaTeX":
#    display(Latex(r"\scriptsize"))
#    display(
#        Latex(
#            deriv_stats_df.to_latex(
#                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
#                multirow=True,
#                multicolumn=True,
#                index=False,
#                caption="Derived confusion matrix results for schizophrenia, from q-value results at a 5\% significance level",
#                label="tab:sz_deriv_drugs",
#                position="htbp",
#                escape=False,
#            )
#        )
#    )
#    display(Latex(r"\normalsize"))
# -

# \nopagebreak
# ### Overlap of significant drugs found per annotation for bipolar disorder

bp_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, 
                                           gwas, run_id, SIGNIF_THRESH)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(bp_drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=100)
usp_figure.set_size_inches(14,21)

# ### Drug enrichment across annotation results for bipolar disorder

bp_top_drugs_df = rdisp.get_top_drugs_ch5(bp_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.bp_novel_drugs)

if formatting == "HTML":
    display(HTML(bp_top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            bp_top_drugs_df.to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5.5cm}p{0.8cm}R{6cm}",                                                                                                                               
                caption=r"Drugs ordered by frequency of appearance as significantly "
                        r"associated with bipolar disorder across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated and if the association is novel according to the criteria "
                        r"listed in the aims for chapter 5, with drug name and drug indications taken from "
                        r"DrugBank.",
                label="tab:bp_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# \newpage
# ## Full results -  Alzheimer's disease

gwas = "AD2022"
run_id = "32annot"

# +
dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 2

# +
stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            dti_df_sets[dtis.index(dti)][annots.index(annot)], gwas
        )
        stats_df = stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, 0.05
            )
        )

# Some formatting
stats_df["Sens."] = stats_df["Sens."].map("{:,.3f}".format)
stats_df["Spec."] = stats_df["Spec."].map("{:,.3f}".format)
stats_df["FDR"] = stats_df["FDR"].map("{:,.3f}".format)
stats_df.drop(["DTI DB"], axis=1, inplace=True)
stats_df['Annotation'] = all_annots

cmatrix_stats_df = stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            cmatrix_stats_df.to_latex(
                column_format="p{4.5cm}p{1cm}R{1.5cm}R{1.5cm}p{1.5cm}R{1.5cm}R{1.5cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Confusion matrix for Alzheimer's disease, from q-value results at a 5\% significance level",
                label="tab:al_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
deriv_stats_df = stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

#if formatting == "HTML":
#    display(HTML(deriv_stats_df.to_html(index=False)))
#if formatting == "LaTeX":
#    display(Latex(r"\scriptsize"))
#    display(
#        Latex(
#            deriv_stats_df.to_latex(
#                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
#                multirow=True,
#                multicolumn=True,
#                index=False,
#                caption="Derived confusion matrix results for schizophrenia, from q-value results at a 5\% significance level",
#                label="tab:sz_deriv_drugs",
#                position="htbp",
#                escape=False,
#            )
#        )
#    )
#    display(Latex(r"\normalsize"))
# -

# ### Overlap of significant drugs found per annotation for Alzheimer's disease

al_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, 
                                           gwas, run_id, SIGNIF_THRESH)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(al_drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=100)

# \newpage
# ### Drug enrichment across annotation results for Alzheimer's disease

al_top_drugs_df = rdisp.get_top_drugs_ch5(al_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.al_novel_drugs)

if formatting == "HTML":
    display(HTML(al_top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            al_top_drugs_df.to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5.5cm}p{0.8cm}R{6cm}",
                caption=r"Drugs ordered by frequency of appearance as significantly "
                        r"associated with Alzheimer's disease across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated and if the association is novel according to the criteria "
                        r"listed in the aims for chapter 5, with drug name and drug indications taken from "
                        r"DrugBank.",
                label="tab:al_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# \newpage
# ## Full results - Parkinson's disease

gwas = "PD2_no23am"
run_id = "ch5"

# +
dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 2

# +
stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            dti_df_sets[dtis.index(dti)][annots.index(annot)], gwas
        )
        stats_df = stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, 0.05
            )
        )

# Some formatting
stats_df["Sens."] = stats_df["Sens."].map("{:,.3f}".format)
stats_df["Spec."] = stats_df["Spec."].map("{:,.3f}".format)
stats_df["FDR"] = stats_df["FDR"].map("{:,.3f}".format)
stats_df.drop(["DTI DB"], axis=1, inplace=True)
stats_df['Annotation'] = all_annots

cmatrix_stats_df = stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            cmatrix_stats_df.to_latex(
                column_format="p{4.5cm}p{1cm}R{1.5cm}R{1.5cm}p{1.5cm}R{1.5cm}R{1.5cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Confusion matrix for Parkinson's disease, from q-value results at a 5\% significance level",
                label="tab:pd_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
deriv_stats_df = stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

#if formatting == "HTML":
#    display(HTML(deriv_stats_df.to_html(index=False)))
#if formatting == "LaTeX":
#    display(Latex(r"\scriptsize"))
#    display(
#        Latex(
#            deriv_stats_df.to_latex(
#                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
#                multirow=True,
#                multicolumn=True,
#                index=False,
#                caption="Derived confusion matrix results for schizophrenia, from q-value results at a 5\% significance level",
#                label="tab:sz_deriv_drugs",
#                position="htbp",
#                escape=False,
#            )
#        )
#    )
#    display(Latex(r"\normalsize"))
# -

# ### Overlap of significant drugs found per annotation for Parkinson's disease

pd_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, 
                                           gwas, run_id, SIGNIF_THRESH)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(pd_drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=100)

# \newpage
# ### Drug enrichment across annotation results for Parkinson's disease

top_drugs_df = rdisp.get_top_drugs_ch5(pd_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.pd_novel_drugs)

if formatting == "HTML":
    display(HTML(top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            top_drugs_df.to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5.5cm}p{0.8cm}R{6cm}",
                caption=r"Drugs ordered by frequency of appearance as significantly "
                        r"associated with Parkinson's disease across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated and if the association is novel according to the criteria "
                        r"listed in the aims for chapter 5, with drug name and drug indications taken from "
                        r"DrugBank.",
                label="tab:pd_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# \newpage
# ## Full results - Huntington's disease

gwas = "HD"
run_id = "32annot"

# +
dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 2

# +
stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            dti_df_sets[dtis.index(dti)][annots.index(annot)], gwas
        )
        stats_df = stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, 0.05
            )
        )

# Some formatting
stats_df["Sens."] = stats_df["Sens."].map("{:,.3f}".format)
stats_df["Spec."] = stats_df["Spec."].map("{:,.3f}".format)
stats_df["FDR"] = stats_df["FDR"].map("{:,.3f}".format)
stats_df.drop(["DTI DB"], axis=1, inplace=True)
stats_df['Annotation'] = all_annots

cmatrix_stats_df = stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            cmatrix_stats_df.to_latex(
                column_format="p{4.5cm}p{1cm}R{1.5cm}R{1.5cm}p{1.5cm}R{1.5cm}R{1.5cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Confusion matrix for Huntington's disease, from q-value results at a 5\% significance level",
                label="tab:hd_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
deriv_stats_df = stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

#if formatting == "HTML":
#    display(HTML(deriv_stats_df.to_html(index=False)))
#if formatting == "LaTeX":
#    display(Latex(r"\scriptsize"))
#    display(
#        Latex(
#            deriv_stats_df.to_latex(
#                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
#                multirow=True,
#                multicolumn=True,
#                index=False,
#                caption="Derived confusion matrix results for schizophrenia, from q-value results at a 5\% significance level",
#                label="tab:sz_deriv_drugs",
#                position="htbp",
#                escape=False,
#            )
#        )
#    )
#    display(Latex(r"\normalsize"))
# -

# ### Overlap of significant drugs found per annotation for Huntington's disease

hd_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, 
                                           gwas, run_id, SIGNIF_THRESH)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(hd_drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=100)

# \newpage
# ### Drug enrichment across annotation results for Huntington's disease

hd_top_drugs_df = rdisp.get_top_drugs_ch5(hd_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.hd_novel_drugs)

if formatting == "HTML":
    display(HTML(hd_top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hd_top_drugs_df.to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5.5cm}p{0.8cm}R{6cm}",
                caption=r"Drugs ordered by frequency of appearance as significantly "
                        r"associated with Huntington's disease across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated and if the association is novel according to the criteria "
                        r"listed in the aims for chapter 5, with drug name and drug indications taken from "
                        r"DrugBank.",
                label="tab:hd_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))
