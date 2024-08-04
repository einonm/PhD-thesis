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

# Also generate tables in latex format, with formatting opts
pd.set_option("display.latex.repr", True)
pd.set_option("display.latex.longtable", True)
pd.set_option("display.latex.multirow", True)
pd.set_option("display.latex.multicolumn", True)
pd.set_option("display.max_colwidth", 2000)

# For plots
fig_unit = 17
fig_width = 3
fig_height = 2

# Switch when converting to PDF ('LaTeX') or ODT ('HTML')
# formatting = 'HTML'
formatting = "LaTeX"

if formatting == "HTML":
    print(
        "FORMATTED FOR REVIEW PURPOSES ONLY. SEE CORRESPONDING PDF FOR CORRECT LAYOUT."
    )
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


# # Appendix: Magma permutation issues
#
# ## Large variations in results
#
# During subsequent runs of the processing pipeline as part of the development process it was quickly noticed that gene set analysis results would vary considerably with no significant changes to the input data or pipeline itself, sometimes the difference being between no significant results and more than 10 significant gene sets on repeat runs.
#
# Noting that a default adaptive permutation MAGMA setting is involved when using the multi=snp-wise joint model, and that the MAGMA manual v1.06<citet data-cite="noauthor_magma_nodate"><sup>magmadoc</sup></citet> alludes to a statistical effect of the fixed permutation setting with:
#
# '*...this simple approach is not very efficient however, since a large number of permutations only has added value if the p-value is very low...*'
#
# without being specific as to the meaning of 'low' or 'added value', an exploratory analysis of the effect of the permutation settings on the result was undertaken.
#
# MAGMA allows two *gene-settings* flag modifiers that affect the computation of the multiple testing correction of the empirical p-value which is amalgamated subsequently into a gene's Z-score, and thus affecting the gene set analysis results. These modifiers are *adap-permp* and *fixed-permp*. *adap-permp* is the default taking optional maximum permutations, minimum permutations and stopping criteria parameters which have default values of 1,000,000, 1,000 and 10 respectively; whilst *fixed-permp* takes an optional parameter of the fixed number of permutations to run that has a default of 1,000.

# ## Methods
#
# Initially, 30 runs of an identical pipeline with identical inputs were run using the default adaptive permutation setting on the summary results from a 2014 Parkinson's disease GWAS<citet data-cite="nalls_large-scale_2014"><sup>nalls</sup></citet>. The disease genes were positionally annotated from the GWAS summary SNPs using MAGMA and drug target genes were obtained solely from Drugbank with the gene sets created based on the ATC code of the drug at all 5 levels of the ATC classification system to test with a wide range of gene set sizes.
#
# This scenario was then repeated three more times, identical to the first other than varying the permutation settings by using the *fixed-permp* flag with a value of 5,000, 10,000 and 20,000 fixed permutations for each. The results, consisting of the set of p-values obtained for each run expressed as a negative log of the gene set p-value are shown in figure \ref{fig:pd_sigperms} as box plots for three of the significant gene sets (J06, V10ZA and L04AA23) and three non-significant gene sets (C03, V03AB and A01AB09), chosen to explore the effect of different p-values.
#
# To see if any effect found was an artefact of the Parkinson's GWAS used, three other comparative sets of results were obtained, one for a related Parkinson's meta-analysis<citet data-cite="international_parkinsons_disease_genomics_consortium_meta-analysis_2017"><sup>chang</sup></citet> from 2017, another from a low powered GWAS of Huntington's disease<citet data-cite="genetic_modifiers_of_huntingtons_disease_gem-hd_consortium_electronic_address_gusellahelixmghharvardedu_cag_2019"><sup>hunt</sup></citet> and two more of significant and not significant gene sets, from a well powered Schizophrenia GWAS<citet data-cite="schizophrenia_working_group_of_the_psychiatric_genomics_consortium_biological_2014"><sup>pgc2</sup></citet>.
# The results from these are shown in Figures \ref{fig:pd2_sigperms}, \ref{fig:hd_perms}, \ref{fig:sz_sigperms} and \ref{fig:sz_notsigperms} respectively.
#
# Also, to further query if any discrepancy found would be due to the properties of the gene sets used, the experiment was run once more with the 2014 Parkinson's GWAS but this time instead of using the ATC defined gene sets a set of 500 randomly chosen gene sets of sizes between 10 and 1000 genes were used, where the constituent genes are selected randomly from the ATC set used previously. Three of the gene set results, with a range of set sizes are shown in figure \ref{fig:pd_randperms} for comparison.

# \newpage

# ## Results
# ### Parkinson's disease (2014) significant gene sets

# +
data_path = "../../../../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"

gwas = "PD"
qtls = "brain"
run_id = "-no_kaviar"

## Summary table functions
# Returns ATC data (ATC code & description) in a dataframe, given an ATC table file of the same
def get_atc_file(filename):
    result = pd.DataFrame()
    with open(filename, "r") as file:
        lines = file.readlines()
        for line in lines:
            code = line.split(" ")[0]
            desc = " ".join(line.split(" ")[1:]).strip()
            result = result.append(
                {"atc_code": code, "Description": desc}, ignore_index=True
            )

    return result


# Returns ATC data from all 5 ATC levels in one dataframe
def get_atc_levels():
    atc_level_files = [
        data_path + "atc/level-" + str(i) + "-atc.txt" for i in range(1, 6)
    ]

    # convert each atc_level_file to a 2-column data frame, and cat all of them together
    atc_levels = pd.concat(
        [get_atc_file(atc_level_file) for atc_level_file in atc_level_files]
    )
    atc_levels.set_index("atc_code", inplace=True)
    return atc_levels


# Read all results from a magma geneset analysis for this gwas and annotation type
def read_results_files(gwas, annot_type):
    prefix = (
        magma_results_dir
        + gwas
        + "/results/"
        + "magma_geneset_result-"
        + gwas
        + "-"
        + annot_type
        + run_id
    )
    results_fileset = [prefix + "-atc" + str(i) + ".sets.out" for i in range(1, 6)]

    return [
        pd.read_csv(file, comment="#", delim_whitespace=True)
        for file in results_fileset
        if os.path.exists(file)
    ]


# Summarise results by applying significance and number of gene thresholds, and store in a file
# Generates results for both P and Q values.
def summarise_drug_results(
    gwas, annot_type, atc_levels, n_gene_thresh=5, signif_thresh=0.05
):
    results = read_results_files(gwas, annot_type)

    # only consider classes/drugs with Qval < QVAL_THRESH and NGENES >= N_GENE_THRESH
    significants = [result[result["NGENES"] >= n_gene_thresh] for result in results]

    if not significants:
        print("returning - no significants for " + gwas + " " + annot_type)
        return

    # Calculate the q-values for the local results per GWAS (per ATC code)
    for result in significants:
        result.loc[:, "Q"] = qvalue.estimate(np.array(result["P"]))

    # Put all individual ATC results into one big dataframe
    all_significant = pd.concat(significants)

    q_significant = all_significant[all_significant["Q"] < signif_thresh]
    q_final = pd.merge(
        q_significant, atc_levels, right_index=True, left_on="SET"
    ).sort_values("Q")
    q_final.to_csv(
        os.path.join(
            summary_results_dir,
            "drugs_found-" + gwas + run_id + "-" + annot_type + "_qvals.tsv",
        ),
        sep="\t",
        index=False,
    )

    p_significant = all_significant[all_significant["P"] < signif_thresh]
    p_final = pd.merge(
        p_significant, atc_levels, right_index=True, left_on="SET"
    ).sort_values("P")
    p_final.to_csv(
        os.path.join(
            summary_results_dir,
            "drugs_found-" + gwas + run_id + "-" + annot_type + "_pvals.tsv",
        ),
        sep="\t",
        index=False,
    )


def summarise_gopath_results(gwas, annot_type, n_gene_thresh=5, signif_thresh=0.05):
    prefix = magma_results_dir + gwas + "/results/"
    file = (
        prefix
        + "magma_geneset_result-"
        + gwas
        + "-"
        + annot_type
        + run_id
        + "-gopaths.sets.out"
    )

    if not os.path.exists(file):
        display("File not found: " + file)
        return

    result = pd.read_csv(file, comment="#", delim_whitespace=True)

    # only consider classes/drugs with Qval < QVAL_THRESH and NGENES >= N_GENE_THRESH
    significants = result[result["NGENES"] >= n_gene_thresh]

    if significants.empty:
        return

    significants.rename(columns={"SET": "SHORT_NAME", "FULL_NAME": "SET"}, inplace=True)

    # Calculate the q-values for the local results per GWAS (per ATC code)
    significants.loc[:, "Q"] = qvalue.estimate(np.array(result["P"]))

    q_significant = significants[significants["Q"] < signif_thresh]
    q_significant.to_csv(
        os.path.join(
            summary_results_dir,
            "pathways_found-" + gwas + run_id + "-" + annot_type + "_qvals.tsv",
        ),
        sep="\t",
        index=False,
    )

    p_significant = significants[significants["P"] < signif_thresh]
    p_significant.to_csv(
        os.path.join(
            summary_results_dir,
            "pathways_found-" + gwas + run_id + "-" + annot_type + "_pvals.tsv",
        ),
        sep="\t",
        index=False,
    )


def __main__():
    annot_type_list = ["func-" + qtls, "prox", "hybrid-" + qtls]
    N_GENE_THRESH = 2
    SIGNIF_THRESH = 0.05


# ATC_LEVELS = get_atc_levels()
# for annot_type in annot_type_list:
#     summarise_drug_results(gwas, annot_type, ATC_LEVELS, N_GENE_THRESH, SIGNIF_THRESH)
#     summarise_gopath_results(gwas, annot_type, N_GENE_THRESH, SIGNIF_THRESH)
if __name__ == "__main__":
    __main__()


def get_annot_results(annot, gwas):
    file = os.path.join(
        magma_results_dir,
        gwas,
        "results",
        "magma_gene_result-" + gwas + "-" + annot + run_id + ".genes.out",
    )
    df = None

    if os.path.exists(file):
        df = pd.read_csv(file, delim_whitespace=True)
    else:
        display("File not found: " + file)

    return df


# All ATC genesets
atc_genesets_df = pd.concat(
    [
        pd.read_csv(
            data_path + "drugbank/dbank_gene_set-atc" + str(x + 1) + ".txt",
            header=None,
            sep="\t",
            index_col=0,
        )
        for x in range(5)
    ]
)
# All GO genesets
go_genesets_df = pd.read_csv(
    data_path + "/magma_analysis/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt",
    header=None,
    sep="\t",
    index_col=0,
)
# All GO + ATC genesets
genesets_df = pd.concat([atc_genesets_df, go_genesets_df])
# Emsembl to HGNC translation table
gene_hgnc_df = pd.read_csv(
    data_path + "/wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    header=None,
    names=["Ensembl", "1", "chrom", "chromStart", "chromEnd", "HGNC"],
)[["Ensembl", "HGNC"]]
# Hybrid (func + prox) gene analysis results
# hybriddf = get_annot_results('hybrid-' + qtls, gwas)
# Functional gene analysis results
# funcdf   = get_annot_results('func-' + qtls, gwas)
# Proximity gene analysis results
# proxdf   = get_annot_results('prox', gwas)
# -

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison between different gene analysis permutation settings of the variance in association p-value with Parkinson's disease, for three ATC codes that show a significant association.}"
        )
    )

# +
perms = ["adapt", "5k", "10k", "20k"]
idx = range(0, 30)


def plot_adapt_fixed_boxplot(setname, gwas, plotnum, min_ylim, max_ylim):
    set_df = pd.DataFrame(index=idx, columns=perms)
    for perm in perms:
        set_df[perm] = pd.read_csv(
            "~/source/magma-perms-analysis/"
            + gwas
            + "/"
            + "magma-perms-"
            + gwas
            + "-"
            + perm
            + "-"
            + setname
            + ".txt",
            header=None,
            delim_whitespace=True,
        )
    set_df = set_df.apply(np.log10).abs()

    fig.add_subplot(fig_height, fig_width, plotnum)
    graph = sns.boxplot(data=set_df).set(
        xlabel="Permutation setting",
        ylabel="-log(p-value)",
        title=gwas + " p-value variation for gene set " + setname,
        ylim=(min_ylim, max_ylim),
    )
    return plotnum + 1


fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("J06", gwas, plotnum, 0, 6)
plotnum = plot_adapt_fixed_boxplot("V10XA", gwas, plotnum, 0, 6)
plotnum = plot_adapt_fixed_boxplot("L04AA23", gwas, plotnum, 0, 6)
# -

if formatting == "LaTeX":
    display(Latex(r"\label{fig:pd_sigperms}"))
    display(Latex(r"\end{figure}"))

# The box plots in figure \ref{fig:pd_sigperms} show a marked discrepancy between the adaptive permutation and other fixed permutation results - in the case of the ATC gene set *J06*, the range of values between the adaptive method and the 20k method are exclusive with no overlap. Assuming that a larger number of fixed permutations gives a more accurate result, this implies that the adaptive setting can give incorrect results, as well as having a wider variance.
#
# Looking at the QQ plots for ATC gene set *J06* in figure \ref{fig:pd_qq_sigperms}, it appears that the discrepancy between the adaptive and fixed permutation results is not due to skewed distributions, as they are all approximately normal.

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison of QQ plots of various permutation settings for Parkinson's disease association with ATC gene set J06. \normalfont{Each blue dot represents a p-value association of the gene set with the disease resulting from one gene analysis run. The red line indicates the normal distribution.}}"
        )
    )

# +
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))

sets_df = pd.DataFrame(index=idx, columns=perms)
for perm in perms:
    sets_df[perm] = pd.read_csv(
        "~/source/magma-perms-analysis/PD/magma-perms-PD-" + perm + "-J06.txt",
        header=None,
        delim_whitespace=True,
    )
    sets_df[perm] = sets_df[perm].apply(np.log10).abs()

ax = fig.add_subplot(fig_height, fig_width, 1)
stats.probplot(sets_df["adapt"], dist="norm", plot=pylab)
ax.set_title("PD/J06 Adaptive p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 2)
stats.probplot(sets_df["10k"], dist="norm", plot=pylab)
ax.set_title("PD/J06 10k p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 3)
stats.probplot(sets_df["20k"], dist="norm", plot=pylab)
ax.set_title("PD/J06 20k p-value QQ plot against Normal")
pylab.show()
# -

if formatting == "LaTeX":
    display(Latex(r"\label{fig:pd_qq_sigperms}"))
    display(Latex(r"\end{figure}"))

# ### Parkinson's disease (2014) non-significant gene sets

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison between different gene analysis permutation settings of the variance in association p-value with Parkinson’s disease, for three ATC codes that show no significant association..}"
        )
    )

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("C03", gwas, plotnum, 0, 0.3)
plotnum = plot_adapt_fixed_boxplot("V03AB", gwas, plotnum, 0, 0.3)
plotnum = plot_adapt_fixed_boxplot("A01AB09", gwas, plotnum, 0, 0.3)

if formatting == "LaTeX":
    display(Latex(r"\label{fig:pd_notsigperms}"))
    display(Latex(r"\end{figure}"))

# Both plots, of significant gene sets (figure \ref{fig:pd_sigperms}) and non-significant gene sets (figure \ref{fig:pd_notsigperms}) show a large increase in variance of the adaptive results compared to any fixed result, with this variance being greater for significant gene sets. A comparison between the two groups of plots shows an greater variance and difference in means between adaptive and fixed permutation methods for significant results.

# ### Parkinson's disease (2017) significant gene sets
#
# As Figure \ref{fig:pd2_sigperms} shows, the difference in variance between the results of the adaptive and fixed permutation methods still exists when a different GWAS is analysed. The difference in means is not as prominent, which may be an effect of the 2017 Parkinson's GWAS being of higher power. The results from analysing the low powered Huntington's GWAS and high powered schizophrenia GWAS may indicate if this effect is related to power.

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison between different gene analysis permutation settings of the variance in association p-value with the 2017 Parkinson's disease GWAS, for three ATC codes that show a significant association.}"
        )
    )

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("J06", "PD2", plotnum, 1, 5)
plotnum = plot_adapt_fixed_boxplot("V10XA", "PD2", plotnum, 1, 5)
plotnum = plot_adapt_fixed_boxplot("L04AA23", "PD2", plotnum, 1, 5)

if formatting == "LaTeX":
    display(Latex(r"\label{fig:pd2_sigperms}"))
    display(Latex(r"\end{figure}"))

# ### Huntington's disease gene sets

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison between different gene analysis permutation settings of the variance in association p-value with Huntington's disease, for three random ATC codes.}"
        )
    )

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("M05", "HD", plotnum, 0, 5)
plotnum = plot_adapt_fixed_boxplot("V03AF", "HD", plotnum, 0, 5)
plotnum = plot_adapt_fixed_boxplot("A11AA02", "HD", plotnum, 0, 5)

if formatting == "LaTeX":
    display(Latex(r"\label{fig:hd_perms}"))
    display(Latex(r"\end{figure}"))

# ### Schizophrenia significant gene sets

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison between different gene analysis permutation settings of the variance in association p-value with Schizophrenia, for three ATC codes that show a significant association.}"
        )
    )

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("C08", "SZ", plotnum, 1, 10)
plotnum = plot_adapt_fixed_boxplot("C08DA", "SZ", plotnum, 1, 10)
plotnum = plot_adapt_fixed_boxplot("C08DA51", "SZ", plotnum, 1, 10)

if formatting == "LaTeX":
    display(Latex(r"\label{fig:sz_sigperms}"))
    display(Latex(r"\end{figure}"))

# ### Schizophrenia non-significant gene sets

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison between different gene analysis permutation settings of the variance in association p-value with Schizophrenia, for three ATC codes that show no significant association.}"
        )
    )

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("N05", "SZ", plotnum, 0, 6)
plotnum = plot_adapt_fixed_boxplot("N05CD", "SZ", plotnum, 0, 6)
plotnum = plot_adapt_fixed_boxplot("N05CD06", "SZ", plotnum, 0, 6)

if formatting == "LaTeX":
    display(Latex(r"\label{fig:sz_notsigperms}"))
    display(Latex(r"\end{figure}"))

# There was no substantial difference in means found between the adaptive and fixed permutation results from analysing either Huntington's disease or schizophrenia, indicating that the large difference in means could be specific to the 2014 Parkinson's GWAS data set. However, the effect of power on the variance of the adaptive method result is still observed - indicating that the adaptive method is better with higher powered analyses, but still not as accurate as any of the fixed permutation results.

# ### Parkinson's disease (2014) with randomly allocated gene sets
#
# The final experiment involves gene sets consisting of a random selection of genes analysed against the positionally annotated 2014 Parkinson's disease GWAS summary results. Here, the p-values are generally lower, which may be expected for a randomly selected gene set, but the increased variance of the adaptive results is still seen. Also observed is a more prominent difference in the adaptive mean compared to any fixed permutation mean. This difference may be a unique property of the 2014 Parkinson's GWAS result set not shared by any other GWAS analysed as part of this exercise.

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Comparison between different gene analysis permutation settings of the variance in association p-value with Parkinson's disease, for three random genesets with the highest significance.}"
        )
    )

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("set_3_177", gwas, plotnum, 0., 5.5)
plotnum = plot_adapt_fixed_boxplot("set_92_218", gwas, plotnum, 0, 5.5)
plotnum = plot_adapt_fixed_boxplot("set_445_499", gwas, plotnum, 0, 5.5)

if formatting == "LaTeX":
    display(Latex(r"\label{fig:pd_randperms}"))
    display(Latex(r"\end{figure}"))
    display(Latex(r"\newpage"))

# ## Conclusions
#
# With a limited exploration into the characteristics of MAGMA's gene analysis adaptive permutation method, it can still be concluded that the method is capable of giving unreliable results with a large variance that may verge on being incorrect. This has not been the case to as great a degree with any of the fixed permutation settings, with the downside to using a large number of fixed permutations being that the resources required to run an analysis are greatly increased.
#
# However, even running with 20,000 fixed permutations on the HPC computing resources employed for this project still results in a full pipeline run being completed in less than 24 hours. Given the accuracy of the result is considerably improved over the adaptive and lower fixed permutation settings, it was decided to use this setting for all analyses undertaken for the project.
