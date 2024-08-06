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

# #!pip install rpy2

# %matplotlib inline
import matplotlib as mpl
import matplotlib.pyplot as plt

from IPython.display import display, Latex, display_latex, HTML

import warnings

warnings.filterwarnings("ignore")

import subprocess

import seaborn as sns
import pylab
import scipy.stats as stats
from scipy.stats import chi2, norm

import statsmodels.api as sm

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score

import re

import venn

sys.path.append("../../../../lib/")

import results_display as rdisp

# Also generate tables in latex format, with formatting opts
pd.set_option("display.latex.repr", True)
pd.set_option("display.latex.longtable", True)
pd.set_option("display.latex.multirow", True)
pd.set_option("display.latex.multicolumn", True)
pd.set_option("display.max_colwidth", 2000)
pd.options.display.max_rows = 1000
pd.options.display.float_format = "{:.3e}".format
mpl.rcParams.update({"font.size": 16})

# Switch when converting to PDF ('LaTeX') or ODT ('HTML')
#formatting = 'HTML'
formatting = "LaTeX"

if formatting == "HTML":
    print(
        "FORMATTED FOR REVIEW PURPOSES ONLY. SEE CORRESPONDING PDF FOR CORRECT LAYOUT."
    )
# -

if formatting == "LaTeX":
    display(Latex(r"\doublespacing"))
    display(Latex(r"\setlength{\parskip}{5mm}"))
    display(Latex(r"\newpage"))

# +
#gitspec = subprocess.check_output(
#    ["git", "describe", "--all", "--dirty", "--long"]
#).strip()
#print("This chapter was generated from git revision " + str(gitspec, "utf-8"))
# -


# # Appraisal of the pipeline performance with hypertension and hypercholesterolemia

# ## Introduction
#
# Diseases with well understood aetiology and pathogenesis, alongside an available familiarity with the disease genetic architecture and an established set of drug treatment options offer a solid basis for evaluating the analysis pipeline. Hence to explore the outcomes of modifying the pipeline, primary hypertension<citet data-cite="padmanabhan_genomics_2021"><sup>ht_genomics</sup></citet>, ICD-10 (international classification of diseases, 10th revision) code I10, and lipoprotein metabolism disorders and other lipidaemias<citet data-cite="sharifi_genetic_2017"><sup>hc_genomics</sup></citet>, ICD-10 code E78 - which will be referred to here by its most prominent disorder, hypercholesterolemia, were chosen.

# Hypertension is defined as an abnormally high blood pressure, and is an important risk factor for potentially fatal cardiovascular disease such as stroke and myocardial infarction. Primary hypertension is the presence of high blood pressure with no known cause, contributing to up to 90% of cases with a heritability of 30-50%<citet data-cite="ehret_6_2018"><sup>ehret_6_2018</sup></citet>. The genetics of primary hypertension is governed by common variants, as rare variants can be identified as a cause of the disease and therefore not classed as primary hypertension. There are around 90 genes identified as having an involvement with primary hypertension<citet data-cite="ehret_6_2018"><sup>ehret_6_2018</sup></citet>.
#
# Hypercholesterolemia is an asymptomatic complex condition defined by elevated levels of blood cholesterol, particularly non-high density lipoproteins (non-HDL) cholesterol. Hypercholesterolemia is associated with various environmental exposures including diet and smoking, common genetic variants and with genetic diseases such as familial hypercholesterolemia and is associated with a higher risk of atherosclerosis, hardening of the arteries, as well as coronary heart disease. Cholesterol levels have been estimated to have a heritability of 40-60% <citet data-cite="tada_hayato_multiple_2014"><sup>tada_2014</sup></citet>.

# In order to ascertain the effectiveness of DUGGIE, the results from a pipeline execution using DUGGIE as the drug-gene interaction data source were compared to those from an identical pipeline execution but using STITCH as the drug-gene interaction source, the largest previously available and largest single contributing database to DUGGIE. A comparison of these two sets of results should highlight whether the extra drugs and targets DUGGIE provides in addition to those present in STITCH result in a greater ability to associate more drugs with a disease, rejecting non-treatment drugs and selecting treatment drugs amongst those found to be significantly associated. Only results using the proximity annotation method are used to make this comparison. The effect of other annotation methods are explored in chapter \ref{using-functional-qtl-gene-annotation-to-better-identify-genes-associated-with-disease}.
#
# Of interest is the classification of drugs into treatment drug or non-treatment drug sets, particularly as a novel treatment would appear in the analysis as a non-treatment drug that is significantly associated with the disease. The ability of the pipeline to identify drugs used to treat the disease in clinical practice is also a strong affirmation that the analysis is effective, provided the drugs target the genetic causes of the disease and not just the symptoms. For well-understood diseases such as primary hypertension, the pool of treatment drugs available will be considerable as will be the number of non-treatment drugs that have a side effect that affects blood pressure, positively or negatively. It is expected that a gene-based analysis as used in the analysis pipeline discussed here would implicate drugs having such side effects with the disease as readily as it would implicate a treatment drug, so the numbers of non-treatment drugs significantly associated with a disease may be large. However, the non-treatment drugs associated with the disease, as approved drugs, would be expected to have any side effects documented so drug labels and other documentation sources were explored to confirm any known side effects of interest. These metrics are also confounded by a lack of definitive data on what constitutes an approved treatment drug - often qualitative factors such as price and availability of similar alternatives are used as rejection criteria.

# ## Aims of the chapter
#
# This chapter describes and discusses results from the gene set analysis pipeline run with gene-disease associations obtained from proximally annotating GWAS summary statistics for hypertension and hypercholesterolemia. For each disease two sets of results were obtained from two separate executions of the pipeline - once with DUGGIE as the input drug-gene interaction dataset and once with STITCH, the largest single contributing database to DUGGIE. Results and derived metrics are presented and discussed to ascertain the differences between DUGGIE and STITCH in identifying drugs and treatment drugs that are significantly associated with each disease.
#
# The key result presented for each pipeline execution was a count of drugs whose target genes are significantly enriched for genetic association with the disease. Further analysis was made to gauge the ability of the pipeline with each drug-gene interaction dataset to distinguish between the disease treatment and non-treatment groups of drugs, where a drug is classified as a treatment drug if it is known to be currently used in a clinical setting to treat the disease. A confusion matrix was constructed on which a Fisher exact test was performed to measure the relationship between the treatment/non-treatment status and significant/not-significant association analysis result for the drugs. A Mann-Whitney test between the set of drug p-values of treatment and non-treatment drugs gave a measure of the analysis capability to discern between the two groups without imposing a pre-defined significance criterion. Receiver-operator characteristic (ROC) curves were also plotted to provide a visual representation corresponding to the Mann-Whitney test, which is an equivalent to the area under the curve (AUC) measurement of the ROC.
#
# Using the data captured in the confusion matrix the ability of the pipeline with each drug-gene interaction database to identify treatment drugs was considered by calculating the sensitivity, otherwise known as the true positive rate. Similarly, the ability to reject non-treatment drugs was measured using specificity, the true negative rate. Also as novel treatment drugs would be classed as significant non-treatment drugs, the proportion of non-treatment drugs counted as significant was recorded by computing the False Discovery Rate (FDR).
#
# Following this, the effect of the increment of gene targets added by DUGGIE over and above those present in STITCH was examined by looking at the relationship between gene count difference and p-value difference for drugs present in both drug-target interaction databases.
#
# To complete the presentation of evidence used to compare the two drug-gene interaction databases, for each disease the significant drugs identified were listed along with their respective DUGGIE and STITCH pipeline run association q-values, treatment drug indication and whether evidence exists of known drug side effects that affect blood pressure or cholesterol for hypertension and hypercholesterolemia respectively.

# ## Methods and materials

# ### Datasets
#
# #### MAGMA proximity annotation reference
#
# The MAGMA tool provides a simple position-based method and dataset annotating SNPs to genes. This method assigns a SNP to a gene if the SNPs location is inside, or in close proximity to the gene according to an adjustable region around the gene. This region is specified as an extension to the transcript start and stop points of the gene with different possible upstream and downstream values. A symmetrical value of 10 kilobases for this window was used, previously shown to yield lower p-values compared to using no such window<citet data-cite="leeuw_magma_2015"><sup>leeuw</sup></citet>. The corresponding annotation dataset was created using this method and the 1000 Genomes project phase 3 European sub-population dataset, available from the MAGMA website<citet data-cite="noauthor_magma_nodate-1"><sup>magmaweb</sup></citet>. When the terms 'proximally' or 'proximity' mapping and annotation are used in this study, they are referring to this method.
#
# #### Geneatlas GWAS summary statistics
#
# The hypercholesterolemia and hypertension disease GWAS datasets used for this chapter were sourced from GeneATLAS<citet data-cite="canela-xandri_atlas_2018"><sup>geneatlas</sup></citet>, a database of genome-wide associations between hundreds of traits and millions of variants using the UK Biobank cohort. These SNP-based GWAS results were annotated using the proximity annotation reference file provided with MAGMA described above and used as one of the pair of input gene sets in the gene-set analyses.
#
# #### Drug-gene interaction databases
#
# The curated drug-target gene interaction database described in chapter \ref{curation-of-a-drug-target-interaction-database}, DUGGIE, and a similarly formatted comparison dataset to DUGGIE from its largest contributor STITCH provide the second of the input gene sets into the gene-set analyses.
#
# #### Disease treatment drug sets
#
# There are a wide range of drugs employed to lower blood pressure and treat hypertension. Broadly, drugs classified under the ATC top-level code 'C', cardiovascular system drugs, make up the treatment drug set. The DrugBank website<citet data-cite="law_drugbank_2014"><sup>db</sup></citet> was used to create this list using the search mechanism to search for drugs marked as treatment for an indication of 'High blood pressure (hypertension)'. All ATC codes listed by DrugBank for each drug which includes combination and mixture drugs are incorporated. It is expected that other drugs not approved for hypertension lower or raise blood pressure, and these are acknowledged where they are found to be significantly associated with the disease, to give the most accurate comparison of drugs identified by the gene set analysis against those used for therapeutic treatment. Furthermore, drugs in the treatment set were removed where evidence exists indicating that they are not currently used in practice or withdrawn.
#
# As for hypertension treatment drugs, the DrugBank site search facility was used to find drugs listed as treatment for the indication of 'Heterozygous Familial Hypercholesterolaemia' and any ATC codes listed for these drugs marked as a treatment drug for hypercholesterolemia. The majority of these drugs have an ATC code starting with 'C10', lipid modifying agents that work on the cardiovascular system. Most other drugs not in this category have acetylsalicylic acid in combinations or mixtures with other drugs. As before, drugs in the treatment set with existing evidence that they are not used in practice or have been withdrawn were removed from the set.
#

# ### Statistical analysis
#
# #### Drug-disease association
#
# The principal metric used to compare the two analyses differing in only the input drug-gene interaction database was the number of drugs found to be significantly associated with the disease. This was obtained from the pipeline output results of drug-disease associations after applying a Storey False Discovery Rate (FDR) transformation<citet data-cite="storey_direct_2002"><sup>storey</sup></citet> to each drug's disease association p-value to obtain a q-value that is taken to be a significant result if less than 0.05, again an arbitrary but common choice.
#
# When correcting Type I (false positive) error rates in multiple statistical testing, the more conservative Bonferroni Family-Wise Error Rate correction controls the probability of obtaining at least one Type I error, whereas FDR corrections control the proportion of false positives. This is often preferable as allowing for more than one Type I error in this way can give an analysis increased power, otherwise understood as a greater probability of a true positive.
#
# The Storey FDR correction takes a set of p-value results and corrects each for multiple testing by transforming to a q-value, where the q-value represents the expected proportion of false positives in the result set at, or at a more extreme value than, the p-value, mirroring the p-value which represents the probability of obtaining a result at least as extreme as that observed. This way a FDR threshold can be chosen and any q-value with a value below the chosen threshold marked as a significant result. A cursory check of the Storey FDR correction and underlying p-value distribution was made by plotting a p-value density histogram and QQ-plot for each analysis.

# #### Confusion matrices, specificity, sensitivity, FDR and Fisher's exact test
#
# How accurately does the set of drugs identified reflect the drugs used to treat the disease in a clinical setting? An ideal analysis would identify all treatment drugs as significant and all non-treatment drugs as not significant, but this is assuming that all non-treatment drugs do not have unrecognised effects on the disease aetiology and all drugs having a genetic based therapeutic effect are classed as treatment drugs. For other outcomes a table of the classification attempt can be constructed and comparison metrics for the analysis derived - in this case the Fisher exact test p-value, sensitivity, specificity and FDR. 
#
# This form of table is often referred to as a contingency table when used to show a distribution of samples between two groups following a test, which is often the correct and incorrect assignment of samples to null and alternate hypothesis groups. A confusion matrix is a term used for a specialised contingency table that similarly presents the ability of a classification method to discern successfully between groups. For a binary confusion matrix with two groups as used in this study, one dimension heading of the table, either rows or columns, gives the predicted condition split into positive and negative result counts. The alternate dimension heading has tallies of actual positive and negative conditions. The crossover squares of the table are then filled in with the counts of true and false positives and negatives - for example, with true negatives counting the results where a negative prediction matches an actual negative classification. A descriptive table is shown in table \ref{tab:example_confusion}.
#
# Given the application of this example confusion matrix to the results of a method classifying between treatment and non-treatment drugs for a disease, the term 'actual positive' refers to a drug that is used to treat a disease and 'actual negative' to a drug not used to treat the disease. Similarly 'predicted positive'/'predicted negative' refers to the classification method labelling a drug as a treatment drug or non-treatment drug respectively.
#
# In turn:
#
# * 'True Positives' would be the count of classifications where an actual treatment drug is correctly labelled as a treatment drug
# * 'False Positives' would be the count of classifications where an actual non-treatment drug is incorrectly labelled as a treatment drug
# * 'False Negatives' would be the count of classifications where an actual treatment drug is incorrectly labelled as a non-treatment drug
# * And finally, 'True Negatives' would be the count of classifications where an actual non-treatment drug is correctly labelled as a non-treatment drug.

# +
confusion_df = pd.DataFrame(
    [
        ["True Positives (TP)", "False Negatives (FN)"],
        ["False Positives (FP)", "True Negatives (TN)"],
    ],
    columns=["Predicted Positives (PP)", "Predicted Negatives (PN)"],
    index=["Actual Positives (P)", "Actual Negatives (N)"],
)

if formatting == "HTML":
    display(HTML(confusion_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            confusion_df.to_latex(
                column_format="p{5cm}p{5cm}p{5cm}",
                multirow=True,
                multicolumn=True,
                caption=(r"Confusion matrix descriptive layout of actual and predicted positive/negative "
                         r"values, with the crossover squares showing the counts of true and false positives "
                         r"and negatives."),
                
                label="tab:example_confusion",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# Using the terminology from table \ref{tab:example_confusion}, the specificity, the proportion of actual negatives correctly identified can be calculated as:

display(Latex(r"$$"))
display(Latex(r"Specificity = \frac{TN} {TN + FP}"))
display(Latex(r"$$"))

# The sensitivity, the proportion of actual positives correctly identified, can be calculated as:

display(Latex(r"$$"))
display(Latex(r"Sensitivity = \frac{TP} {TP + FN}"))
display(Latex(r"$$"))

# And the False Discovery Rate (FDR), the proportion of non-treatment drugs that are significant can be found using:

display(Latex(r"$$"))
display(Latex(r"FDR = \frac{FP} {FP + TP}"))
display(Latex(r"$$"))

# Fisher's exact test is used to determine whether or not there is a significant relationship between two categorical variables, with the null hypothesis stating that no relationship exists and the variables are independent. The alternate hypothesis is therefore one of a relationship between the variables, and that they are not independent. Because this relationship may be in either direction, a two-tailed Fisher's exact test is used, with a chosen significance level of 0.05.
#
# For each disease, a table was created containing the data from confusion matrices for both DUGGIE and STITCH along with a second table showing metrics derived from the confusion matrices. The listed *significant* and *non-significant* counts are the equivalent of predicted positive and predicted negative test results respectively from table \ref{tab:example_confusion}, similarly the *treatment* and *non treatment* drugs are the equivalent of actual positives and actual negatives. The number of *significant treatment drugs* can therefore be considered true positives and the tally of *significant non-treatment drugs* would be false positives. Similarly the *non-significant treatment drugs* and *non-significant non-treatment drugs* identified by the analysis are respectively the false negatives and true negatives of the confusion matrix. The specificity, sensitivity, FDR and Fisher's exact test p-values are calculated using these definitions. It should be noted that a potential novel therapeutic drug would appear as a significant non-treatment drug, or false positive in the confusion matrix and the FDR gives a metric that is useful for measuring the relative ability of the analyses to find novel candidate therapeutic drugs.
#
# #### Mann-Whitney test and the Receiver-Operator Characteristic
#
# The final column in the derived metrics table lists the Mann-Whitney p-value from a test of all significant and non-significant q-values between the treatment and non-treatment drug groups. This measures the ability of each analysis to discern treatment from non-treatment drugs without imposing a significance threshold, a measurement also explored visually using plotted receiver-operator characteristic (ROC) curves. Also plotted was a boxplot graph of the data sets used to calculate the Mann-Whitney p-values, highlighting the characteristics of each drug target interaction database's analysis in differentiating between treatment and non-treatment drug groups.

# #### P-value correlations
#
# The correlation between the negative log p-values of drugs from the DUGGIE and STITCH pipeline results was investigated in a correlation plot for each disease, which by definition involve only drugs present in both the STITCH and DUGGIE analyses with each marker representing one drug. Because of the large overlap in target genes for each drug, with the STITCH gene targets always a subset of the gene targets recorded in DUGGIE, a high degree of correlation is expected. But also, any difference in p-value for a drug between the analyses of the two drug-gene interaction databases can be attributed to the extra gene targets assigned in DUGGIE.
#
# Two versions of the correlation plot were presented, the first plot coloured to highlight treatment/non-treatment drugs and unique markers used to differentiate between drugs found to be significant/non-significant for DUGGIE and/or STITCH. The second plot uses uniform markers, but differentiates each marker by colour according to which quartile of the drug's gene count difference between the DUGGIE and STITCH datasets.
#
# #### P-value difference and drug gene count difference correlations
#
# The size of effect of each extra target gene assigned to a drug in DUGGIE was explored by first plotting the relationship between the number of extra target genes DUGGIE assigns to a drug over STITCH and the relative difference in drug association negative log p-values between analyses using the two drug target databases.
#
# The negative log p-value difference is plotted on the Y-axis such that a positive value indicates DUGGIE has a more significant (lower) p-value than STITCH, and a negative value indicates STITCH has the more significant p-value. Drugs which are used in the treatment of the disease are coloured differently to non-treatment drugs. Two regression lines are added for the linear regression for all drugs and treatment drugs only.
#
# #### Significantly associated drugs table
#
# All significant drugs found by both DUGGIE and STITCH were listed in a table, ordered by DUGGIE q-value, along with Venn diagrams showing the overlap between the sets of DUGGIE identified significantly associated drugs and STITCH significantly associated drugs. The treatment/non-treatment status of each drug is listed, as well as if the drug has a known side effect of interest, e.g. affecting blood pressure for hypertension.

# ## Hypertension results

# +
data_path = "../../../../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"

ht_gwas = "GA_I10"
qtls = "all"
run_id = "108b_both"
MAGMA_VER = 108

N_GENE_THRESH = 4
SIGNIF_THRESH = 0.05

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
annots = ["prox", "qtl-" + qtls, "hybrid-" + qtls, "hybridboth-" + qtls]

ht_dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            ht_gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        ht_dti_df_sets[dtis.index(dti)].append(newlist)

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 2
# -
# For selecting drugs that may be therapeutic in the treatment of hypertension, the analysis pipeline using proximally annotated genes for both DUGGIE and STITCH produced many significant results - 59 and 42 drugs respectively using a q-value significance threshold of 0.05. Hence DUGGIE is more successful than STITCH in the principal metric of the number of drugs identified as significant.

# ### Hypertension - proximity annotated drug-disease association p-values
#
# Figures \ref{fig:ht_duggie_drug_pvals} and \ref{fig:ht_stitch_drug_pvals} show the p-value density histograms and QQ-plots for DUGGIE and STITCH. A density histogram shows a histogram of p-values along the x-axis with the area of all histogram bars normalised to one and the y-axis labelled accordingly. Both of the histograms figure \ref{fig:ht_duggie_drug_pvals} (a) and figure \ref{fig:ht_stitch_drug_pvals} (a) show a large proportion of true positives after the Storey q-value false discovery rate (FDR) correction, and a broadly asymptotic progression towards higher p-values, meaning that the Storey FDR correction is applicable. The QQ plots for both in (b) of p-value Z-scores, normally distributed for an expected null-hypothesis uniform p-value distribution also show a break from expected values, equivalent to a skewed peak at low p-values as seen in the histograms.

bold_caption = (
    r"Normalised density histogram (a) and QQ plot (b) of association p-values of drug "
    r"target sets with hypertension, for DUGGIE with proximally annotated genes."
)
normal_caption = (
    r"The area under the histogram (a) sums to 1 for comparison and the significant p-values are "
    r"coloured orange. The red dashed line indicates the Storey pi0 value, an estimate of the "
    r"false discovery rate (FDR)."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic(
    [[1, 2]]
)
ax1 = plt.gca()
dti = "duggie"

# plot histogram
plt.subplot(fig_height, fig_width, 1, sharex=ax1, sharey=ax1)
rdisp.plot_density_histogram(
    ht_gwas,
    "prox",
    "drug",
    ht_dti_df_sets[dtis.index(dti)][annots.index("prox")],
    "P prox " + dti,
)
plt.xlabel("DUGGIE Drug-disease association p-value")
plt.title("(a)")

# plot QQ
ax2 = plt.subplot(fig_height, fig_width, 2)
ht_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti] = stats.norm.ppf(
    ht_dti_df_sets[dtis.index(dti)][annots.index("prox")]["P prox " + dti]
)
sm.qqplot(
    ht_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti],
    line="45",
    ax=ax2,
)
plt.ylabel("Measured DUGGIE drug-disease association Z-score")
plt.xlabel("Expected normal quantiles")
dummy = plt.title("(b)")
# -

rdisp.end_caption(formatting, r"fig:ht_duggie_drug_pvals")

bold_caption = (
    r"Normalised density histogram (a) and QQ plot (b) of association p-values of drug "
    r"target sets with hypertension, for STITCH with proximally annotated genes."
)
normal_caption = (
    r"The area under the histogram (a) sums to 1 for comparison and the significant p-values are "
    r"coloured orange. The red dashed line indicates the Storey pi0 value, an estimate of the "
    r"false discovery rate (FDR)."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic(
    [[1, 2]]
)
ax1 = plt.gca()
dti = "stitch"

# plot histogram
plt.subplot(fig_height, fig_width, 1, sharex=ax1, sharey=ax1)
rdisp.plot_density_histogram(
    ht_gwas,
    "prox",
    "drug",
    ht_dti_df_sets[dtis.index(dti)][annots.index("prox")],
    "P prox " + dti,
)
plt.xlabel("STITCH Drug-disease association p-value")
plt.title("(a)")

# plot QQ
ax2 = plt.subplot(fig_height, fig_width, 2)
ht_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti] = stats.norm.ppf(
    ht_dti_df_sets[dtis.index(dti)][annots.index("prox")]["P prox " + dti]
)
sm.qqplot(
    ht_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti],
    line="45",
    ax=ax2,
)
plt.ylabel("Measured STITCH drug-disease association Z-score")
plt.xlabel("Expected normal quantiles")
dummy = plt.title("(b)")

# TODO - get axes of QQ plot to match
# -

rdisp.end_caption(formatting, r"fig:ht_stitch_drug_pvals")

# \newpage
# ### Hypertension - confusion matrices and derived metrics
#
# Table \ref{tab:ht_c3_cmatrix_drugs} contains results data extracted from confusion matrices for both DUGGIE and STITCH when analysing hypertension, shown together for ease of comparison. Table \ref{tab:ht_c3_deriv_drugs} displays the metrics derived from these confusion matrices for each of the two drug-gene interaction databases, of the sensitivity, specificity, false discovery rate and the Mann-Whitney and Fisher p-values.

# +
ht_stats_df = pd.DataFrame()
for dti in dtis:
    (
        ht_dti_df_sets[dtis.index(dti)][annots.index("prox")],
        group_labels,
        colours,
    ) = rdisp.group_drugs(
        ht_dti_df_sets[dtis.index(dti)][annots.index("prox")], ht_gwas
    )
    ht_stats_df = ht_stats_df.append(
        rdisp.show_treatment_drug_stats(
            list([ht_dti_df_sets[dtis.index(dti)][annots.index("prox")]]), 3, 0.05
        )
    )

ht_stats_df.drop(["Annotation"], axis=1, inplace=True)

# Some formatting
ht_stats_df["Sens."] = ht_stats_df["Sens."].map("{:,.3f}".format)
ht_stats_df["Spec."] = ht_stats_df["Spec."].map("{:,.3f}".format)
ht_stats_df["FDR"] = ht_stats_df["FDR"].map("{:,.3f}".format)
ht_stats_df["DTI DB"] = [x.upper() for x in ht_stats_df["DTI DB"]]

ht_cmatrix_stats_df = ht_stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(ht_cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_cmatrix_stats_df.to_latex(
                column_format="p{2cm}p{1.2cm}p{2cm}R{1.2cm}p{1.2cm}p{1.4cm}R{3cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=(r"Confusion matrix results for hypertension, obtained from a q-value cutoff at a "
                        r"0.05 significance level. \normalfont{Covering significant and non-significant "
                        r"results for treatment and non-treatment groups of drugs.}"),
                label="tab:ht_c3_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# +
ht_deriv_stats_df = ht_stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

if formatting == "HTML":
    display(HTML(ht_deriv_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_deriv_stats_df.to_latex(
                column_format="p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}p{3.5cm}p{1cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=(r"Derived confusion matrix results for hypertension, obtained from a q-value cutoff "
                         r"at a 0.05 significance level. \normalfont{Covering significant and non-significant "
                          r"results for treatment and non-treatment groups of drugs.}"),
                label="tab:ht_c3_deriv_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# Here temporarily to prevent text going off end of PDF page after 2 tables
if formatting == "LaTeX":
    display(Latex(r"\newpage"))

# Both analyses are picking up comparable numbers of treatment drugs, with STITCH identifying a greater proportion of treatment drugs giving a higher sensitivity, attributable to a smaller pool of treatment drugs involved in the analysis - indicating that the gene targets present in STITCH may be of a higher quality than those in DUGGIE as a whole. The proportion of non-treatment drugs being rejected is also very similar as seen in the specificity values, but DUGGIE is showing an advantage over STITCH with twice the number of non-treatment drugs available to be analysed. The FDR from the DUGGIE analysis is marginally higher than that for STITCH, as can be seen in the higher proportion and number of false positives / significant non-treatment drugs found. These may be genuine false positives but could also be novel undiscovered treatment drugs.
#
# Overall, the Fisher p-values indicate that both analyses have a useful capacity to detect treatment drugs, with STITCH having a better ability but selecting treatment drugs from a smaller pool. The increased Mann-Whitney p-value of DUGGIE is considered to be due to the fact that DUGGIE has almost twice the number of non-treatment drugs but is rejecting them correctly at a similar ratio to STITCH.

# #### Hypertension - ROC curves
#
# Plotting the ROC curves for DUGGIE and STITCH in figure \ref{fig:hc_prox_roc} reveals that both are broadly identical at identifying treatment drugs correctly, with STITCH gaining a small advantage of higher positive identifications at the more desirable lower false positive rates. The effect of this can be seen at the black dashed markers representing the significant treatment and significant non-treatment rates at a 0.05 significance level for both DUGGIE and STITCH, where STITCH has a higher significant treatment rate than DUGGIE, and subsequently a smaller Fisher test p-value. However, as the DUGGIE analysis involves more drugs overall, it does indicate that DUGGIE is capable of identifying more treatment drugs without a significant loss of accuracy.

bold_caption = (
    r"Receiver-Operator Characteristic (ROC) curves from gene set analyses of hypertension "
    r"using the DUGGIE and STITCH drug-gene interaction databases and proximity gene annotation."
)
normal_caption = (
    r"The rate of correct classification of drugs (significant treatment and non-significant non-treatment) "
    r"is compared against incorrect classification (significant non-treatment and non-significant treatment). "
    r"Area Under the Curve (AUC) values are shown against each annotation method in the key. "
    r"Dashed black lines show the significant treatment and significant non-treatment values for "
    r"both DUGGIE and STITCH at the 0.05 significance level."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
roc_fig_unit = 10
plt.figure(
    figsize=(roc_fig_unit, roc_fig_unit * fig_height / fig_width)
).subplot_mosaic([[1]])
lw = 2
for dti in ["duggie", "stitch"]:
    for annot in ["prox"]:
        df = ht_dti_df_sets[dtis.index(dti)][annots.index(annot)]
        if df["colour_id"].any():
            fpr, tpr, thresholds = roc_curve(
                (df["colour_id"] != 0), 1 - df.iloc[:, [3]]
            )
            roc_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, lw=lw, label=dti.upper() + " (AUC = %0.3f)" % roc_auc)

plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.plot([0.026, 0.026], [0, 0.366], color="black", lw=lw, linestyle="dashed")
plt.plot([0, 0.026], [0.366, 0.366], color="black", lw=lw, linestyle="dashed")
plt.plot([0.021, 0.021], [0, 0.226], color="black", lw=lw, linestyle="dashed")
plt.plot([0, 0.021], [0.226, 0.226], color="black", lw=lw, linestyle="dashed")
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.05])
plt.xlabel("Significant Non-treatment & Non-significant Treatment Rate")
plt.ylabel("Significant Treatment & Non-significant Non-treatment Rate")
plt.legend()
plt.show()
# -

rdisp.end_caption(formatting, r"fig:ht_prox_roc")

# Figure \ref{fig:ht_boxplot_treat_nontreat} shows a boxplot graph of the data set used to calculate the Mann-Whitney p-values, highlighting the ability of each drug target interaction database's analysis to differentiate between the treatment and non-treatment drug groups. It can be seen that although the STITCH difference in means between the treatment and non-treatment drug sets is greater than that for DUGGIE, DUGGIE's improved Mann-Whitney p-value also stems from the lower variation in association values obtained for the treatment drug set.

bold_caption = r"Box plot of hypertension association negative log p-values for treatment and non-treatment drug sets."
normal_caption = (
    r"A comparison between DUGGIE and STITCH using proximally annotated "
    r"genes associated with hypertension."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
fig, ax = plt.subplots(figsize=(10, 6))
duggie_df = ht_dti_df_sets[dtis.index("duggie")][annots.index("prox")]
stitch_df = ht_dti_df_sets[dtis.index("stitch")][annots.index("prox")]

duggie_df["-log(p)"] = -np.log10(duggie_df["P prox duggie"])
stitch_df["-log(p)"] = -np.log10(stitch_df["P prox stitch"])

plot_df = pd.DataFrame(
    {
        "DUGGIE treatment": duggie_df[duggie_df["colour_id"] != 0].loc[:, "-log(p)"],
        "DUGGIE nontreatment": duggie_df[duggie_df["colour_id"] == 0].loc[:, "-log(p)"],
        "STITCH treatment": stitch_df[stitch_df["colour_id"] != 0].loc[:, "-log(p)"],
        "STITCH nontreatment": stitch_df[stitch_df["colour_id"] == 0].loc[:, "-log(p)"],
    }
)

g = sns.boxplot(data=plot_df).set(xlabel="Drug group", ylabel="-log10(p-value)")
dummy = ax.set_xticklabels(
    labels=[
        "DUGGIE treatment",
        "DUGGIE non-treatment",
        "STITCH treatment",
        "STITCH non-treatment",
    ],
    rotation=30,
)
# -

rdisp.end_caption(formatting, r"fig:ht_boxplot_treat_nontreat")

# ### Hypertension - correlation of STITCH and DUGGIE proximally annotated drug p-values
#
# The correlation between the negative logs of the p-values for DUGGIE and STITCH was investigated in Figure \ref{fig:ht_correlation_qvals}. With only drugs that are present in both STITCH and DUGGIE plotted, a high degree of correlation can be seen with an r coefficient of 0.81 as would be expected from the large overlap in target genes for each drug. Treatment drugs are also prevalent at the higher plotted values, encouragingly showing that both databases include treatment drugs as a large proportion of the highly significant drugs. The criterion for significance is similar for both DUGGIE and STITCH, at a -log(p) value of ~2.1, despite DUGGIE testing more drugs so the general significance increase in the STITCH results for treatment drugs may indicate that the extra targets provided by DUGGIE are reducing the power of the DUGGIE analysis, for those drugs that are present in both databases.

bold_caption = (
    r"Hypertension correlation plot of DUGGIE vs STITCH -log(p-value) for proximally "
    r"annotated drugs present in both databases."
)
normal_caption = (
    r"Blue markers are drugs listed in DrugBank as a treatment for the indication of essential "
    r"hypertension. Other drugs are shown as pink markers. Drugs with significant "
    r"p-values for STITCH only, DUGGIE only, both DUGGIE and STITCH and neither are shown as "
    r"down triangles, up triangles, diamonds and points respectively."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic([[1]])
axlimit = 7

dti1 = "stitch"
dti2 = "duggie"

plt.subplot(fig_height, fig_width, 1)
rdisp.plot_drug_scatter(
    ht_gwas,
    "prox " + dti1,
    ht_dti_df_sets[dtis.index(dti1)][annots.index("prox")],
    "prox " + dti2,
    ht_dti_df_sets[dtis.index(dti2)][annots.index("prox")],
    "P",
    axlimit,
)
plt.xlabel("Hypertension STITCH drug association -log10(p-value)")
plt.ylabel("Hypertension DUGGIE drug association -log10(p-value)")
dummy = plt.title("")
# -

rdisp.end_caption(formatting, r"fig:ht_correlation_qvals")

# Exploring the effect of the extra DUGGIE target genes in the shared drug subset, the target gene difference was investigated by looking at the relationship between the number of extra target genes DUGGIE assigns to a drug over STITCH and the difference in drug association absolute p-value between analyses using the two drug target databases, plotted in figure \ref{fig:ht_correlation_gene_counts}.

bold_caption = (
    r"Hypertension correlation plot of DUGGIE vs STITCH -log10(p-value) for proximally "
    r"annotated drugs present in both databases."
)
normal_caption = (
    r"Dots are coloured differently based "
    r"on the size of the difference between drug gene target counts for DUGGIE and STITCH."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic([[1]])
axlimit = 7

dti1 = "stitch"
dti2 = "duggie"

diffs = []

for annot in ["prox"]:
    plt.subplot(fig_height, fig_width, annots.index(annot) + 1)
    df = rdisp.plot_drug_genediff_scatter(
        ht_gwas,
        annot + " " + dti1,
        ht_dti_df_sets[dtis.index(dti1)][annots.index(annot)],
        annot + " " + dti2,
        ht_dti_df_sets[dtis.index(dti2)][annots.index(annot)],
        "P",
        axlimit,
    )
    diffs.append(df)
    plt.xlabel("Hypertension STITCH drug association -log10(p-value)")
    plt.ylabel("Hypertension DUGGIE drug association -log10(p-value)")
    plt.title("")

# -

rdisp.end_caption(formatting, r"fig:ht_correlation_gene_counts")

# This correlation plot highlights a broad trend where the most significant drugs with the highest gene count difference between their DUGGIE and STITCH target gene sets are furthest away from the central equivalence line. There also appears to be a relationship between the groups of colour where a smaller difference in gene count places the drug closer to the equivalence line, so in effect the more genes added by DUGGIE amplifies the p-value result, positively or negatively - but most often negatively. Again, as with the sensitivity metric comparison, this indicates that the extra DUGGIE drug targets are of a lower quality than those already in STITCH.

# ### Hypertension - correlation of p-value difference and drug gene count difference between DUGGIE and STITCH
#
# To investigate further the effect of the extra drug targets in DUGGIE, a correlation plot of the p-value difference between DUGGIE and STITCH and drug target count difference between DUGGIE and STITCH for each drug present in both databases was plotted in figure \ref{fig:ht_prox_diffs}. Drugs which are used in the treatment of hypertension are coloured blue with all other drugs coloured pink. The green and blue regression lines are the linear regressions for all drugs and treatment drugs respectively.

bold_caption = (
    r"Hypertension correlation plot of DUGGIE vs STITCH -log10(p-value) difference vs gene "
    r"count difference for proximally annotated drugs present in both databases."
)
normal_caption = (
    r"Blue dots are anti-hypertensive treatment drugs. The green regression line is for all drugs, "
    r"whilst the blue regression line is for anti-hypertensive treatment drugs only."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
# Diff on x-axis, p-val on y
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic([[1]])

for i in range(0, 1):
    df = diffs[i]

    # -ve diff => STITCH better, +ve DUGGIE better
    df["p-diffs"] = df["prox duggie -log(P)"] - df["prox stitch -log(P)"]

    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df["diff"], df["p-diffs"]
    )

    df, group_labels, colours = rdisp.group_drugs(df, ht_gwas)

    groups = df.groupby("colour_id")

    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(
        df[df["colour_id"] == 1]["diff"], df[df["colour_id"] == 1]["p-diffs"]
    )

    plt.subplot(fig_height, fig_width, i + 1)

    for name, group in groups:
        plt.scatter(
            group["diff"],
            group["p-diffs"],
            label=group_labels[int(name)],
            c=colours[int(name)],
        )

    plt.xlabel("DUGGIE/STITCH gene target count difference")
    plt.ylabel("-log(p-value) DUGGIE/STITCH difference")

    ax2 = plt.gca()
    plt.text(
        0,
        -0.2,
        "r_value = "
        + "{0:.5f}".format(r_value)
        + ", slope = "
        + "{0:.5f}".format(slope)
        + ", err "
        + "{0:.5f}".format(std_err),
        transform=ax2.transAxes,
    )

    plt.text(
        0,
        -0.3,
        "r_value(treatment) = "
        + "{0:.5f}".format(r_value2)
        + ", slope = "
        + "{0:.5f}".format(slope2)
        + ", err "
        + "{0:.5f}".format(std_err2),
        transform=ax2.transAxes,
    )

    legend = ax2.legend([])
    plt.legend()
    plt.plot(df["diff"], slope2 * df["diff"] + intercept2, i)
    dummy = plt.plot(df["diff"], slope * df["diff"] + intercept, i)
# -

rdisp.end_caption(formatting, r"fig:ht_prox_diffs")

# Overall, the all-drug regression line slope indicates that the number of extra drug targets DUGGIE contributes is not significantly related to the p-value difference. However, for some individual drugs there is a marked difference and this is particularly prevalent for treatment drugs which display the highest p-value discrepancies. Taking only the treatment drugs, the regression confirms that there is a heavy downward trend indicating that in general the extra targets contributed by DUGGIE over those in STITCH have a negative effect on the association of treatment drugs with the disease. This could be due to the fact that approved drugs present in STITCH would already contain the known significant therapeutic targets and hence any targets added by DUGGIE are likely to only reduce the significance of the drug's gene target set association with the disease.

# ### Significant hypertension associated  drugs
#
# All significant drugs found by both DUGGIE and STITCH are listed in table \ref{tab:ht_sig_drugs}, ordered by DUGGIE q-value. DUGGIE contains the three drugs with the most significant associations from the two analyses, which are not present in the STITCH analysis as STITCH assigns less than the applied threshold of 4 target genes. At the bottom of the table where drugs that are significant in the STITCH analysis but not in the DUGGIE analysis reside, there is a larger difference in the number of target genes assigned to each drug which would be reducing the association p-values. Venn diagrams showing the overlap between the sets of DUGGIE identified significantly associated drugs and STITCH significantly associated drugs are in figure \ref{fig:ht_all_venn}.
#
# The hypertension treatment/non-treatment status of each drug is listed, as well as if blood pressure changes are a known effect of the drug either as part of DrugBank's 'agents that produce hypertension' list or from a literature search. Only 8 drugs of the 78 that are listed as significant by either the DUGGIE or STITCH analyses do not have a well known relationship to blood pressure. It would be interesting to produce summary statistics including this broader definition, but as this would involve manually and  consistently recording the effect of each drug on blood pressure, was not undertaken due to time constraints.

# \newpage

# +
proxdf = rdisp.get_annot_results(ht_gwas, "prox", magma_results_dir, run_id)
ATC_LEVELS = rdisp.get_atc_levels(data_path)

# We want all results, so we set gene count and q-value theshold to 1.
for dti in dtis:
    for annot in ["prox"]:
        rdisp.summarise_drug_results(
            ht_gwas,
            annot,
            ATC_LEVELS,
            magma_results_dir,
            summary_results_dir,
            dti,
            run_id,
            MAGMA_VER,
            1,
            1,
        )
# -

rdisp.display_joint_table_thesis(
    formatting,
    ht_gwas,
    run_id,
    "prox",
    "duggie_qvals.tsv",
    "stitch_qvals.tsv",
    "prox duggie",
    "prox stitch",
    proxdf,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
    r"Significant hypertension drugs found by DUGGIE and STITCH. ",
    r"q-values in brackets are non-significant and included for completeness. "
    r"Drugs with no STITCH values are only present in the DUGGIE analysis. "
    r"Columns are also included indicating if the drug is a known treatment drug "
    r"and if the drug is known to affect blood pressure.",
    "tab:ht_sig_drugs",
)

bold_caption = (
    r"Venn diagrams of overlap of (a) all significant drugs and (b) significant treatment drugs for "
    r"hypertension."
)
normal_caption = r""
rdisp.start_caption(formatting, bold_caption, normal_caption)

rdisp.display_dgi_venn_thesis(
    ht_gwas,
    run_id,
    "prox",
    "duggie_qvals.tsv",
    "stitch_qvals.tsv",
    "prox duggie",
    "prox stitch",
    proxdf,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
)

rdisp.end_caption(formatting, r"fig:ht_all_venn")

# ## Hypercholesterolemia results
#
# Mirroring the hypertension result when a q-value significance threshold of 0.05 is applied, DUGGIE has a marginal increase over STITCH in the number of significant drugs found - 19 compared to 18 from the STITCH analysis. Given the larger number of both treatment and non-treatment drugs present in the DUGGIE analysis this marginal increase highlights once more that the extra drug target genes contributed by DUGGIE in addition to those by STITCH may not be as effective as those already present in STITCH in differentiating between treatment and non-treatment drugs.

hc_gwas = "GA_E78"

# +
hc_dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            hc_gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        hc_dti_df_sets[dtis.index(dti)].append(newlist)
# -

# ### Hypercholesterolemia - proximity annotated drug-disease association p-values
#
# The drug p-value density histogram and QQ-plot for DUGGIE is shown in figure \ref{fig:hc_duggie_drug_pvals} and the corresponding plots for STITCH are in figure \ref{fig:hc_stitch_drug_pvals}. Both histograms in figure \ref{fig:hc_duggie_drug_pvals} (a) and \ref{fig:hc_stitch_drug_pvals} (a) show a large peak of true positives after FDR correction, but at a somewhat lower level than for hypertension which is reflected in the lower number of significant drugs.

bold_caption = (
    r"Normalised density histogram (a) and QQ plot (b) of association p-values of drug "
    r"target sets with hypercholesterolemia, for DUGGIE with proximally annotated genes."
)
normal_caption = (
    r"The area under the histogram (a) sums to 1 for comparison and the significant "
    r"p-values are coloured orange. The red dashed line indicates the Storey pi0 value, an estimate "
    r"of the false discovery rate (FDR)."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic(
    [[1, 2]]
)
ax1 = plt.gca()
dti = "duggie"

# plot histogram
plt.subplot(fig_height, fig_width, 1, sharex=ax1, sharey=ax1)
rdisp.plot_density_histogram(
    hc_gwas,
    "prox",
    "drug",
    hc_dti_df_sets[dtis.index(dti)][annots.index("prox")],
    "P prox " + dti,
)
plt.xlabel("DUGGIE Drug-disease association p-value")
plt.title("(a)")

# plot QQ
ax2 = plt.subplot(fig_height, fig_width, 2)
hc_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti] = stats.norm.ppf(
    hc_dti_df_sets[dtis.index(dti)][annots.index("prox")]["P prox " + dti]
)
sm.qqplot(
    hc_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti],
    line="45",
    ax=ax2,
)
plt.ylabel("Measured DUGGIE drug-disease association Z-score")
plt.xlabel("Expected normal quantiles")
dummy = plt.title("(b)")
# -

rdisp.end_caption(formatting, r"fig:hc_duggie_drug_pvals")

bold_caption = (
    r"Normalised density histogram (a) and QQ plot (b) of association p-values of drug "
    r"target sets with hypercholesterolemia, for STITCH with proximally annotated genes."
)
normal_caption = (
    r"The area under the histogram (a) sums to 1 for comparison and the significant "
    r"p-values are coloured orange. The red dashed line indicates the Storey pi0 value, an estimate "
    r"of the false discovery rate (FDR)."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic(
    [[1, 2]]
)
ax1 = plt.gca()
dti = "stitch"

# plot histogram
plt.subplot(fig_height, fig_width, 1, sharex=ax1, sharey=ax1)
rdisp.plot_density_histogram(
    hc_gwas,
    "prox",
    "drug",
    hc_dti_df_sets[dtis.index(dti)][annots.index("prox")],
    "P prox " + dti,
)
plt.xlabel("STITCH Drug-disease association p-value")
plt.title("(a)")

# plot QQ
ax2 = plt.subplot(fig_height, fig_width, 2)
hc_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti] = stats.norm.ppf(
    hc_dti_df_sets[dtis.index(dti)][annots.index("prox")]["P prox " + dti]
)
sm.qqplot(
    hc_dti_df_sets[dtis.index(dti)][annots.index("prox")]["Z prox " + dti],
    line="45",
    ax=ax2,
)
plt.ylabel("Measured STITCH drug-disease association Z-score")
plt.xlabel("Expected normal quantiles")
dummy = plt.title("(b)")

# TODO - get axes of QQ plot to match
# -

rdisp.end_caption(formatting, r"fig:hc_stitch_drug_pvals")

# ### Hypercholesterolemia - confusion matrices and derived metrics
#
# Table \ref{tab:hc_cmatrix_drugs_c3} contains results data from confusion matrices for both DUGGIE and STITCH when
# analysing hypercholesterolemia, shown together for ease of comparison. Table \ref{tab:hc_deriv_drugs_c3} lists the metrics derived from these confusion matrices or each of the two drug-gene interaction databases, of the sensitivity, pecificity, false discovery rate and the Mann-Whitney and Fisher p-values.

# +
hc_stats_df = pd.DataFrame()
for dti in dtis:
    (
        hc_dti_df_sets[dtis.index(dti)][annots.index("prox")],
        group_labels,
        colours,
    ) = rdisp.group_drugs(
        hc_dti_df_sets[dtis.index(dti)][annots.index("prox")], hc_gwas
    )
    hc_stats_df = hc_stats_df.append(
        rdisp.show_treatment_drug_stats(
            list([hc_dti_df_sets[dtis.index(dti)][annots.index("prox")]]), 3, 0.05
        )
    )

hc_stats_df.drop(["Annotation"], axis=1, inplace=True)

# Some formatting
hc_stats_df["Fisher p-value"] = (
    hc_stats_df["Fisher p-value"].astype(float).map("{:,.3e}".format)
)
hc_stats_df["Sens."] = hc_stats_df["Sens."].map("{:,.3f}".format)
hc_stats_df["Spec."] = hc_stats_df["Spec."].map("{:,.3f}".format)
hc_stats_df["FDR"] = hc_stats_df["FDR"].map("{:,.3f}".format)
hc_stats_df["DTI DB"] = [x.upper() for x in hc_stats_df["DTI DB"]]

hc_cmatrix_stats_df = hc_stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(hc_cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_cmatrix_stats_df.to_latex(
                column_format="p{2cm}p{1.2cm}p{2cm}R{1.2cm}p{1.2cm}p{1.4cm}R{3cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Confusion matrix for hypercholesterolemia, from q-value results at a 0.05 significance level",
                label="tab:hc_cmatrix_drugs_c3",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
hc_deriv_stats_df = hc_stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

if formatting == "HTML":
    display(HTML(hc_deriv_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_deriv_stats_df.to_latex(
                column_format="p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}p{3cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Derived confusion matrix results for hypercholesterolemia, from q-value results at a 0.05 significance level",
                label="tab:hc_deriv_drugs_c3",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# Here temporarily to prevent text going off end of PDF page after 2 tables
if formatting == "LaTeX":
    display(Latex(r"\newpage"))

# As with the hypertension results, the DUGGIE analysis increase over the STITCH analysis in significant drugs identified does not translate into an increase of significant treatment drugs found. The STITCH analysis again has a higher sensitivity as it is more successful in picking out treatment drugs. However, the DUGGIE gene-set analysis is  successful at marking most of the greater number of non-treatment drugs it covers as non-significant, but as the number of non-significant treatment drugs is also higher, the DUGGIE analysis FDR is slightly higher that the FDR for the STITCH analysis. Further work would be needed to discern if these false positives are genuine or potential novel treatment drugs.
#
# The Fisher p-values, along with the sensitivity and specificity are affected by the very small number of true positives found for DUGGIE and STITCH, 5 and 7 respectively. Although this gives DUGGIE a lower significant Fisher p-value, the ROC curve in figure \ref{fig:hc_prox_roc} shows that over all significance thresholds, the metrics are closer with an AUC for STITCH of 0.93 and 0.88 for DUGGIE. The large improvement seen in the Mann-Whitney p-value for DUGGIE over STITCH is reflected in the box plot of figure \ref{fig:hc_boxplot_treat_nontreat}, where the peak for DUGGIE treatment drug q-values shows a narrower width covering a larger number of drugs compared to STITCH treatment drugs.

# #### Hypercholesterolemia - ROC curves
#
# The hypercholesterolemia ROC curves in figure \ref{fig:hc_prox_roc} shows both the DUGGIE and STITCH analyses classifying treatment and non-treatment drugs with a high accuracy, but with STITCH having an advantage in identifying treatment drugs correctly over the entire range of true and false positive rates. However, the difference is marginal and the DUGGIE analysis involves over twice the number of non-treatment drugs and a greater count of treatment drugs, making it the analysis with the most potential for finding novel treatment drugs. The reason for the smaller STITCH Fisher p-value in table \ref{tab:hc_deriv_drugs_c3} is also highlighted by the black dotted lines, showing the significant treatment rates and significant non-treatment rates for DUGGIE and STITCH at the 0.05 significance level. Here, STITCH is better at identifying significant treatment drugs at very small significant non-treatment rates.

bold_caption = (
    r"Receiver-Operator Characteristic (ROC) curves from gene set analyses of hypercholesterolemia "
    r"using the DUGGIE and STITCH drug-gene interaction databases and proximity gene annotation."
)
normal_caption = (
    r"The rate of correct classification of drugs (i.e. significant treatment and non-significant non-treatment) "
    r"is compared against incorrect classification (i.e. significant non-treatment and non-significant treatment). "
    r"Area Under the Curve (AUC) values are shown against each annotation method in the key. "
    r"Dashed black lines show the significant treatment and significant non-treatment values for "
    r"both DUGGIE and STITCH at the 0.05 significance level."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
roc_fig_unit = 10
plt.figure(
    figsize=(roc_fig_unit, roc_fig_unit * fig_height / fig_width)
).subplot_mosaic([[1]])
lw = 2
for dti in ["duggie", "stitch"]:
    for annot in ["prox"]:
        df = hc_dti_df_sets[dtis.index(dti)][annots.index(annot)]
        if df["colour_id"].any():
            fpr, tpr, thresholds = roc_curve(
                (df["colour_id"] != 0), 1 - df.iloc[:, [3]]
            )
            roc_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, lw=lw, label=dti.upper() + " (AUC = %0.3f)" % roc_auc)

plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.plot([0.011, 0.011], [0, 0.125], color="black", lw=lw, linestyle="dashed")
plt.plot([0, 0.011], [0.125, 0.125], color="black", lw=lw, linestyle="dashed")
plt.plot([0.018, 0.018], [0, 0.269], color="black", lw=lw, linestyle="dashed")
plt.plot([0, 0.018], [0.269, 0.269], color="black", lw=lw, linestyle="dashed")
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.05])
plt.xlabel("Significant Non-treatment & Non-significant Treatment Rate")
plt.ylabel("Significant Treatment & Non-significant Non-treatment Rate")
plt.legend()
plt.show()
# -

rdisp.end_caption(formatting, r"fig:hc_prox_roc")

# Figure \ref{fig:hc_boxplot_treat_nontreat} shows the boxplot graph of the negative log of p-values, where the p-values were used to calculate the Mann-Whitney metrics, highlighting the ability of each drug target interaction database’s analysis to differentiate between the treatment and non-treatment drug groups. Again DUGGIE’s improved Mann-Whitney p-value is aided by the lower variation in association values obtained for the treatment drug set.

bold_caption = r"Box plot of hypercholesterolemia association p-values for treatment and non-treatment drug sets."
normal_caption = (
    r"A comparison between DUGGIE and STITCH using proximally annotated "
    r"genes associated with hypercholesteremia."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
fig, ax = plt.subplots(figsize=(10, 6))
duggie_df = hc_dti_df_sets[dtis.index("duggie")][annots.index("prox")]
stitch_df = hc_dti_df_sets[dtis.index("stitch")][annots.index("prox")]

duggie_df["-log(p)"] = -np.log10(duggie_df["P prox duggie"])
stitch_df["-log(p)"] = -np.log10(stitch_df["P prox stitch"])

plot_df = pd.DataFrame(
    {
        "DUGGIE treatment": duggie_df[duggie_df["colour_id"] != 0].loc[:, "-log(p)"],
        "DUGGIE nontreatment": duggie_df[duggie_df["colour_id"] == 0].loc[:, "-log(p)"],
        "STITCH treatment": stitch_df[stitch_df["colour_id"] != 0].loc[:, "-log(p)"],
        "STITCH nontreatment": stitch_df[stitch_df["colour_id"] == 0].loc[:, "-log(p)"],
    }
)

g = sns.boxplot(data=plot_df).set(xlabel="Drug group", ylabel="-log10(p-value)")
dummy = ax.set_xticklabels(
    labels=[
        "DUGGIE treatment",
        "DUGGIE non-treatment",
        "STITCH treatment",
        "STITCH non-treatment",
    ],
    rotation=30,
)
# -

rdisp.end_caption(formatting, r"fig:hc_boxplot_treat_nontreat")

# ### Hypercholesterolemia - correlation of STITCH and DUGGIE proximally annotated drug p-values 
#
# As for hypertension in section \ref{hypertension---correlation-of-stitch-and-duggie-proximally-annotated-drug-p-values}, figure \ref{fig:hc_correlation_qvals} shows the correlation of negative logs of the p-values of drugs present in both DUGGIE and STITCH is high with an r-value of 0.89, which is expected from the large number of gene targets shared between them. Similarly, treatment drugs often feature at the higher plotted values where there is a marked increase in significance of STITCH p-values, which is assumed to be a result of STITCH drugs already containing the most relevant targets with additional DUGGIE targets adding to more noise in the result.

bold_caption = (
    r"Hypercholesterolemia correlation plot of DUGGIE vs STITCH -log(p-value) for "
    r"proximally annotated drugs present in both databases."
)
normal_caption = (
    r"Blue markers are drugs listed in DrugBank as a treatment for the indication of "
    r"hypercholesterolemia. Other drugs are shown as pink markers. Drugs with significant "
    r"p-values for STITCH only, DUGGIE only, both DUGGIE and STITCH and none are shown as points, "
    r"down triangles, up triangles and diamonds respectively."
)
# TODO - key for these!
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic([[1]])
axlimit = 8

dti1 = "stitch"
dti2 = "duggie"

plt.subplot(fig_height, fig_width, 1)
rdisp.plot_drug_scatter(
    hc_gwas,
    "prox " + dti1,
    hc_dti_df_sets[dtis.index(dti1)][annots.index("prox")],
    "prox " + dti2,
    hc_dti_df_sets[dtis.index(dti2)][annots.index("prox")],
    "P",
    axlimit,
)
plt.title("")
plt.xlabel("Hypercholesterolemia STITCH drug assoc. -log10(p-value)")
dummy = plt.ylabel("Hypercholesterolemia DUGGIE drug assoc. -log10(p-value)")

# -

rdisp.end_caption(formatting, r"fig:hc_correlation_qvals")

# The relationship between drug gene target count and drug-disease association p-value is plotted in figure \ref{fig:hc_correlation_gene_counts}, where the number of target genes DUGGIE assigns in addition to those present in STITCH are split into four ranges and coloured independently. Again, a trend can be seen where a high gene count difference results in a larger disparity between drug p-values from the STITCH and DUGGIE analyses. Depending on the drug, this disparity can either be favouring the DUGGIE or STITCH analysis and in drugs which have a more significant result, the disparity favours STITCH.

bold_caption = (
    r"Hypercholesterolemia correlation plot of DUGGIE vs STITCH -log(p-value) for "
    r"proximally annotated drugs present in both databases."
)
normal_caption = (
    r"Dots are coloured differently based on the size of the difference between drug gene "
    r"target counts for DUGGIE and STITCH."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic([[1]])
axlimit = 8

dti1 = "stitch"
dti2 = "duggie"

diffs = []

for annot in ["prox"]:
    plt.subplot(fig_height, fig_width, annots.index(annot) + 1)
    df = rdisp.plot_drug_genediff_scatter(
        hc_gwas,
        annot + " " + dti1,
        hc_dti_df_sets[dtis.index(dti1)][annots.index(annot)],
        annot + " " + dti2,
        hc_dti_df_sets[dtis.index(dti2)][annots.index(annot)],
        "P",
        axlimit,
    )
    diffs.append(df)
    plt.xlabel("Hypercholesterolemia STITCH drug assoc. -log10(p-value)")
    plt.ylabel("Hypercholesterolemia DUGGIE drug assoc. -log10(p-value)")
    plt.title("")
# -

rdisp.end_caption(formatting, r"fig:hc_correlation_gene_counts")

# ### Hypercholesterolemia - correlation of p-value difference and drug gene count difference between DUGGIE and STITCH
#
# Looking at the variation of the negative log of p-value difference between the DUGGIE and STITCH analyses against drug target gene count difference between the two databases in Figure \ref{fig:hc_prox_diffs} shows that for hypercholesterolemia, even small differences in the gene counts added to DUGGIE have an adverse effect on the measured association with the disease. The situation appears to be improved at the higher numbers of genes added to DUGGIE, but still results in lower overall average p-values compared to STITCH.

bold_caption = (
    r"Hypercholesterolemia correlation plot of DUGGIE vs STITCH -log(p-value) difference vs gene "
    r"count difference for proximally annotated drugs present in both databases."
)
normal_caption = (
    r"Blue dots are anti-hypercholesterolemia drugs. The green regression line is for all drugs, "
    r"whilst the blue regression line is for anti-hypercholesterolemia drugs only."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

# +
# Diff on x-axis, p-val on y
plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width)).subplot_mosaic([[1]])

for i in range(0, 1):
    df = diffs[i]

    # -ve diff => STITCH better, +ve DUGGIE better
    df["p-diffs"] = df["prox duggie -log(P)"] - df["prox stitch -log(P)"]

    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df["diff"], df["p-diffs"]
    )

    df, group_labels, colours = rdisp.group_drugs(df, hc_gwas)

    groups = df.groupby("colour_id")

    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(
        df[df["colour_id"] == 1]["diff"], df[df["colour_id"] == 1]["p-diffs"]
    )

    plt.subplot(fig_height, fig_width, i + 1)

    for name, group in groups:
        plt.scatter(
            group["diff"],
            group["p-diffs"],
            label=group_labels[int(name)],
            c=colours[int(name)],
        )

    plt.xlabel("DUGGIE/STITCH drug gene difference")
    plt.ylabel("-log(p-value) DUGGIE/STITCH difference")

    ax2 = plt.gca()
    plt.text(
        0,
        -0.2,
        "r_value = "
        + "{0:.5f}".format(r_value)
        + ", slope = "
        + "{0:.5f}".format(slope)
        + ", err "
        + "{0:.5f}".format(std_err),
        transform=ax2.transAxes,
    )

    plt.text(
        0,
        -0.3,
        "r_value(treatment) = "
        + "{0:.5f}".format(r_value2)
        + ", slope = "
        + "{0:.5f}".format(slope2)
        + ", err "
        + "{0:.5f}".format(std_err2),
        transform=ax2.transAxes,
    )

    legend = ax2.legend([])
    plt.legend()
    plt.plot(df["diff"], slope2 * df["diff"] + intercept2, i)
    dummy = plt.plot(df["diff"], slope * df["diff"] + intercept, i)
# -

rdisp.end_caption(formatting, r"fig:hc_prox_diffs")

# ### Significant hypercholesterolemia associated drugs
#
# Table \ref{tab:hc_sig_drugs} lists all 24 significant drugs found by either DUGGIE or STITCH, ordered by DUGGIE q-value. All drugs listed are either treatment drugs or have a known effect on cholesterol, obtained from DrugBank's lists of anticholesteremic or lipid-modifying agents. All 5 treatment drugs found to be significantly associated with the disease by DUGGIE are also found significant by STITCH, and there are 6 DUGGIE and 5 STITCH significant drugs found not significant by the other's database analysis. Of these 6 DUGGIE significant drugs, only 1 is analysed by STITCH with the remaining 5 not analysed due to having less than the threshold of 4 gene targets assigned. The Venn diagrams showing these overlaps can be seen in figure \ref{fig:hc_all_venn}.

# +
proxdf = rdisp.get_annot_results(hc_gwas, "prox", magma_results_dir, run_id)
ATC_LEVELS = rdisp.get_atc_levels(data_path)

for dti in dtis:
    for annot in ["prox"]:
        rdisp.summarise_drug_results(
            hc_gwas,
            annot,
            ATC_LEVELS,
            magma_results_dir,
            summary_results_dir,
            dti,
            run_id,
            MAGMA_VER,
            1,
            1,
        )
# -

rdisp.display_joint_table_thesis(
    formatting,
    hc_gwas,
    run_id,
    "prox",
    "duggie_qvals.tsv",
    "stitch_qvals.tsv",
    "prox duggie",
    "prox stitch",
    proxdf,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
    r"Significant hypercholesterolemia drugs found by DUGGIE and STITCH. ",
    r"q-values in brackets are non-significant and included for completeness. "
    r"Drugs with no STITCH values are only present in the DUGGIE analysis.",
    "tab:hc_sig_drugs",
)

bold_caption = (
    r"Venn diagrams of overlap of (a) all significant drugs and (b) significant treatment drugs for "
    r"hypercholesterolemia."
)
normal_caption = r""
rdisp.start_caption(formatting, bold_caption, normal_caption)

rdisp.display_dgi_venn_thesis(
    hc_gwas,
    run_id,
    "prox",
    "duggie_qvals.tsv",
    "stitch_qvals.tsv",
    "prox duggie",
    "prox stitch",
    proxdf,
    summary_results_dir,
    genesets_df,
    gene_hgnc_df,
)

rdisp.end_caption(formatting, r"fig:hc_all_venn")

# ## Discussion
#
# The drug-gene interaction dataset DUGGIE has been shown to have an advantage over its largest contributing database, STITCH, in identifying drugs significantly associated with hypercholesterolemia and hypertension. Both are broadly comparable in having significant Fisher exact test p-values when finding treatment drugs and rejecting non-treatment drugs at a 0.05 significance level. Over all significance levels as displayed in the ROC curves, AUC values range between 0.74 and 0.94 with corresponding Mann-Whitney p-values found to be significant with exponents varying between -13 and -19 showing that both analyses are adept at correctly classifying drugs listed by DrugBank as treatments for the disease and rejecting non-treatment drugs.
#
# Of further interest is the observation that DUGGIE identifies more drugs as significantly associated than STITCH but does not include a greater number of treatment drugs among these significant drugs. One explanation for this could be that STITCH, being the older data set, has a bias towards including known and previously studied successful drug targets for diseases for which they are approved. That is, a successfully proven target gene is more likely to appear only in STITCH, as opposed to only in DUGGIE data not sourced from STITCH.
#
# Given this observation, the fact that STITCH is a subset of DUGGIE and that the comparison is a measurement of the difference made by extra targets added by DUGGIE in addition to those already in STITCH, these extra targets are of lower quality and would not be as effective in identifying treatment drugs for already extensively researched diseases. But conversely for diseases for which there are few known targets such as the neuropsychiatric diseases schizophrenia, Parkinson's disease and Alzheimer's disease that are studied later, this bias may well be less pronounced.
#
# Also of interest is the overlap in significantly associated drugs, treatment or non-treatment, between the DUGGIE and STITCH analysis results. There are many drugs which are only identified as significantly associated by one of the interaction database analyses, and not the other. It is thought that some DUGGIE results are reduced in significance due to the extra statistical noise introduced from both a larger number of genes annotated to each drug, and the larger set of drugs in DUGGIE able to take part in the competitive gene set analysis stage. Conversely, the extra coverage of the possible drug and gene target space by DUGGIE has led to some other drugs being available to DUGGIE which STITCH cannot analyse.
#
# Because of this difference in the nature of the two analyses - with STITCH being more targeted towards successful targets and DUGGIE having a wider range of drugs, targets and interactions, both could be useful in analysing neuropsychiatric diseases and will take part in the final analysis where two sets of pipeline runs will be performed, one using DUGGIE and another using STITCH.
