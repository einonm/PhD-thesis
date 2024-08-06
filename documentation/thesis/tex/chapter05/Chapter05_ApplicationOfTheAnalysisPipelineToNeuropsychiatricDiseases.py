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

# #!pip install pySankey
from pySankey.sankey import sankey

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

def count_drug_disease_genes(file_postfix, gwas, run_id, annot, dti):
    pval_file = [
                file
                for file in os.listdir(summary_results_dir)
                if file.endswith(file_postfix)
                and "found-" + gwas + "-" + run_id + "-" + annot + "-" in file
                ]

    # Take the first summary file, as the gene count is a property of the DGI-DB, not the annotation
    table = pd.read_csv(os.path.join(summary_results_dir, pval_file[0]), sep="\t")
    table.sort_values('Q ' + annot + ' ' + dti, inplace=True)
    table.set_index("ATC_CODE", inplace=True)

    # Dummy df
    gene_twas_df = pd.DataFrame()

    rdisp.add_genesets(table, zscores_df, genesets_df, gene_hgnc_df, gene_twas_df)

    return table


# # Application of the analysis pipeline to neuropsychiatric and neurodegenerative diseases
#
# ## Introduction
#
# Now that a gene set analysis pipeline has been validated and optimised to produce the best performing results, the pipeline is capable of being applied to neurodegenerative and psychiatric diseases that have little or no treatment drug options available with the aim of identifying novel repurposing opportunities.
#
# To recap, for each disease a two-stage gene set analysis pipeline was employed to produce a summary set of significant drug-disease associations, consisting firstly of a disease-gene annotation stage assigning GWAS summary SNP p-values to genes using one of a number of gene annotation schemes, followed by a competitive gene-set enrichment analysis stage between the set of annotated disease gene p-values and many drug gene sets taken from a drug-gene interaction database, this second stage outputting a disease association p-value for each drug.
#
# In optimising the pipeline, performance was measured principally in terms of the number of drugs found to be significantly associated with the disease and secondly by the number of treatment drugs found to be significantly associated with the disease. Following this, lower priority metrics measuring the performance of the pipeline in discerning between treatment and non-treatment drugs were used, in the form of Mann-Whitney and Fisher p-values - see table \ref{tab:metrics_hierarchy}. Using this approach, the best performing pipeline parameters were chosen by studying the effects of different drug-gene interaction databases (see chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}) and different gene annotation schemes (in chapter \ref{using-functional-qtl-gene-annotation-to-better-identify-genes-associated-with-disease}) on the number of drugs whose gene targets show significant enrichment for association with the disease. This resulted in the conclusion that the results of running a battery of pipeline executions using 16 annotation schemes, once each  with both the DUGGIE and STITCH drug-gene interaction databases giving 32 analyses in total, was the most successful in obtaining the maximum number of significantly-associated candidate drugs for repurposing. Furthermore, collating and ranking these drugs by the number of combinations of annotations and interaction databases in which the drugs were discovered to be significantly associated with the disease was found in the previous chapter to prioritise treatment over non-treatment drugs.

# ## Aims of the chapter
#
# This chapter presents the results of applying the analysis pipeline on several neuropsychiatric diseases, consisting of three psychiatric diseases - major depressive disorder, bipolar disorder and schizophrenia, and three neurodegenerative diseases - Parkinson's disease, Alzheimer's disease and Huntington's disease age of onset of motor symptoms. The drugs showing significant enrichment for association with disease were investigated to ascertain if any of them could be likely novel repurposing candidates, and for at least one identified novel repurposing candidate for each disease (the drug at the top of the ranked list), if the drug-gene interaction evidence behind its drug-disease associations is valid.
#
# A total of 32 gene set analysis pipeline executions were carried out for each disease, covering the combination of 16 different gene SNP annotations explored in chapter \ref{using-functional-qtl-gene-annotation-to-better-identify-genes-associated-with-disease} for both the DUGGIE and STITCH drug-gene interaction databases. Drugs which were found to be significantly associated with each disease were further classified as novel or not novel, depending on if a published study exists investigating the drug and disease interaction, or DrugBank listed the drug as an indication affected by the drug - in these cases the drug-disease relationship was deemed not novel. 
#
# After listing any novel drugs thus remaining, at least the top ranked novel drug from each disease list, that found to be significantly associated with the disease appearing in the greatest number of gene set analysis pipeline executions (out of the 32), was selected to validate the genes driving the association by examining the evidence supporting their link to the drug. This list could also be used to identify drug-relevant genes for further biological study - to test for biological pathways among these genes, for example.
#
# Further to this, for each disease the drugs found to be significantly associated by the most successful annotation and drug-gene interaction database combination were selected, i.e. the one giving the greatest number of drugs showing significant enrichment, and a tally of all gene targets of these drugs compiled. This tally was sorted by the count of drugs targeting the gene and a brief discussion given as to the relation of the genes to known disease mechanisms.
#
# Finally, the shared genetic architecture of the psychiatric diseases examined - schizophrenia, major depressive disorder and bipolar disorder<citet data-cite="grotzinger_shared_2021"><sup>grotzinger</sup></citet>, suggests that there would be an overlap of the significantly associated drugs identified for each, due to the potential existence of pleiotropic gene targets present within this group of diseases. Given that a drug which has been shown to be effective against one disease and shares a pleiotropic gene target with another disease may be prioritised as a drug repurposing candidate for this second disease<citet data-cite="omara_editorial_2019"><sup>omara</sup></citet>, the set of drugs that are identified for each disease was displayed on an UpSet plot, along with the overlaps between these sets. Any drugs which have evidence of being used to treat one of the diseases in this group will be noted as a repurposing candidate of interest for the other diseases. Similarly, the overlaps of drug targets from the most successful annotation for each disease was also shown on another UpSet plot in order to confirm that the expected overlap in targeted genes exists.
#
# ### Diseases
#
# #### Schizophrenia
#
# Schizophrenia is a condition primarily identified by the presence of psychosis, where experience or beliefs do not match reality, outside of cultural norms. The Diagnostic and Statistical Manual of Mental Disorders (5th Ed.)(DSM-5)<citet data-cite="noauthor_diagnostic_nodate"><sup>dsm5</sup></citet> lists the diagnosic criteria as being one of delusions, hallucinations or disorganised speech being present for at least six months, alongside an additional criterion from this list or one of disorganised/catatonic behaviour or another negative symptom such as depression, avolition or reduced emotional expression.
#
# Lifetime prevalence estimates of schizophrenia are given by the DSM-5 as being between 0.3% and 0.7%, and heritability has been assessed to be in the range 60%-80%<citet data-cite="bobo_diagnosis_2017"><sup>bpbobo</sup></citet>.
#
# #### Major depressive disorder (MDD)
#
# Major depressive disorder (MDD) is a complex, common mood disorder that the DSM-5<citet data-cite="noauthor_diagnostic_nodate"><sup>dsm5</sup></citet> defines as having five or more mainly negative symptoms from a given list of nine for more than two weeks, including depressed mood and loss of interest or pleasure. Other possible symptoms include significant weight change, continuous insomnia or hypersomnia, restlessness or slower movements, fatigue, feelings of worthlessness or inappropriate guilt, loss of concentration and recurrent thoughts of death.
#
# The DSM-5 also reports the twelve-month prevalence in the US as being approximately 7%, with significant variations across age and gender and a heritability estimate of around 40%. 
#
# #### Bipolar disorder
#
# Bipolar disorder is a mood and behavioural disorder characterised by distinct periods of mood instability, either where the mood is abnormally and persistently elevated or there is a period of depression, as described by the DSM-5<citet data-cite="noauthor_diagnostic_nodate"><sup>dsm5</sup></citet>, which further differentiates between three subtypes of bipolar disorder, referred to as bipolar I disorder, bipolar II disorder and cyclothymia, thus:
#
# Type I bipolar disorder is characterised by manic episodes, with the mood disturbance causing marked impairment in social or occupational functioning or with psychotic features. Type II bipolar disorder includes hypomanic episodes where there is an uncharacteristic change in behaviour during a period of abnormally elevated, expansive or irritable mood. Both bipolar I and II disorders involve major depressive episodes. Cyclothymia conversely, is diagnosed by an absence of either a major depressive or hypomanic episode but where there have been numerous hypomanic and depressive symptoms without reaching the full diagnostic criteria of either hypomanic or major depressive disorders.
#
# Lifetime prevalence of these bipolar disorders is estimated to be between 3% and 7% with heritability estimates of between 73% and 93%<citet data-cite="bobo_diagnosis_2017"><sup>bpbobo</sup></citet>.
#
# #### Alzheimer's disease
#
# Alzheimer's disease is a progressive neurodegenerative disorder and the leading cause of dementia, thought to account for 60%-70% of all cases<citet data-cite="noauthor_dementia_nodate"><sup>whodem</sup></citet>. At a rate over and above that which is expected due to ageing, starting with minor memory problems, the disease can progress to involve confusion and disorientation, difficulty planning, speech and language problems, issues in moving about and personality changes<citet data-cite="noauthor_alzheimers_2018"><sup>nhsdem</sup></citet>. 
#
# Alzheimer's is estimated to affect 1 in 14 people over the age of 65, and 1 in 6 people over the age of 80. Heritability is estimated to be between 60% and 80%<citet data-cite="bellenguez_new_2022"><sup>bellenguez22</sup></citet>.
#
# #### Parkinson's disease
#
# Parkinson's disease is a neurodegenerative disorder that symptomatically affects the motor system through cell death in the substantia nigra brain region, with involuntary tremors, slow movement and inflexible muscles being the main symptoms. Other symptoms include depression, balance issues, anosmia (loss of sense of smell), insomnia and memory problems<citet data-cite="noauthor_parkinsons_2017"><sup>nhs_parkinsons</sup></citet>.
#
# The disease prevalence of Parkinson's is estimated to be between 0.5% and 2%, with heritability estimated to be in the range 16% - 36%<citet data-cite="nalls_identification_2019"><sup>nalls19</sup></citet>.
#
# #### Huntington's disease
#
# Huntington's disease is an inherited neurodegenerative condition caused by variation of a single gene, the huntingtin gene (HTT), where the length of a CAG nucleotide trio repeat of 36 repeats or above gives a chance of the disease occurring. The age of onset depends on the CAG repeat length above this threshold with a CAG repeat length of 39 or above guaranteeing onset of the disease over a natural lifetime. Symptoms usually start between the ages of 30 and 50, with approximately half of the age of onset variation accounted for by CAG repeat length and approximately half of the remainder of the variation being heritable<citet data-cite="wexler_venezuelan_2004"><sup>wexler</sup></citet>. Symptoms progress from subtle mood and mental ability problems through coordination issues and involuntary movements to dementia and problems swallowing, speaking and breathing, with death typically occurring 15–20 years from when the disease was first detected<citet data-cite="noauthor_huntingtons_2022,noauthor_huntingtons_2017"><sup>huntingtons22,huntingtons17</sup></citet>.
#
# Prevalence of Huntington's disease is 10-15 cases per 100,000 persons amongst those of European decent. Whilst an inherited condition, up to 10% of cases are thought to be due to a new HTT CAG expansion into the disease range<citet data-cite="dayalu_huntington_2015"><sup>huntingtons_treatment</sup></citet>.
#
#

# ## Methods and materials
#
# ### Datasets
#
# ####  Schizophrenia GWAS summary statistics
#
# Schizophrenia GWAS summary statistics used to perform the gene set analysis identifying drugs significantly associated with the disorder<citet data-cite="sullivan_bip_2021"><sup>bpgwas</sup></citet> were taken from a recent PGC GWAS study<citet data-cite="trubetskoy_mapping_2022"><sup>trubetskoy</sup></citet>, involving up to 76,755 cases and 243,649 controls. This study associated 287 common variant genomic loci at the level of genome-wide significance with schizophrenia. The summary statistics arising from what the study identifies as the core meta-analysis for autosomal SNPs were used, obtained using 67,290 cases and 94,015 controls.
#
# #### Major depressive disorder GWAS summary statistics
#
# Major depressive disorder genome wide association study (GWAS) summary results, used as an input to the gene set analysis pipeline for identifying MDD associated drugs, were taken from the PGC Major Depression 2 GWAS results public data release<citet data-cite="patrick_sullivan_mdd_2018"><sup>mddgwas</sup></citet>. This GWAS analysis identified 44 independent genome-wide significantly associated genetic loci<citet data-cite="wray_genome-wide_2018"><sup>wray</sup></citet>, based on 135,458 cases and 244,901 controls. However, the GWAS summary dataset used for this study did not include 75,607 cases and 231,747 controls from 23andMe, which are not publicly available. 
#
# ####  Bipolar disorder GWAS summary statistics
#
# The bipolar disease gene set analysis GWAS summary results used<citet data-cite="sullivan_scz_2022"><sup>sczgwas</sup></citet> were taken from a recent GWAS study<citet data-cite="mullins_genome-wide_2021"><sup>mullins</sup></citet>, identifying 64 genomic loci as genome-wide significantly associated with bipolar disorder. The study used 41,917 bipolar I disorder and bipolar II disorder cases and 371,549 controls of European ancestry.
#
# #### Alzheimer's disease GWAS summary statistics
#
# Genome-wide summary statistics taken from a recent Alzheimer's study<citet data-cite="bellenguez_new_2022"><sup>bellenguez22</sup></citet> were used to identify drugs associated with Alzheimer's disease. The study took 111,326 Alzheimer's cases and 677,663 controls, identifying 75 genomic risk loci for Alzheimer's disease.
#
# ####  Parkinson's disease GWAS summary statistics
#
# Genome wide summary statistics from a recent meta analysis of Parkinson's disease<citet data-cite="nalls_identification_2019"><sup>nalls19</sup></citet> were used to investigate drug associations with the disease. Whilst the full study, involving 37,688 cases, 18,618 'proxy' cases of individuals with a first relative that has a Parkinson's diagnosis and 1,400,000 controls. From these, 90 independent genome-wide risk loci were identified. The GWAS summary statistics used for this study was the freely available version excluding data obtained from 23andMe.
#
# #### Huntington's disease age of onset GWAS summary statistics
#
# Huntington's is a monogenic disease, and whilst the HTT gene CAG repeat size is thought to explain ~50% of the variation in the age at onset of the disease, approximately half of the remainder is heritable<citet data-cite="wexler_venezuelan_2004"><sup>wexler</sup></citet>. Therefore, GWAS have been performed on age at onset to identify the genes and variants that account for this heritability. The summary statistics from the most recent GWAS study of Huntington's age at onset<citet data-cite="genetic_modifiers_of_huntingtons_disease_gem-hd_consortium_electronic_address_gusellahelixmghharvardedu_cag_2019"><sup>hunt2019</sup></citet> were used as input to a gene set analysis to identify drugs that may delay the age of onset of Huntington's disease.

# ### Analysis methods
#
# #### Generating drug repurposing candidates and ascertaining novelty
#
# The gene-set analysis pipeline utilising the MAGMA tool, as described in section \ref{a-gene-set-enrichment-analysis-pipeline} was once more employed, with 32 pipelines executed for each disease differing in the input drug-gene interaction database used - either DUGGIE or STITCH, which are both discussed in chapter \ref{curation-of-a-drug-target-interaction-database} - or the SNP to gene annotation method used, covering positional, functional and hybrid approaches which are listed in table \ref{tab:db_annotations}. Any drugs found to be significantly associated with the disease in any of the pipeline executions were pooled and ranked by the number of analyses in which they appear as significant - a method which when performed for hypercholesterolemia and hypertension in chapter \ref{using-functional-qtl-gene-annotation-to-better-identify-genes-associated-with-disease}, was seen to favour known treatment drugs.
#
# Following this generation of significantly associated drugs for each disease, a classification of each drug as a novel or non-novel repurposing candidate was made. The information to drive this classification was gathered from the DrugBank web site<citet data-cite="law_drugbank_2014"><sup>law</sup></citet> and literature searches, conducted using the search functionality of the PubMed web site<citet data-cite="noauthor_pubmed_nodate"><sup>pubmed</sup></citet> using the boolean search string ‘(drug) AND (disease)’, replacing the 'drug' and 'disease' terms with the relevant items as necessary. A lack of evidence arising from these activities identified drugs that were classed as novel. Drugs are classified from this search by finding evidence for each drug as:
#
# * a known repurposing candidate for the disease, where a study or randomised control trial (RCT) has directly investigated the drug's effectiveness against the disease, not including case studies.
# * a known ineffective drug against the disease, where studies (excluding case studies) or RCTs have investigated the drug's effectiveness against the disease and few or no supportive findings have been made compared to negative findings.
# * a known existing treatment drug for the disease, where DrugBank lists the disease as an 'indication' or 'associated condition' for the drug. DrugBank also lists indications/associated conditions of withdrawn drugs, which would also appear under this classification.
# * a new, novel drug whose effectiveness in treating the disease is yet to be investigated if no link was found between the disease and the drug from literature searches and DrugBank. 
#
# If a published study exists investigating the drug and disease interaction or DrugBank lists the drug as an indication affected by the drug, then the drug-disease relationship was deemed not novel. This method is not an exhaustive search for evidence supporting the novelty of a repurposing candidate, but it does provide a quick method of discarding a large number of drugs that have been considered previously as repurposing candidates.
#
# #### Gene targets of significantly associated drugs
#
# The association Z-scores assigned to each gene are dependent on the particular annotation method, and the set of genes assigned to a drug dependent on the drug-gene interaction database used. This means that compiling a list of targets of successful drugs is not straightforward, as the list of potential drug repurposing candidates is compiled from results using sixteen different annotation methods and two drug-target interaction databases. Therefore, for each disease, the annotation and drug-target interaction database combination - which showed the highest number of significantly disease-associated drugs was chosen, and the list of targets for these drugs was used to compile the target list for that disease, ordered in descending count of the number of significantly associated drugs which target that gene. The gene's chromosome, start position and Z-score association with the disease were also shown.
#
# #### Comparing overlaps of significantly associated drugs and their gene targets 
#
# UpSet plots<citet data-cite="lex_upset_2014"><sup>upset</sup></citet> showing the overlaps, intersections and aggregates of the sets of drugs identified as significantly associated with each disease were plotted. This plot used the drug results from all 32 pipeline executions before this list was reduced to just the potential novel drug repurposing candidates following the outcome of the literature search. 
#
# A second UpSet plot of target genes was created, this time using the gene lists compiled to create the gene target tables for each disease described in the previous section - from the most successful annotation and drug-target database combination for each disease. This means that the data used for this second gene target UpSet plot is not identical to that used for the drug UpSet plot, but it does give the largest number of targets available for each disease that are assigned and associated in a consistent manner for that disease, i.e. consistently using the same annotation and drug-gene interaction database.
#
# A drug or target found to be significantly associated with two or more diseases which are understood to have a shared genetic architecture may be basing this overlap upon causal genes and pathways shared between the diseases. If one of these drugs is a treatment drug for one disease in this set, then the likelihood that the drug is effective against the other diseases is greater, and hence would have a more compelling reason to be ranked higher than other drugs or targets which are not linked to the diseases in this way.

# ## Results

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


# All STITCH + ATC genesets, remove STITCH duplicates
genesets_df = pd.concat([atc_genesets_df, stitch_genesets_df])
genesets_df.reset_index(inplace=True)
genesets_df = genesets_df.drop_duplicates(subset=[0])
genesets_df.set_index(0, inplace=True)

# Emsembl to HGNC translation table
gene_hgnc_df = pd.read_csv(
    data_path + "wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    header=None,
    names=["Ensembl", "1", "chrom", "chromStart", "chromEnd", "HGNC"],
)[["Ensembl", "HGNC", "chrom"]]

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
        
ATC_LEVELS = rdisp.get_atc_levels(data_path)

db_indications_df = pd.DataFrame.from_dict(rdisp.indications, orient='index')
# -

# ### Results - Schizophrenia

# +
gwas = "PGC3_SZ"
run_id = "3qtl"

dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        dti_df_sets[dtis.index(dti)].append(newlist)
# -

sz_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, 
                                           gwas, run_id, SIGNIF_THRESH)

# #### Schizophrenia - drug enrichment across annotation results

# In total, 109 drugs were identified as having a significant association with schizophrenia identified genes in at least one pipeline analysis run. The full list can be seen in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results}, table \ref{tab:sz_most_frequent_drugs}. From these drugs, 26 of which were classed as novel following the literature search, table \ref{tab:sz_most_frequent_novel_drugs} lists the top 20 novel drugs by number of analyses in which they were found to be significantly associated with schizophrenia. A majority of these drugs are approved to treat hypertension and related indications, showing the prominence of the CACNA1C gene as a drug target that is highly genome-wide associated with schizophrenia. Despite CACNA1C being a frequently investigated target for schizophrenia, many drugs with this gene target do not have studies directly investigating their effect on the disease and so appear as novel drugs in these results.

sz_top_drugs_df = rdisp.get_top_drugs_ch5(sz_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.sz_novel_drugs)

# Save drug names for publication search script
sz_top_drugs_df['Description'].to_csv(data_path + "wrangling/sz_drugnames.txt", index=False, header=None)

sz_top_novel_drugs = sz_top_drugs_df[sz_top_drugs_df['Novel'] == 'yes'].drop(['Novel'], axis=1)
sz_top_novel_drugs.rename(columns={'Description':'Drug name'}, inplace=True)
sz_top_novel_drugs.index.rename('ATC', inplace=True)

if formatting == "HTML":
    display(HTML(sz_top_novel_drugs.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            sz_top_novel_drugs.head(20).to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5cm}R{8cm}",
                caption=r"Top 20 novel drugs, ordered by frequency of appearance as significantly "
                        r"associated with schizophrenia across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated, with drug name and drug indications taken from DrugBank.",
                label="tab:sz_most_frequent_novel_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### The evidence behind an example novel drug for schizophrenia
#
# To select an example novel drug for schizophrenia, the list in table \ref{tab:sz_most_frequent_novel_drugs} was traversed from the highest analysis count downwards and the evidence previously supporting each drug-gene target association was reviewed. The majority of drugs listed were found to have a top gene target (the gene target with the largest Z-score association with the disease) of CACNA1C, a calcium channel gene <citet data-cite="trubetskoy_mapping_2022"><sup>trubetskoy</sup></citet> and also a widely studied drug target for treating the disease <citet data-cite="harrison_voltage-gated_nodate"><sup>harrison</sup></citet>. The first CACNA1C targeting novel drug is magnesium sulfate, found to be associated with schizophrenia in 17 of the 32 analysis runs using both the DUGGIE and STITCH drug-gene interaction databases. Both DUGGIE and STITCH list the same set of six gene targets, obtained from the STITCH and DrugBank databases - CACNA1C, CACNB2, CACNA2D1, CACNG1, CACNA1S and CACNB1, all members of the voltage-gated calcium channel superfamily of genes and listed in table \ref{tab:sz_top_drug_genes1}. The evidence for these targets are listed on the relevant DrugBank page<citet data-cite="noauthor_magnesium_nodate"><sup>db_magnesiumsulfate</sup></citet>, which also states that magnesium sulfate is used to control convulsions caused by eclampsia in pregnancy and treating acute nephritis in children. It is also associated with the treatment of constipation, hypomagnesia and Torsades de Pointes heart arrhythmias.
#
# In additiion to CACNA1C, of the five other calcium channel targets listed, one has also been implicated directly from schizophrenia GWAS studies, CACNB2<citet data-cite="schizophrenia_working_group_of_the_psychiatric_genomics_consortium_biological_2014"><sup>pgc2</sup></citet>. Other studies have discussed evidence linking CACNA2D1 <citet data-cite="andrade_genetic_2019"><sup>andrade</sup></citet>, CACNB1 and CACNA1S<citet data-cite="heyes_genetic_2015"><sup>heyes</sup></citet>, whilst no literature was found linking CACNG1 to schizophrenia.

# +
annot = 'prox'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
sz_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')

# +
sz_top_drug1_df = rdisp.add_genesets_thesis(sz_table[sz_table.index == 'A06AD04'], 
                                           zscores_df, 
                                           genesets_df, 
                                           gene_hgnc_df)

sz_top_drug1_df.drop(['Num genes', 'Q', 'Description', 'ATC code'], axis=1, inplace=True)
sz_top_drug1_df.rename(columns={'ZSTAT':'Z'}, inplace=True)
sz_top_drug1_df.set_index('HGNC', inplace=True)
sz_top_drug1_df['Z'] = sz_top_drug1_df['Z'].astype('float').round(3)
# -

sz_top_drug1_df['protein'] = [
    'Calcium Voltage-Gated Channel Subunit Alpha1 C',
    'Calcium Voltage-Gated Channel Auxiliary Subunit Beta 2',
    'Calcium Voltage-Gated Channel Subunit Alpha1 S',
    'Calcium Voltage-Gated Channel Auxiliary Subunit Alpha2delta 1',
    'Calcium Voltage-Gated Channel Auxiliary Subunit Beta 1',
    'Calcium Voltage-Gated Channel Auxiliary Subunit Gamma 1',
]
sz_top_drug1_df['evidence sources'] = [
    'DrugBank, STITCH, T3DB',
    'DrugBank, STITCH',
    'DrugBank, STITCH',
    'DrugBank, STITCH',
    'DrugBank, STITCH',
    'DrugBank, STITCH',
]

if formatting == "HTML":
    display(HTML(sz_top_drug1_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            sz_top_drug1_df.to_latex(
                column_format=r"p{1.5cm}p{1.5cm}R{7cm}R{4cm}",
                caption=r"Genes interacting with the drug magnesium sulfate according to the DUGGIE "
                        r"drug-gene interaction database with gene Z-score associations with schizophreina "
                        r"given for the "
                        r"'Proximity' annotation and validated sources of evidence used to compile DUGGIE.",
                label="tab:sz_top_drug_genes1",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# Of interest also are the non-CACNA1C targeting drugs in the list which may highlight novel drug targets for schizophrenia, the first three of which are cisplatin, minaprine and gutathione. Cisplatin is a chemotherapy drug approved to treat many cancers - for example testicular, ovarian and bladder cancers, acting on cells that divide quickly, for which there are significant gene-set analysis associations obtained using both DUGGIE and STITCH. The top 10 gene targets of cisplatin by Z-score out of 24 stored in DUGGIE are listed in table \ref{tab:sz_top_drug_genes2}. 
#
# Cisplatin's gene target HIST2H4A has the most significant Z-score schizophrenia association, with the evidence for the target originating from the STITCH drug-target interaction database. STITCH assigns an experimental evidence based score of 0.8 (out of a maximum of 1, STITCH indicates that a score > 0.7 has a high confidence) to the HIST1H4H gene expressed protein histone H4 in binding with cisplatin. Histone H4 is encoded by several paralog genes, of which two are HIST2H4A and HIST1H4H. The gene set analysis pipeline assigns the HIST2H4A gene to cisplatin in this case, as the Entrez ID for the gene (8370) is used in the conversion, which forward maps by Entrez to the HGNC symbol HIST2H4A only, but reverse maps several of the histone genes to the one Entrez ID, 8370. Schizophrenia is thought to be associated with some epigenetic histone modifications<citet data-cite="chen_research_2021"><sup>chen</sup></citet>, but the interaction of a histone targeting drug with methylation would need further investigation. 
#
# However, the therapeutic mechanism of action of cisplatin is the prevention of DNA replication, including of cancer cells<citet data-cite="dasari_cisplatin_2014"><sup>dasari</sup></citet>, and along with common serious side effects such as liver damage (nephrotoxicity), bone marrow failure and arrhythmias<citet data-cite="noauthor_cisplatin_nodate"><sup>nice-sisplatin</sup></citet> that make this drug unsuitable for repurposing to treat schizophrenia. The analyses undertaken in this study do not use such evidence in identifying drug repurposing opportunities and it it likely that most candidate drugs may be discounted following further consideration, an exercise which is beyond the scope of this study but an essential future step in gauging the feasibility of a repurposing opportunity.

# \newpage

# +
annot = 'alleqtlbrainhybridboth'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
sz_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')

# +
sz_top_drug2_df = rdisp.add_genesets_thesis(sz_table[sz_table.index == 'L01XA01'], 
                                           zscores_df, 
                                           genesets_df, 
                                           gene_hgnc_df)

sz_top_drug2_df.drop(['Num genes', 'Q', 'Description', 'ATC code'], axis=1, inplace=True)
sz_top_drug2_df.rename(columns={'ZSTAT':'Z'}, inplace=True)
sz_top_drug2_df.set_index('HGNC', inplace=True)
sz_top_drug2_df['Z'] = sz_top_drug2_df['Z'].astype('float').round(3)
# -

sz_top_drug2_df.drop('A2M', inplace=True)

sz_top_drug2_df['protein'] = [
    'H4 Clustered Histone 14',
    'Mitogen-Activated Protein Kinase 3', 
    'Death Domain Associated Protein',
    'FYN Proto-Oncogene, Src Family Tyrosine Kinase',
    'Thromboxane A Synthase 1',
    'DNA Polymerase Eta',
    'Melanocortin 5 Receptor',
    'Lysozyme',
    'ATP Binding Cassette Subfamily C Member 5',
    'ATP Binding Cassette Subfamily C Member 2',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
]
sz_top_drug2_df['evidence sources'] = [
    'STITCH',
    'DrugCentral', 
    'STITCH',
    'DrugCentral',
    'DrugCentral',
    'STITCH',
    'DrugCentral',
    'STITCH',
    'DrugCentral',
    'DrugCentral',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
]

if formatting == "HTML":
    display(HTML(sz_top_drug2_df.head(10).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            sz_top_drug2_df.head(10).to_latex(
                column_format=r"p{1.5cm}p{1.5cm}p{7cm}p{4cm}",
                caption=r"Top 10 genes by Z-score interacting with the drug cisplatin according to "
                        r"the DUGGIE drug-gene interaction database with gene Z-score associations with "
                        r"schizophrenia given for the "
                        r"'Brain-combined' annotation and validated sources of evidence used to compile DUGGIE.",
                label="tab:sz_top_drug_genes2",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### Gene target prevalence across drugs significantly associated with schizophrenia
#
# From Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results} table \ref{tab:sz_cmatrix_drugs}, the annotation giving the largest number of significant drugs (79) is 'Brain-mQTL-comb. DUGGIE', which uses psychEncode eQTL and mQTL data with position data added in the combination hybrid approach (Described as 'Brain-mQTL-combined' in table \ref{tab:db_annotations}), with the DUGGIE drug-gene interaction database. Using this drug-gene interaction database as a reference set of drugs and drug targets, the top 20 most commonly occurring gene targets along with their Z-score association with schizophrenia and count of drugs targeting each gene from the results using this annotation set are listed in table \ref{tab:sz_drug_gene_targets}. Calcium channel genes dominate this list, with CACNA1C appearing as the fourth most frequent target, below three alpha adrenergic receptor genes (ADRAxx), which predominantly co-occur as a target for drugs that also target CACNA1C, as well as being targeted by a few other non-CACNA1C targeting drugs.

# +
annot = 'allbrainhybridboth'
file_postfix = 'duggie_qvals.tsv'
zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)

sz_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')
sz_table2 = rdisp.tally_genesets_thesis(sz_table, zscores_df, genesets_df, gene_hgnc_df)
# Remove any X-chromosome targets (with no Z-score)
sz_table2.dropna(inplace=True)
sz_table2 = sz_table2.sort_values('count', ascending=False)
sz_table2 = sz_table2[['count', 'chrom', 'start_pos', 'Z']]
sz_table2['Z'] = sz_table2['Z'].astype('float').round(3)
# -

if formatting == "HTML":
    display(HTML(sz_table2.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            sz_table2.head(20).to_latex(
                column_format=r"p{2cm}p{1cm}p{1cm}p{3cm}p{2cm}",
                caption=r"Top 20 drug gene targets for drugs significantly associated with "
                        r"schizophrenia in results from a gene set analysis using the brain-mQTL-combined "
                        r"annotation and DUGGIE drug-gene interaction database. "
                        r"Each gene target is listed along with the gene chromosome and position, count of "
                        r"significantly associated drugs which target the gene and gene-disease association "
                        r"Z-score resulting from the brain-mQTL-combined MAGMA SNP annotation stage.",
                
                label="tab:sz_drug_gene_targets",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
# Add all genes to a set, with counts
#sankey_df = pd.DataFrame(columns=['drug', 'gene'])

#for index, row in bp_table.iterrows():
#    for gene in row['HGNC'].split('\n'):
#        atc_group = index[0:3]
#        sankey_df = sankey_df.append({'drug':atc_group, 'gene':gene}, ignore_index=True)

# +
#sankey(sankey_df["drug"], sankey_df["gene"], aspect=10)
#fig = plt.gcf()
#fig.set_size_inches(30, 50)
# -

# ### Results - Major depressive disorder

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
# -

mdd_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, gwas, run_id, SIGNIF_THRESH)

# #### Major depressive disorder - drug enrichment across annotation results
#
# Over all annotations used to perform gene set analyses with major depressive disorder, 49 drugs were found to have a significant association with the disease. The complete list is given in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results}, table \ref{tab:mdd_most_frequent_drugs}. Of these, 20 were found to not have any dedicated studies linking their effectiveness in treating major depressive disorder and are listed in table \ref{tab:mdd_most_frequent_novel_drugs} as novel drugs. 

mdd_top_drugs_df = rdisp.get_top_drugs_ch5(mdd_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.mdd_novel_drugs)

mdd_top_novel_drugs = mdd_top_drugs_df[mdd_top_drugs_df['Novel'] == 'yes'].drop(['Novel'], axis=1)

# Save drug names for publication search script
mdd_top_drugs_df['Description'].to_csv(data_path + "wrangling/mdd_drugnames.txt", index=False, header=None)

if formatting == "HTML":
    display(HTML(mdd_top_novel_drugs.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            mdd_top_novel_drugs.head(20).to_latex(
                column_format=r"p{1cm}p{0.8cm}p{5cm}p{8cm}",                                                                                                                               
                caption=r"Top 20 novel drugs, ordered by frequency of appearance as significantly "
                        r"associated with major depressive disorder across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated, with drug name and drug indications taken from DrugBank.",
                label="tab:mdd_most_frequent_novel_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### The evidence behind an example novel drug for major depressive disorder
#
# The top novel drug, that which appears most frequently in the different annotation analyses, found to be associated with major depressive disorder listed in table \ref{tab:mdd_most_frequent_novel_drugs} is triflupromazine, which DrugBank states is used mainly in the management of psychoses and for the control of nausea and vomiting<citet data-cite="noauthor_triflupromazine_nodate"><sup>db_triflupromazine</sup></citet>, as a dopamine D1/D2 receptor antagonist. Receptor antagonists bind to a receptor, interfering with the receptor function and thus dampens the normal biological response of the receptor.
#
# DUGGIE lists this association as driven by a set of eleven genes which are shown in table \ref{tab:mdd_top_drug_genes} in descending order of the gene-disease association Z-score, with the verified sources of evidence also indicating that the drug acts as an antagonist to the majority of proteins coded by these genes.

# +
annot = 'alltissues'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
mdd_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')

# +
mdd_top_drug_df = rdisp.add_genesets_thesis(mdd_table[mdd_table.index == 'N05AA05'], 
                                            zscores_df, 
                                            genesets_df, 
                                            gene_hgnc_df)

mdd_top_drug_df.drop(['Num genes', 'Q', 'Description', 'ATC code'], axis=1, inplace=True)
mdd_top_drug_df.rename(columns={'ZSTAT':'Z'}, inplace=True)
mdd_top_drug_df.set_index('HGNC', inplace=True)
mdd_top_drug_df['Z'] = mdd_top_drug_df['Z'].astype('float').round(3)
# -

mdd_top_drug_df['protein'] = [
    'Cholinergic Receptor Muscarinic 2',
    'Serotonin Receptor 2A',
    'Tumor Protein P53',
    'Dopamine Receptor D2',
    'Cytochrome P450 1A2',
    'Serotonin Receptor 2B',
    'ATP Binding Cassette Subfamily B Member 1',
    'Cholinergic Receptor Muscarinic 1',
    'Cytochrome P450 2D6',
    'Dopamine Receptor D1',
    'G Protein-Coupled Receptor 55 '
]
mdd_top_drug_df['evidence sources'] = [
    'DrugBank, STITCH, DGIdb',
    'PDSP (via DrugCentral, STITCH)',
    'DsigDB',
    'DrugBank, STITCH, DGIdb, DSigDB, DrugCentral',
    'DsigDB',
    'DrugBank, STITCH, DGIdb, TTD, T3DB',
    'DrugCentral',
    'DrugBank, STITCH, DGIdb ',
    'DsigDB',
    'DsigDB',
    'DGIdb',
]

if formatting == "HTML":
    display(HTML(mdd_top_drug_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            mdd_top_drug_df.to_latex(
                column_format=r"p{1.2cm}p{1.2cm}p{5.6cm}p{6.5cm}",
                caption=r"Genes interacting with the drug triflupromazine according to the DUGGIE "
                        r"drug-gene interaction database with gene Z-score associations with major depressive "
                        r"disorder given for the "
                        r"'AlleQTLs-mQTL' annotation and validated sources of evidence used to compile DUGGIE.",
                label="tab:mdd_top_drug_genes",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### Gene target prevalence across drugs significantly associated with major depressive disorder
#
# The annotation giving the largest number of significant drugs (21), using positional data ('Proximity' from Table \ref{tab:db_annotations}) and the DUGGIE drug-gene interaction database is used as a source of drug-target and gene-disease association Z-score information for MDD, with the full list of statistics for all annotations available in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results} table \ref{tab:mdd_cmatrix_drugs}. The top 20 genes from this group of drugs and targets are shown in table \ref{tab:mdd_drug_gene_targets}, along with frequency counts of drugs targeting the genes and association Z-scores obtained with this annotation. Serotonin receptor genes (also known as 5-hydroxytryptamine receptors, hence their identification as HTRxx) appear heavily, and also appearing are the DRDx genes which encode dopamine receptor proteins. Both serotonin and dopamine are well known to have a role in depression<citet data-cite="esposito_serotonindopamine_2008"><sup>esposito</sup></citet>.

# \newpage

# +
annot = 'prox'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
mdd_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')
mdd_table2 = rdisp.tally_genesets_thesis(mdd_table, zscores_df, genesets_df, gene_hgnc_df)
# Remove any X-chromosome targets (with no Z-score)
mdd_table2.dropna(inplace=True)
mdd_table2 = mdd_table2.sort_values('count', ascending=False)
mdd_table2 = mdd_table2[['count', 'chrom', 'start_pos', 'Z']]
mdd_table2['Z'] = mdd_table2['Z'].astype('float').round(3)
# -

if formatting == "HTML":
    display(HTML(mdd_table2.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            mdd_table2.head(20).to_latex(
                column_format=r"p{2cm}p{1cm}p{1cm}p{3cm}p{2cm}",
                caption=r"Top 20 drug gene targets for drugs significantly associated with "
                        r"major depressive disorder in results from a gene set analysis using the "
                        r"proximity annotation and DUGGIE drug-gene interaction database. "
                        r"Each gene target is listed along with the gene chromosome and position, count of "
                        r"significantly associated drugs which target the gene and gene-disease association "
                        r"Z-score resulting from the proximity MAGMA SNP annotation stage.",
                label="tab:mdd_drug_gene_targets",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Results - Bipolar disorder

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
# -

bp_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, gwas, run_id, SIGNIF_THRESH)

# #### Bipolar disorder - drug enrichment across annotation results
#
# The bipolar disorder gene set analysis found the largest number of significantly associated drugs, 116, out of all six diseases. The full list of 49 novel and 67 non-novel significant drugs can be found in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results}, table \ref{tab:bp_most_frequent_drugs}. Table \ref{tab:bp_most_frequent_novel_drugs} gives as a subset of these the top 20 novel drugs by frequency of appearance in the results of all 32 pipeline executions, differing by the gene annotation and drug-gene interaction database used. As for schizophrenia, magnesium sulfate is the top novel drug with a top target of CACNA1C, having the highest Z-score of association with the disease. And also mirroring the schizophrenia results, CACNA1C targeting drugs for hypertension and related indications feature heavily in the results suggesting an overlap of gene sets associated with both diseases.

bp_top_drugs_df = rdisp.get_top_drugs_ch5(bp_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.bp_novel_drugs)
bp_top_novel_drugs_df = bp_top_drugs_df[bp_top_drugs_df['Novel'] == 'yes'].drop(['Novel'], axis=1)

# Save drug names for publication search script
bp_top_drugs_df['Description'].to_csv(data_path + "wrangling/bp_drugnames.txt", index=False, header=None)

if formatting == "HTML":
    display(HTML(bp_top_novel_drugs_df.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            bp_top_novel_drugs_df.head(20).to_latex(
                column_format=r"p{1cm}p{0.8cm}p{5cm}p{8cm}",                                                                                                                               
                caption=r"Top 20 novel drugs, ordered by frequency of appearance as significantly "
                        r"associated with bipolar disorder across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated, with drug name and drug indications taken from DrugBank.",
                label="tab:bp_most_frequent_novel_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### The evidence behind an example novel drug for bipolar disorder
#
# The most frequent novel drug found to be significantly associated with bipolar disorder using both the DUGGIE and STITCH drug-gene interaction databases was magnesium sulfate, a drug approved to treat convulsions and the evidence behind which is discussed in the schizophrenia results above. As other novel non-CACNA1C drugs may also yield interesting targets, the first four such drugs in table \ref{tab:bp_most_frequent_novel_drugs} are ethchlorvynol, pentobarbital, secobarbital and hexobarbital.
#
# The most frequently associated drug in the set of 32 gene set analyses, ethchlorvynol, targets a protein group of Gamma-aminobutyric acid, or GABA, receptors. The drug is significantly associated with bipolar disorder in pipeline runs involving only the STITCH drug-gene interaction database. For analyses involving DUGGIE, a few more gene targets are assigned and also the competitive nature of the gene-set analyses result in less than significant Q-values being given to the ethchlorvynol drug-disease association. Eight GABA genes as targets of ethchlorvynol - GABRA2, GABRA4, GABRB2, GABRA6, GABRB1, GABRA1, GABRB3, GABRA5 as listed in table \ref{tab:bp_top_drug_genes}. The other three non-CACNA1C targeting drugs, all barbituates, also mainly target this GABA receptor protein group and were found to be significantly associated with bipolar disorder using both DUGGIE and STITCH. GABA abnormalities have previously been implicated in bipolar disorder<citet data-cite="brady_brain_2013"><sup>brady</sup></citet>.

# +
annot = 'alleqtlbrainhybrid'
file_postfix = 'stitch_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
bp_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'stitch')

# +
bp_top_drug_df = rdisp.add_genesets_thesis(bp_table[bp_table.index == 'N05CM08'], 
                                            zscores_df, 
                                            stitch_genesets_df, 
                                            gene_hgnc_df)

bp_top_drug_df.drop(['Num genes', 'Q', 'Description', 'ATC code'], axis=1, inplace=True)
bp_top_drug_df.rename(columns={'ZSTAT':'Z'}, inplace=True)
bp_top_drug_df.set_index('HGNC', inplace=True)
bp_top_drug_df['Z'] = bp_top_drug_df['Z'].astype('float').round(3)
# -

bp_top_drug_df['protein'] = [
    'Gamma-Aminobutyric Acid Receptor Subunit Alpha-6',
    'Gamma-Aminobutyric Acid Receptor Subunit Alpha-2',
    'Gamma-Aminobutyric Acid Receptor Subunit Beta-1',
    'Gamma-Aminobutyric Acid Receptor Subunit Beta-3',
    'Gamma-Aminobutyric Acid Receptor Subunit Alpha-1',
    'Gamma-Aminobutyric Acid Receptor Subunit Alpha-4',
    'Gamma-Aminobutyric Acid Receptor Subunit Alpha-5',
    'Gamma-Aminobutyric Acid Receptor Subunit Beta-2',
]
bp_top_drug_df['evidence sources'] = [
    'STITCH, DrugBank, DGIdb',
    'STITCH, DrugBank, DGIdb',
    'STITCH, DrugBank, DGIdb',
    'STITCH, DrugBank, DGIdb',
    'STITCH, DrugBank, DGIdb',
    'STITCH, DrugBank, DGIdb',
    'STITCH, DrugBank, DGIdb',
    'STITCH, DrugBank, DGIdb',
]

if formatting == "HTML":
    display(HTML(bp_top_drug_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            bp_top_drug_df.to_latex(
                column_format=r"p{1.2cm}p{1.2cm}p{7.2cm}p{5cm}",
                caption=r"Genes interacting with the drug ethchlorvynol according to the STITCH "
                        r"drug-gene interaction database with gene Z-score associations with bipolar disorder "
                        r" given for the 'Brain-missing' annotation and validated sources of evidence used.",
                label="tab:bp_top_drug_genes",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### Gene target prevalence across drugs significantly associated with bipolar disorder
#
# The DUGGIE drug-gene interaction database alongside the psychEncode eQTL annotation with position data added in the missing hybrid approach ('Brain-missing' from Table \ref{tab:db_annotations}) resulted in the largest number of significant drugs found (74) and is used as a source of drug-target and gene-disease Z-score association information for bipolar disorder. The complete list of such results for all annotations can be seen in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results} table \ref{tab:bp_cmatrix_drugs}. The top 20 most frequently occurring gene targets from this results set is shown in table \ref{tab:bp_drug_gene_targets}, and closely matches the equivalent target list obtained for schizophrenia (table \ref{tab:sz_drug_gene_targets}) showing the dominance in the results of calcium channel genes, particularly CACNA1C, and co-occurring gene targets.

# +
annot = 'alleqtlbrainhybrid'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
bp_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')
bp_table2 = rdisp.tally_genesets_thesis(bp_table, zscores_df, genesets_df, gene_hgnc_df)
# Remove any X-chromosome targets (with no Z-score)
bp_table2.dropna(inplace=True)
bp_table2 = bp_table2.sort_values('count', ascending=False)
bp_table2 = bp_table2[['count', 'chrom', 'start_pos', 'Z']]
bp_table2['Z'] = bp_table2['Z'].astype('float').round(3)
# -

if formatting == "HTML":
    display(HTML(bp_table2.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            bp_table2.head(20).to_latex(
                column_format=r"p{2cm}p{1cm}p{1cm}p{3cm}p{2cm}",
                caption=r"Top 20 drug gene targets for drugs significantly associated with "
                        r"bipolar disorder in results from a gene set analysis using the brain-missing "
                        r"annotation and DUGGIE drug-gene interaction database. "
                        r"Each gene target is listed along with the gene chromosome and position, count of "
                        r"significantly associated drugs which target the gene and gene-disease association "
                        r"Z-score resulting from the brain-missing MAGMA SNP annotation stage.",
                label="tab:bp_drug_gene_targets",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Results - Alzheimer's disease

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
# -

al_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, gwas, run_id, SIGNIF_THRESH)

# #### Alzheimer's disease - drug enrichment across annotation results
#
# The gene set analysis for Alzheimer's disease gave the fewest drugs, 8, associated with the disease resulting in just 2 which were be classed as novel. All 8 are shown in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results} table \ref{tab:al_most_frequent_drugs}, with table \ref{tab:al_most_frequent_novel_drugs} listing the 2 novel drugs for which no studies investigating their link to Alzheimer's disease were found.

al_top_drugs_df = rdisp.get_top_drugs_ch5(al_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.al_novel_drugs)
al_top_novel_drugs_df = al_top_drugs_df[al_top_drugs_df['Novel'] == 'yes'].drop(['Novel'], axis=1)

al_top_novel_drugs = al_top_drugs_df[al_top_drugs_df['Novel'] == 'yes'].drop(['Novel'], axis=1)

# Save drug names for publication search script
al_top_drugs_df['Description'].to_csv(data_path + "wrangling/al_drugnames.txt", index=False, header=None)

if formatting == "HTML":
    display(HTML(al_top_novel_drugs_df.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            al_top_novel_drugs_df.head(20).to_latex(
                column_format=r"p{1cm}p{0.8cm}p{5cm}p{8cm}",                                                                                                                               
                caption=r"All novel drugs, ordered by frequency of appearance as significantly "
                        r"associated with Alzheimer's disease across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated, with drug name and drug indications taken from DrugBank.",
                label="tab:al_most_frequent_novel_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### The evidence behind an example novel drug for Alzheimer's disease
#
#  The drug atracurium is the first novel drug appearing in table \ref{tab:al_most_frequent_novel_drugs}, of the two results. Atracurium is only found to be significantly associated with Alzheimer's disease using the DUGGIE drug gene interaction database, which lists six gene targets and in descending Z-score order when annotated using the 'AlleQTLs-mQTL-missing' annotation are CHRNA2, AMH, CHRNA3, CHRNA1, CHRNA4 and CHRNA7, as shown in table \ref{tab:al_top_drug_genes}, with the CHRNA proteins also referred to as acetylcholine receptors. Atracurium is approved as a neuromuscular blocker to relax muscles during mechanical ventilation under general anesthesia or intubation.

# +
annot = 'alltissueshybrid'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
al_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')

# +
al_top_drug_df = rdisp.add_genesets_thesis(al_table[al_table.index == 'M03AC04'], 
                                            zscores_df, 
                                            genesets_df, 
                                            gene_hgnc_df)

al_top_drug_df.drop(['Num genes', 'Q', 'Description', 'ATC code'], axis=1, inplace=True)
al_top_drug_df.rename(columns={'ZSTAT':'Z'}, inplace=True)
al_top_drug_df.set_index('HGNC', inplace=True)
al_top_drug_df['Z'] = al_top_drug_df['Z'].astype('float').round(3)
# -

al_top_drug_df['protein'] = [
    'Cholinergic Receptor Nicotinic Alpha 2 Subunit',
    'Anti-Mullerian Hormone',
    'Cholinergic Receptor Nicotinic Alpha 3 Subunit',
    'Cholinergic Receptor Nicotinic Alpha 1 Subunit',
    'Cholinergic Receptor Nicotinic Alpha 4 Subunit',
    'Cholinergic Receptor Nicotinic Alpha 7 Subunit',
]
al_top_drug_df['evidence sources'] = [
    'STITCH, DGIdb, DrugCentral',
    'DGIdb',
    'GtoPdb',
    'DGIdb, GtoPdb',
    'DGIdb, GtoPdb',
    'DGIdb, GtoPdb',
]

if formatting == "HTML":
    display(HTML(al_top_drug_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            al_top_drug_df.to_latex(
                column_format=r"p{1.5cm}p{1.5cm}p{6cm}p{5cm}",
                caption=r"Genes interacting with the drug atracurium according to the STITCH "
                        r"drug-gene interaction database with gene Z-score associations with Alzheimer's "
                        r"given for the "
                        r"'Brain-missing' annotation and validated sources of evidence used.",
                label="tab:al_top_drug_genes",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### Gene target prevalence across drugs significantly associated with Alzheimer's disease
#
# The annotation giving the largest number of significant drugs (4) from Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results} table \ref{tab:al_cmatrix_drugs}, uses GTEx, psychEncode and eQTLgen eQTLs and mQTL data with position data added in the missing hybrid approach ('AlleQTLs-mQTL-missing' from Table \ref{tab:db_annotations}) with the DUGGIE drug-gene interaction database. Hence DUGGIE is used as a source of drug-target information for Alzheimer's in table \ref{tab:al_drug_gene_targets}, along with the gene-disease association Z-scores and tallies of drugs targeting each gene from the analysis results using this annotation. From these four drugs, two target groups of acetylcholine receptors which have had some interest as potential targets for Alzheimer's treatments<citet data-cite="kihara_alzheimers_2004"><sup>kihara</sup></citet>. Also appearing are genes with a high degree of association with Alzhheimer's - Aurora kinase A (AURKA), polypyrimidine tract binding protein 1 (PTBP1) and anti-Mullerian hormone (AMH).

# +
annot = 'alltissueshybrid'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
al_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')
al_table2 = rdisp.tally_genesets_thesis(al_table, zscores_df, genesets_df, gene_hgnc_df)
# Remove any X-chromosome targets (with no Z-score)
al_table2.dropna(inplace=True)
al_table2 = al_table2.sort_values('count', ascending=False)
al_table2 = al_table2[['count', 'chrom', 'start_pos', 'Z']]
al_table2['Z'] = al_table2['Z'].astype('float').round(3)
# -

if formatting == "HTML":
    display(HTML(al_table2.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            al_table2.head(20).to_latex(
                column_format=r"p{2cm}p{1cm}p{1cm}p{3cm}p{2cm}",
                caption=r"Top 20 drug gene targets for drugs significantly associated with "
                        r"Alzheimer's disease in results from a gene set analysis using the alleQTLs-mQTL-missing "
                        r"annotation and DUGGIE drug-gene interaction database. "
                        r"Each gene target is listed along with the gene chromosome and position, count of "
                        r"significantly associated drugs which target the gene and gene-disease association "
                        r"Z-score resulting from the alleQTLs-mQTL-missing MAGMA SNP annotation stage.",
                label="tab:al_drug_gene_targets",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Results - Parkinson's disease

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
# -

pd_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, gwas, run_id, SIGNIF_THRESH)

# #### Parkinson's disease - drug enrichment across annotation results
#
# After running the gene set analysis pipeline on the 32 different annotation/drug-gene interaction database combinations for the Parkinson's disease GWAS, 50 drugs were found to be significantly associated with the disease as displayed in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results}, table \ref{tab:pd_most_frequent_drugs}. The top 20 drugs from this table found to be novel, out of the 32 novel drugs it contains, are shown in table \ref{tab:pd_most_frequent_novel_drugs}. This list contains gastrointestinal drugs and a large number of anti-infectives, which as a group heavily target cytochrome p450 enzymes, genetic variations in which have been investigated as a potential susceptibility to Parkinson's disease<citet data-cite="ur_rasheed_cytochrome_2017"><sup>rasheed</sup></citet>.

pd_top_drugs_df = rdisp.get_top_drugs_ch5(pd_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.pd_novel_drugs)

pd_top_novel_drugs = pd_top_drugs_df[pd_top_drugs_df['Novel'] == 'yes'].drop(['Novel'], axis=1)

# Save drug names for publication search script
pd_top_drugs_df['Description'].to_csv(data_path + "wrangling/pd_drugnames.txt", index=False, header=None)

if formatting == "HTML":
    display(HTML(pd_top_novel_drugs.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            pd_top_novel_drugs.head(20).to_latex(
                column_format=r"p{1cm}p{0.8cm}R{5cm}R{8cm}",                                                                                                                               
                caption=r"Top 20 novel drugs, ordered by frequency of appearance as significantly "
                        r"associated with Parkinson's disease across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated, with drug name and drug indications taken from DrugBank.",
                label="tab:pd_most_frequent_novel_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### The evidence behind an example novel drug for Parkinson's disease
#
# Pantoprazole is given as one of the novel drugs appearing as significant in the highest number of analyses, using both the STITCH and DUGGIE drug-gene interaction databases with each listing different gene target sets - STITCH targets being a subset of the DUGGIE targets. DUGGIE lists the pantoprazole targets in descending Z-score association with Parkinson's disease as MAPT, DDAH1, CYP2C19, ATP4A, SLC47A2, FASN, CYP2C9, SLC22A2, AHR, VDR, ATP12A, SLC47A1, ATP1A1, NFE2L2, CYP3A4, GCG, with the CYP enzymes and ATP genes also given by STITCH as shown in table \ref{tab:pd_top_drug_genes}.
#
# MAPT is known to be an important susceptibility gene for idiopathic Parkinson's<citet data-cite="coupland_dna_2014"><sup>coupland</sup></citet> and it is listed as a target for pantoprazole in the DsigDB database,which gives a potency measurement (25119) nM but no other mechanism of action, sourced from both PubChem and ChEMBL.

# +
annot = 'alleqtlbrainhybridboth'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
pd_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')

# +
pd_top_drug_df = rdisp.add_genesets_thesis(pd_table[pd_table.index == 'A02BC02'], 
                                            zscores_df, 
                                            genesets_df, 
                                            gene_hgnc_df)

pd_top_drug_df.drop(['Num genes', 'Q', 'Description', 'ATC code'], axis=1, inplace=True)
pd_top_drug_df.rename(columns={'ZSTAT':'Z'}, inplace=True)
pd_top_drug_df.set_index('HGNC', inplace=True)
pd_top_drug_df['Z'] = pd_top_drug_df['Z'].astype('float').round(3)
# -

pd_top_drug_df['protein'] = [
    'Microtubule Associated Protein Tau',
    'Dimethylarginine Dimethylaminohydrolase 1',
    'Cytochrome P450 2C19',
    'ATPase H+/K+ Transporting Subunit Alpha',
    'Solute Carrier Family 47 Member 2',
    'Fatty Acid Synthase',
    'Cytochrome P450 2C9',
    'Solute Carrier Family 22 Member 2',
    'Aryl Hydrocarbon Receptor',
    'Vitamin D Receptor',
    '',
    '',
    '',
    '',
    '',
    '',
]
pd_top_drug_df['evidence sources'] = [
    'DsigDB',
    'DrugBank',
    'STITCH, DGIdb',
    'DrugBank, STITCH, DGIdb',
    'DsigDB, DrugCentral',
    'DrugCentral',
    'STITCH',
    'DsigDB, DrugCentral',
    'DsigDB',
    'DsigDB',
    '',
    '',
    '',
    '',
    '',
    '',
]

if formatting == "HTML":
    display(HTML(pd_top_drug_df.head(10).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            pd_top_drug_df.head(10).to_latex(
                column_format=r"p{1.5cm}p{1.5cm}p{6cm}p{5cm}",
                caption=r"Top 10 genes by Z-score interacting with the drug pantoprazole according "
                        r"to the DUGGIE drug-gene interaction "
                        r"database with gene Z-score associations with Parkinson's given for the "
                        r"'Brain-mQTL-combined' annotation and validated sources of evidence used.",
                label="tab:pd_top_drug_genes",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### Gene target prevalence across drugs significantly associated with Parkinson's disease
#
# Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results} table \ref{tab:pd_cmatrix_drugs} shows that the annotation giving the largest number of significant drugs (40), uses psychEncode eQTLs and mQTL data with position data added in the combined hybrid approach ('Brain-mQTL-combined' from Table \ref{tab:db_annotations}) and the DUGGIE drug-gene interaction database. Hence this combination is used as a source of drug-target and gene association information for Parkinson's and to generate table \ref{tab:pd_drug_gene_targets} of the 20 most commonly occurring drug targets in these results. CYP genes, from the cytochrome P450 superfamily of enzymes fill the top half of the table, many having a moderate association with Parkinson's with a Z-score around 1.

# +
annot = 'alleqtlbrainhybridboth'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
pd_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')
pd_table2 = rdisp.tally_genesets_thesis(pd_table, zscores_df, genesets_df, gene_hgnc_df)
# Remove any X-chromosome targets (with no Z-score)
pd_table2.dropna(inplace=True)
pd_table2 = pd_table2.sort_values('count', ascending=False)
pd_table2 = pd_table2[['count', 'chrom', 'start_pos', 'Z']]
pd_table2['Z'] = pd_table2['Z'].astype('float').round(3)
# -

if formatting == "HTML":
    display(HTML(pd_table2.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            pd_table2.head(20).to_latex(
                column_format=r"p{2cm}p{1cm}p{1cm}p{3cm}p{2cm}",
                caption=r"Top 20 drug gene targets for drugs significantly associated with "
                        r"Parkinson's disease in results from a gene set analysis using the brain-mQTL-combined "
                        r"annotation and DUGGIE drug-gene interaction database. "
                        r"Each gene target is listed along with the gene chromosome and position, count of "
                        r"significantly associated drugs which target the gene and gene-disease association "
                        r"Z-score resulting from the brain-mQTL-combined MAGMA SNP annotation stage.",
                label="tab:pd_drug_gene_targets",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Results - Huntington's disease age of onset

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
# -

hd_drug_sets_df = rdisp.load_drug_sets_ch5(dtis, annots, display_annots, summary_results_dir, gwas, run_id, SIGNIF_THRESH)

# #### Huntington's disease age of onset - drug enrichment across annotation results
#
# The array of gene set analyses for Huntington's disease age of onset produced 24 drugs significantly associated with the disease age of onset. The full list can be seen in Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results}, table \ref{tab:hd_most_frequent_drugs}. No published evidence of studies investigating the drug and Huntington's disease could be found for 18 of these drugs, and these novel repurposing candidates are listed in table \ref{tab:hd_most_frequent_novel_drugs} in order of descending count of the number of analysis results in which they appear. Notable in this list are the DNA mechanism inhibitors daunorubicin, epirubicin and mitoxantrone.

hd_top_drugs_df = rdisp.get_top_drugs_ch5(hd_drug_sets_df, db_indications_df, ATC_LEVELS, rdisp.hd_novel_drugs)

hd_top_novel_drugs = hd_top_drugs_df[hd_top_drugs_df['Novel'] == 'yes'].drop(['Novel'], axis=1)

# Save drug names for publication search script
hd_top_drugs_df['Description'].to_csv(data_path + "wrangling/hd_drugnames.txt", index=False, header=None)

if formatting == "HTML":
    display(HTML(hd_top_novel_drugs.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hd_top_novel_drugs.head(20).to_latex(
                column_format=r"p{1cm}p{0.8cm}p{5cm}p{8cm}",                                                                                                                               
                caption=r"Top 20 novel drugs, ordered by frequency of appearance as significantly "
                        r"associated with Huntington's disease across all annotations analysed. "
                        r"ATC drug codes are listed alongside the count of analyses in which the drug is "
                        r"significantly associated, with drug name and drug indications taken from DrugBank.",
                label="tab:hd_most_frequent_novel_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### The evidence behind an example novel drug for Huntington's disease age of onset
#
# Daunorubicin is found to be significantly associated with Huntington's disease with both DUGGIE and STITCH drug-gene interaction databases. DUGGIE lists 52 gene targets for this drug, including HTT and topoisomerase II enzyme genes, which includes the six targets listed by STITCH - TOP2B, DCXR, TOP2A, ABCG2, ABCB1, TACR1 in ascending association Z-score value. The HTT gene association with Daunorubicin in DUGGIE is obtained from DsigDB, which sites ChEMBL and PubChem as a source and gives a potency of 2239 nM. These sites then give DTC (Drug Target Commons) as the source, which lists the potency measurement as coming from a PubChem bioassay.

# +
annot = 'alleqtlbloodhybrid'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
hd_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')

# +
hd_top_drug_df = rdisp.add_genesets_thesis(hd_table[hd_table.index == 'L01DB02'], 
                                            zscores_df, 
                                            genesets_df, 
                                            gene_hgnc_df)

hd_top_drug_df.drop(['Num genes', 'Q', 'Description', 'ATC code'], axis=1, inplace=True)
hd_top_drug_df.rename(columns={'ZSTAT':'Z'}, inplace=True)
hd_top_drug_df.set_index('HGNC', inplace=True)
hd_top_drug_df = hd_top_drug_df.head(10)
hd_top_drug_df['Z'] = hd_top_drug_df['Z'].astype('float').round(3)
# -

hd_top_drug_df['protein'] = [
    'Huntingtin',
    'Cytochrome P450 1A1',
    'DNA Methyltransferase 3 Alpha',
    'Cytochrome P450 19A1',
    'ATP Binding Cassette Subfamily C Member 1',
    'Tachykinin Receptor 1',
    'Lysine Demethylase 4E',
    'WT1 Transcription Factor',
    'DNA Topoisomerase II Beta',
    'Serine/Threonine Kinase 33',
]
hd_top_drug_df['evidence sources'] = [
    'DsigDB',
    'DGIdb',
    'DGIdb',
    'DsigDB',
    'DGIdb, DrugCentral',
    'TTD',
    'DsigDB',
    'DGIdb',
    'DGIdb',
    'DsigDB',
]

if formatting == "HTML":
    display(HTML(hd_top_drug_df.head(10).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hd_top_drug_df.head(10).to_latex(
                column_format=r"p{1.5cm}p{1.5cm}p{6cm}p{5cm}",
                caption=r"Top 10 genes by Z-score interacting with the drug daunorubicin according "
                        r"to the DUGGIE drug-gene interaction "
                        r"database with gene Z-score associations with Huntington's given for the "
                        r"'Blood-missing' annotation and validated sources of evidence used.",
                label="tab:hd_top_drug_genes",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# #### Gene target prevalence across drugs significantly associated with Huntington's disease age of onset
#
# From Appendix \ref{appendix-full-neuropsychiatric-and-neurodegenerative-disease-results}, table \ref{tab:hd_cmatrix_drugs}, the annotation giving the largest number of significant drugs (12), using GTEx(Blood) eQTLgen eQTL data with position data added in the missing hybrid approach ('Blood' from Table \ref{tab:db_annotations}) and the DUGGIE drug-gene interaction database. This is used as a source of drug-target and gene association information for Huntington's age of onset to produce table \ref{tab:hd_drug_gene_targets} of the top gene targets of drugs significantly associated Huntington's age of onset. The top two entries are both topoisomerase type II genes which are involved in untangling DNA strands during replication.

# \newpage

# +
annot = 'alleqtlblood'
file_postfix = 'duggie_qvals.tsv'

zscores_df = rdisp.get_annot_results(gwas, annot, magma_results_dir, run_id)
hd_table = count_drug_disease_genes(file_postfix, gwas, run_id, annot, 'duggie')
hd_table2 = rdisp.tally_genesets_thesis(hd_table, zscores_df, genesets_df, gene_hgnc_df)
# Remove any X-chromosome targets (with no Z-score)
hd_table2.dropna(inplace=True)
hd_table2 = hd_table2.sort_values('count', ascending=False)
hd_table2 = hd_table2[['count', 'chrom', 'start_pos', 'Z']]
hd_table2['Z'] = hd_table2['Z'].astype('float').round(3)
# -

hd_table2 = hd_table2.drop([np.NaN])

if formatting == "HTML":
    display(HTML(hd_table2.head(20).to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hd_table2.head(20).to_latex(
                column_format=r"p{2cm}p{1cm}p{1cm}p{3cm}p{2cm}",
                caption=r"Top 20 drug gene targets for drugs significantly associated with "
                        r"Huntington's disease in results from a gene set analysis using the blood "
                        r"annotation and DUGGIE drug-gene interaction database. "
                        r"Each gene target is listed along with the gene chromosome and position, count of "
                        r"significantly associated drugs which target the gene and gene-disease association "
                        r"Z-score resulting from the blood MAGMA SNP annotation stage.",
                label="tab:hd_drug_gene_targets",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Overlap of significantly associated drugs across all diseases
#
# Using the drugs found to be significantly associated with each disease pooled from the 32 gene set analysis pipeline executions, the overlap of these drug sets is shown in the UpSet plot in figure \ref{fig:drug_overlap_upset}. This plot shows that there is a large amount of drugs shared between many of the drug sets significantly associated with each disease, particularly between schizophrenia and bipolar disorder (43 unique shared drugs) and to a lesser extent depression, which shares 12 potentially repurposable drugs jointly with schizophrenia and bipolar, as well as 2 drugs with just schizophrenia and 8 with just bipolar disorder. This is consistent with what is understood of the shared genetic architecture between these disorders<citet data-cite="grotzinger_shared_2021"><sup>grotzinger</sup></citet>.

# +
drug_sets_df = pd.DataFrame()
labels = ('SCZ', 'ALZ', 'PD', 'MDD', 'BPD', 'HD')
drugs = (sz_drug_sets_df, al_drug_sets_df, pd_drug_sets_df, mdd_drug_sets_df, bp_drug_sets_df, hd_drug_sets_df)

for index in range(0,6):
        drug_set = drugs[index].index
        # Add series as column to get a matrix of boolean columns
        drug_sets_df[labels[index]] = False
        for drug in drug_set:
            if drug not in drug_sets_df.index:
                drug_sets_df = drug_sets_df.append(pd.Series([False], name=drug))
            drug_sets_df.at[drug, labels[index]] = True
        
drug_sets_df = drug_sets_df.fillna(False)
drug_sets_df = drug_sets_df.drop(0, axis=1)
# -

bold_caption = (
    r"UpSet plot showing the overlap between diseases of significantly associated drugs found for each disease."
)
normal_caption = (
    r"Vertical bars show the details for each disease, with the histogram showing the significant"
    r" drug tallies per disease. The horizontal histogram shows the size of intersection between drug sets"
    r" indicated horizontally in the dot matrix. Single dots represent the tally of drugs unique to"
    r" the named disease."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(drug_sets_df), 
                 subset_size='count', 
                 fig=usp_figure,
                 show_counts=True, 
                 orientation='vertical')

rdisp.end_caption(formatting, r"fig:drug_overlap_upset")

# +
#scz_bpd_list = drug_sets_df[drug_sets_df['SCZ'] & drug_sets_df['BPD']].index

#scz_bpd_df = db_indications_df.loc[scz_bpd_list]

## or look at bp_top_novel_drugs_df
#pd.merge(sz_top_novel_drugs, scz_bpd_df, left_index=True, right_index=True, how="inner") 
# -

# Of the 55 drugs (43 shared between schizophrenia and bipolar disorder plus the 12 common to the drug sets associated with schizophrenia, bipolar disorder and major depressive disorder) from figure \ref{fig:drug_overlap_upset}) significantly associated with both schizophrenia and bipolar disorder, one has an indication for schizophrenia, pipotiazine (ATC code N05AC04), a novel repurposing candidate for bipolar disorder and another, minaprine (ATC code N06AX07), with an indication of depressions and a novel repurposing candidate for schizophrenia. Both of these drugs could be prioritised as novel repurposing candidates for the associated diseases, on account of the shared genetic architecture between the diseases.

# ### Overlap of target genes of significantly associated drugs across all diseases
#
# Similar to the analysis of significantly associated drugs, the drug target gene sets driving the most successful annotations in analysing each disease and their overlaps can be seen in the UpSet plot in figure \ref{fig:gene_overlap_upset}. Here, the gene target sets used to produce the tables of top 20 drug gene targets above (tables \ref{tab:sz_drug_gene_targets}, \ref{tab:mdd_drug_gene_targets}, \ref{tab:bp_drug_gene_targets}, \ref{tab:al_drug_gene_targets}, \ref{tab:pd_drug_gene_targets} and \ref{tab:hd_drug_gene_targets}) are compared against each other. The results are not unsurprisingly complementary to the drug overlaps, with schizophrenia and bipolar disorder having the highest degree of overlap, sharing 147 gene targets uniquely between them, and both sharing a further 58 gene targets uniquely with major depressive disorder. Of note also is the larger unique overlap of genes (41) between schizophrenia, bipolar and depression with Parkinson's disease drug targets, which may also allude to a shared genetic architecture, with evidence that their comorbidity has been previously investigated<citet data-cite="li_overlapping_2021"><sup>overlap</sup></citet>. This confirms the large amount of shared genetic architecture between these diseases.

# +
target_sets_df = pd.DataFrame()
labels = ('SCZ', 'ALZ', 'PD', 'MDD', 'BPD', 'HD')
targets = (sz_table2, al_table2, pd_table2, mdd_table2, bp_table2, hd_table2)

for index in range(0,6):
        target_set = targets[index].index
        # Add series as column to get a matrix of boolean columns
        target_sets_df[labels[index]] = False
        for target in target_set:
            if target not in target_sets_df.index:
                target_sets_df = target_sets_df.append(pd.Series([False], name=target))
            target_sets_df.at[target, labels[index]] = True
        
target_sets_df = target_sets_df.fillna(False)
target_sets_df = target_sets_df.drop(0, axis=1)
# -

bold_caption = (
    r"UpSet plot showing the overlap between gene targets of significantly associated drugs found for each disease."
)
normal_caption = (
    r"Vertical bars show the drug target set details per disease, with histogram showing the significant"
    r" target tallies per annotation. The horizontal histogram shows the size of overlap between the target sets indicated"
    " horizontally in the dot matrix. Single dots represent the tally of targets unique to the named disease."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(target_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical')
usp_figure.set_size_inches(6, 9)

rdisp.end_caption(formatting, r"fig:gene_overlap_upset")

# ## Discussion
#
# This chapter presented the results of testing six sets of neuropsychiatric and neurodegenerative disease GWAS summary statistics for enrichment of association signal, using a selection of different SNP to gene annotations and two curated databases of drug-gene interactions - one of which, DUGGIE, was curated as part of this project and the details of which can be found in chapter \ref{curation-of-a-drug-target-interaction-database}. The six diseases studied were three neuropsychiatric diseases - schizophrenia, major depressive disorder and bipolar disorder, and three neurodegenerative conditions - Alzheimer's disease, Parkinson's disease and Huntington's age of onset. 
#
# All six analyses yielded sets of both novel and known repurposing drug candidates, with some existing treatment drugs also included in the results. The most frequently associated novel drugs for each disease were listed along with their current treatment indications, and these represent the best candidates for potential drug repurposing. The drug-gene interaction evidence behind the drug ranked top in each of these lists was investigated to validate that the evidence was sound and not subject to errors in the dataset curation or analysis pipeline. The top ranked non-CACNA1C targeting drug repurposing candidate for schizophrenia was also considered, as well as the top non-CACNA1C targeting drug repurposing candidate for bipolar disorder, as it shared its top ranked drug with that from the schizophrenia analysis (magnesium sulfate) - all the drug-gene interaction evidence examined was found to be valid when doing this.
#
# For schizophrenia, 26 out of the 109 significantly associated drugs were placed in the final list of potential repurposing candidates following the literature search, with a prevalence of drugs approved for hypertension and related indications that target CACNA1C and associated calcium channel genes, targets known to be shared between hypertension and schizophrenia<citet data-cite="harrison_voltage-gated_nodate"><sup>harrison</sup></citet>, including magnesium sulfate whose supporting drug target evidence was checked. Cisplatin was also examined as a non-CACNA1C targeting drug, which was identified as having its association with schizophrenia driven by a histone gene target. 
#
# Major depressive disorder was significantly associated with 49 drugs via the gene set analyses and following the literature search for studies investigating drug-disease links, 20 of these drugs remained as novel potential drug repurposing candidates. Investigating the targets of significantly associated drugs from the annotation and drug-gene interaction database combination giving the largest number of significant drugs, serotonin and dopamine receptor genes featured extensively in the list of most frequently occurring targets. These receptors are well known to have a role in depression<citet data-cite="moncrieff_serotonin_2022, clausius_relevance_2009"><sup>moncrieff, clausius</sup></citet>.
#
# The set of gene set analyses of bipolar disorder gave 116 drugs associated with the indication, of which 49 were left after discarding drugs found to have been studied previously as repurposing candidates. As was the case for schizophrenia, which is known to have close genetic overlap with bipolar disorder, CACNA1C and calcium channel targets were the type of gene most frequently found among the targets of significantly associated drugs. The drug target evidence supporting a non-calcium channel targeting drug, ethchlorvynol, was examined. This drug targets GABA receptor genes, which have been implicated in the disease<citet data-cite="breuer_independent_2011"><sup>breuer</sup></citet>.
#
# Only two novel drug repurposing candidates were obtained from analysing Alzheimer's disease, from 8 drugs found to be significantly associated. From this small number of drugs and genes, an examination of gene targets across drugs as well as the top ranked novel drug investigated, atracurium, gave as frequent targets acetylcholine receptors which have been previously investigated as a potential gene target for Alzheimer's disease<citet data-cite="kihara_alzheimers_2004"><sup>kihara</sup></citet>.
#
# The gene set analyses executed using the Parkinson's disease GWAS identified 50 drugs as significantly associated with the disease, of which 32 remained after discarding drugs found to have previous studies investigating their interaction with the disease in the literature. Cytochrome P450 enzymes are the most frequent type of gene targeted by drugs that were significantly associated with the disease by the most successful annotation method, despite individually having moderate association Z-scores with the disease. These enzymes also appear as a target of the sampled novel candidate drug, pantoprazole, whose association is otherwise driven by the MAPT gene - a known susceptibility gene for Parkinson's disease. This looks a complex but interesting relationship to be further pursued as studies have been performed investigating the action of foreign substances within the body, known as xenobiotics, which have been linked with Parkinson's disease<citet data-cite="steventon_review_2001,bjorklund_role_2020"><sup>steventon,bjorklund</sup></citet> and which cytochrome P450 enzymes play a part in metabolising<citet data-cite="esteves_central_2021"><sup>esteves</sup></citet>.
#
# The Huntington's disease age of onset gene set analysis found 24 drugs to be significantly associated with the feature, of which 18 remained following a literature search for previous studies, including several drugs affecting DNA repair mechanisms, which are involved in pathways that are significantly enriched in GWAS signal for HD onset and progression<citet data-cite="maiuri_dna_2019"><sup>maiuri</sup></citet>. Similarly, the top two targets taken from all drugs found by the annotation method significantly associating the most drugs are both topoisomerase type II genes which are involved in DNA replication and repair processes<citet data-cite="pommier_human_2022"><sup>pommier</sup></citet>. The connection between DNA repair mechanisms and Huntington's disease age of onset is well known<citet data-cite="bettencourt_dna_2016,goold_fan1_2019"><sup>bettencourt,goold</sup></citet>, making these drugs and targets of interest for further work. The sample novel repurposing candidate chosen, daunorubicin, also had a top target of the huntingtin gene itself.
#
# Finally, a further validation against previously known gene overlaps between the neuropsychiatric diseases (schizophrenia, bipolar disorder and major depressive disorder) was made by comparing the gene targets from the results of the annotation and drug-target interaction databased giving the largest significantly associated set of drugs for each disease. This confirmed the overlap, with the extra observation of a noticeable overlap of gene targets between the three neuropsychiatric diseases and Parkinson's disease. The significantly associated drugs for each disease also exhibited the same magnitude of overlap as the target genes, which also allowed two drug repurposing candidates to be highlighted as potentially based on this shared genetic architecture, pipotiazine for bipolar disorder and minaprine for schizophrenia, and hence these could be prioritised in further repurposing investigations.
