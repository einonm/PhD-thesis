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


# # Discussion and conclusions
#
# ## Overview
#
# The aim of this work was to identify novel drug repurposing opportunities for neuropsychiatric and neurodegenerative diseases using a competitive gene-set analysis pipeline. This pipeline enables the user to statistically associate drugs with a disease by testing their sets of target genes for enrichment of association in a GWAS and, further, rank the drugs based on the relative enrichment of each drug's target gene set against the target gene sets of all other approved drugs involved in the analysis. To ensure this analysis was as comprehensive as possible, two other activities were carried out to choose analysis pipeline input data sources providing the greatest number of drugs significantly associated with the disease:
#
# 1. The curation and validation of the largest approved drug-target interaction database possible using available published drug-target databases.
# 2. The mapping of SNP-disease associations obtained from GWAS summary results to gene-disease associations using a wide range of functional QTL based annotations, measuring their ability to identify drug repurposing opportunities.
#
# The resulting choices of input data sources were used to execute a battery of gene-set analyses for each of the neuropsychiatric and neurodegenerative diseases of interest and rank drugs based on the number of annotation and drug-gene interaction database combinations in which the drug's target genes appear as significantly enriched for GWAS signal. Any ranked drugs thus identified were classified as novel or not novel depending on the discovery of a previous drug approval or existing studies investigating the effect of the drug on the disease.
#
# Chapter \ref{curation-of-a-drug-target-interaction-database}, gave an overview of the curation of a new drug-target gene database, which involved choosing a definition of a drug's 'gene target', a literature search for available drug-target information, selection of eligible data sources taking account of database cross-inclusion, collection and amalgamation of the selected drug-target data, and finally a comparative analysis of the resulting database titled 'DUGGIE' (DrUG-Gene IntEractions) with its constituent data sources including the current largest drug target database, STITCH. The new database, DUGGIE, was found to have 76% more drugs and 17% more targets whilst also capturing 2.2 times the number of drug-target interactions than STITCH.
#
# Chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia} then explored any relative advantage in associating treatment drugs with disease that DUGGIE may have over STITCH. It did this by executing identical analysis pipelines using a positional SNP-gene annotation, save for the different drug target interaction databases, on two well-understood diseases each with an array of approved drug treatments - hypertension and hypercholesterolemia. The abilities of each drug target interaction database to identify drugs approved for the disease (treatment drugs) and discard non-approved (non-treatment) drugs were compared using several direct measurements, principally the number of correctly identified approved drugs, and several statistics - Fisher p-values, Mann-Whitney and ROC curves - measuring the ability of the analyses to differentiate between the two groups of treatment and non-treatment drugs. Further evidence was used of drug side-effects to ascertain if the drug has a negative effect on the disease, such as increasing blood pressure with hypertension, which could account for a non-treatment drug having a significant association with a disease.
#
# This analysis revealed that both DUGGIE and STITCH were adept at identifying treatment drugs and differentiating between treatment and non-treatment drug groups to a statistically significant level. However, despite DUGGIE finding a greater number of drugs to be significantly associated with each disease, STITCH still performed better at identifying a greater number of treatment drugs as significant. Interpreting this finding, it could be because STITCH may contain a greater proportion of well understood gene targets which align well with targets of comprehensively studied approved drugs, and that the extra targets and drug-target interactions provided by DUGGIE are more novel, and hence may be of greater benefit in finding novel repurposing opportunities - further work would be required to confirm this hypothesis. As both DUGGIE and STITCH were effective at associating drugs with disease, with STITCH focussing on approved drugs and DUGGIE being able to identify a larger set of disease-associated drugs, a decision was made to include both in the final analysis of neuropsychiatric and neurodegenerative diseases as well as in the investigation into the use of functional annotations to improve the pipeline outcomes.
#
# Chapter \ref{using-functional-qtl-gene-annotation-to-better-identify-genes-associated-with-disease} covered a study using functional gene annotation to map GWAS-derived SNP-disease association p-values to gene-disease association p-values. Due to the lower gene coverage of functional QTL measurements, several hybrid tissue-specific functional/positional annotation schemes were analysed along with purely functional and positional annotations to ascertain if an advantage in detecting treatment drugs could be discovered. This analysis was performed using the same hypertension and hypercholesterolemia GWAS summary results used in the previous chapter, but differing in the annotation schemes used. Along with the proximity annotation used previously, fifteen other annotations were considered using blood, brain and all-tissues functional QTL data combined with positional annotation data in several hybrid schemes. Each of these schemes was analysed using the gene-set analysis pipeline and both DUGGIE and STITCH, resulting in 32 analysis executions per disease. Again, Fisher p-values and Mann-Whitney test statistics were used to compare the ability of the different annotation schemes to correctly classify treatment and non-treatment drugs, alongside the metric of number of treatment drugs identified by each. This annotation comparison was used to test four hypotheses stated as:
#
# 1. Functionally annotating genes increases the number of significant drugs and treatment drugs in analysis pipeline results
# 2. Augmenting functional annotations with positional data increases the number of significant drugs and treatment drugs in analysis pipeline results
# 3. Using QTL data obtained from larger numbers of tissue samples increases the number of significant drugs and treatment drugs in analysis pipeline results
# 4. Using the specific brain tissue methylation-QTL data alongside expression QTL data increases the number of significant drugs and treatment drugs in analysis pipeline results.
#
# The results of this annotation study firstly showed that overall, the gene-set analysis pipeline was capable of significantly differentiating between the treatment drug and non-treatment drug groups under most types of annotation. The results also confirmed the conclusions from chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}, that DUGGIE overall identifies a greater number of drugs as significantly associated with a disease than STITCH, but STITCH is more adept at finding treatment drugs within those that it identifies as significantly associated with a disease. Conversely, this could be interpreted as DUGGIE identifying a larger number of potential novel repurposing candidates which would be present in the 'significantly associated, non-treatment' drug category.
#
# Regarding the first two hypotheses given above involving annotations, some evidence was found to support the use of functional annotations, and there was stronger evidence for a hybrid functional/positional annotation approach over pure functional annotation, but no discernible pattern was evident to enable a choice of which functional or hybrid approach would be better in any particular situation. Additionally, no evidence was found to support the hypothesis that functional QTL annotations involving greater numbers of tissue samples result in an increased number of significant or treatment drugs, and it was also inconclusive whether the addition of methylation-QTL data increased the amount of significant or treatment drugs identified. 
#
# Hence due to the nature of the outcome sought, i.e. the identification of novel repurposing candidates, and to avoid discarding potentially useful evidence for associating a particular drug with a disease, it was decided that evidence from analyses across the whole range of annotation methods investigated could be usefully employed in compiling a list of repurposing candidates for a disease. Following analysis of the ability of the amalgamated list of significantly identified drugs to prioritise treatment drugs over non-treatment drugs, it was also decided to rank this list by the number of analyses in which the drug was found to be significant in subsequent investigations. This also suggests that the GWAS enrichment is somewhat robust to the choice of annotation and drug-gene interaction database combination.
#
# The final results chapter, chapter \ref{application-of-the-analysis-pipeline-to-neuropsychiatric-and-neurodegenerative-diseases}, used multiple runs of the gene-set analysis pipeline over the battery of annotation schemes with both DUGGIE and STITCH drug-gene target interaction databases to generate a list of drugs found to be significantly associated with a disease. Drugs were ranked by the number of annotation and drug-gene interaction database combinations in which the drug was found to be associated. Six diseases were investigated with this method - schizophrenia, major depressive disorder, bipolar disorder, Alzheimer's disease, Parkinson's disease and the age of onset of Huntington's disease. From the resulting ranked lists of significantly associated drugs for each disease, an initial decision as to which drugs may be novel repurposing candidates was made firstly by removing drugs that have been approved to treat the disease and secondly based on a literature search for the appearance of studies investigating the effectiveness of the drug for that disease. If no evidence of studies investigating a drug with a disease was found, then the drug remained a potential novel repurposing candidate. Following this, the targets of the top ranked drug selected from each resulting list of a disease's novel repurposing candidates were examined and the supporting evidence validated. Also, for each of the six diseases, a list of genes targeted by drugs found to be significantly associated with the disease was compiled, using the annotation and drug-gene interaction database combination that was successful in finding the most significant drugs. This paves the way for further study such as testing for biological pathways among these genes. Lastly, the drug and target overlaps between the neuropsychiatric diseases were investigated to highlight any drugs that are identified as repurposing candidates which are common between diseases that share an underlying genetic architecture, and have been successful in treating at least one of these diseases, with the increased possibility that the drugs chances of success in treating another disease sharing a genetic architecture.
#
# Multiple novel drug repurposing candidates, ranging between 2 and 49 for each disease, were identified in this manner for all six diseases, following the removal of drugs found to be significantly associated with the disease but having evidence in the literature of studies into their interaction with the disease. From examining the supporting drug target experimental evidence of a small sample of novel drugs, all the evidence was found to be present and credible, further validating the pipeline. The gene targets of significantly associated drugs produced a series of gene sets also useful for further study - with some results showing strong alignment in the genes targeted by large numbers of significantly associated drugs, for example calcium channel genes from schizophrenia and bipolar disorder, and cytochrome P450 enzymes for Parkinson's disease. Finally, the drug and target set overlap between each disease was found to concur with current knowledge of the shared genetic architecture of schizophrenia, bipolar disorder and major depressive disorder. 

# ## Updates to this study in context
#
# Following the literature search in section \ref{this-study-in-context} exploring the wider scientific context of the work outlined in this thesis at it's commencement in 2017, and the conclusion that the approach is novel, other research may have been published to contradict or challenge this view. To explore this possibility the same literature search was undertaken again, searching the PubMed<citet data-cite="noauthor_pubmed_nodate"><sup>pubmed</sup></citet> database using the boolean search string '((drug) OR (drugs)) AND (("gene set") OR ("gene sets"))', constraining the results to studies involving the human species and whose publication text is marked as being freely available, but covering the years 2018 to 2023. Following manual curation of the abstracts and titles of the 1,308 results returned, 551 did not involve gene set enrichment and were discarded. A further 7 publications focussed on alternatives to gene set enrichment methods mentioning drugs in passing, 4 were retracted and 740 involved mostly either pathway gene set enrichment or single-sample gene set enrichment analysis to ascertain immune infiltration levels, mainly concerning cancer research. The remaining 8 publications were of interest, involving methods or the application of gene set enrichment techniques to drug target gene sets.
#
# These 6 publications involved 4 unique analysis pipelines, with 2 named pipelines appearing in two publications each. Examining the 2 single-publication pipelines in comparison to the methods employed in this study:
#
# * *Calcium channel blockers as drug repurposing candidates for gestational diabetes: Mining large scale genomic and electronic health records data to repurpose medications* (Goldstein et al.)<citet data-cite="goldstein_calcium_2018"><sup>goldstein</sup></citet> used a bespoke Python pipeline for gene set enrichment, taking drug-gene interaction data from the Psychoactive Drug Screening Programme (PDSP, a contributing source to STITCH) and Drugbank, with a further literature search to identify targets for remaining drugs with no targets. The resulting 129 drugs were grouped into 65 'mechanistic' classes used as gene sets in the analysis. Gene-level disease associations were obtained by taking the strongest SNP per gene using Fisher's method (no clear annotation scheme mentioned), corrected using permutations. 11 of the drug classes were found to have an association with gestational diabetes.
#     
#     
# * *Host - virus - drug interactions as determinants of COVID-19's phenotypes: A data-driven hypothesis* (Vavougios)<citet data-cite="vavougios_host_2020"><sup>vavougios</sup></citet> performs a brief drug - gene set enrichment analysis using the 'Enrichr' web service, performing a simple Fisher's exact test, involving the Gene Expression Omnibus and DsigDB (a contributing source to DUGGIE) as drug-gene interaction data sources, with differentially expressed genes extracted from pathway lists and a second case-control COVID-19 expression dataset used to provide trait-associated gene p-values. Further details of the method are not clear. Several substances were found to have a significant association.
#
# And for the 2 named pipelines appearing in two publications each:
#     
# * dpGSEA - 'Drug perturbation Gene Set Enrichment Analysis' introduced in *Drug perturbation gene set enrichment analysis (dpGSEA): a new transcriptomic drug screening approach* (Fang et al.)<citet data-cite="fang_drug_2021"><sup>fang</sup></citet> is based on the premise that a drug-gene transcription profile negatively correlated with a disease-gene profile can rescue a disease state. To highlight the most promising drugs with this possibility, the method creates a matrix of 10, 20 or 50 top drug targets annotated by cell line and direction of gene regulation with which to replace a drug's gene-set in a gene set enrichment analysis, and producing an extra statistic based on the level of drug-gene interaction and regulation, the 'Target Compatibility Score' (TCS), alongside an enrichment score, to capture the statistical significance of this relationship. Favourable comparisons are made with other methods, but no absolute comparison with known drug/disease outcomes are made. Another recent publication using dpGSEA *Screening Potential Drugs for the Development of NAFLD Based on Drug Perturbation Gene Set* (Gao et al.)<citet data-cite="gao_screening_2022"><sup>gao</sup></citet> uses dpGSEA as a secondary analysis of Nonalcoholic fatty liver disease, but only highlights seven drugs which have significant enrichment score p-values, with no discussion of TCS.
#     
#     
# * CGSEA - 'Chemical-related Gene Set Enrichment Analysis' is a method first described in *CGSEA: A Flexible Tool for Evaluating the Associations of Chemicals with Complex Diseases* (Cheng et al.)<citet data-cite="cheng_cgsea_2020"><sup>cheng</sup></citet>, using the Comparative Toxico-genomic Database (CTD, a contributing source to STITCH) as a source of drug-gene interaction data, where the gene-disease association Z-scores are obtained from the output of a Transcriptome Wide Association Study (TWAS), calculated using the FUSION software tool. A TWAS analysis involves combining GWAS summary statistics with tissue-specific gene expression weights to calculate association statistics between gene expression levels and a trait or disease whic hcan be used in a gene set enrichment analysis - this publication used brain RNA-seq and whole-blood RNA expression weights alongside summary GWAS statistics to study autism spectrum disorder (ASD) and attention deficit hyperactivity disorder (ADHD), and cervical squamous cell carcinoma RNA-seq expression weights along with whole-blood RNA expression weights reference data to study cervical cancer. For each disease, several chemicals were found to be significantly associated. The second CGSEA publication, *Identifying Novel Drug Targets for Epilepsy Through a Brain Transcriptome‑Wide Association Study and Protein‑Wide Association Study with Chemical‑Gene‑Interaction Analysis* (Lu et al.)<citet data-cite="lu_identifying_2023"><sup>lu</sup></citet> applied the same method to epilepsy, using expression weights from multiple brain regions. Multiple associations for each brain region studied were reported.
#
# In addition, a further PubMed search for '((drug or drugs) AND magma)' was made, between the years 2015 (the initial MAGMA publication) and 2023, producing 12 results involving the MAGMA gene set analysis tool after removing publications using the term 'magma' to describe other concepts. None of these publications were found to use the MAGMA tool with drug target gene sets.
#
# Of the two named pipelines CGSEA and dpGSEA, both are of interest as they use relevant drug-gene interaction data and attempt to augment a gene set enrichment analysis with functional information, including the direction of gene expression changes, and a comparison with the pipeline outlined in this work would be illuminating and a suggestion for future work. However, none use QTL data explicitly or employ such a comprehensive drug-gene interaction dataset as DUGGIE, implying that this study still employs a novel approach at the time of publication.

# ## Suggestions for future work
#
# This study identified novel drug repurposing candidates for several neuropsychiatric and neurodegenerative diseases, work which does not provide conclusive evidence to directly enable improved clinical outcomes, but is intended as evidence for further investigations. These investigations would continue the search for evidence for and against the suitability of the identified novel drugs as a therapeutic, or even as harmful to the disease - information which would still be of use and is not captured by the methods employed in this study. Such investigations may involve further studies covering the effect of the specific gene targets assigned to a drug by DUGGIE, such as the direction of gene expression changes the drugs elicit, utilising perturbagen datasets such as the Connectivity Map (CMAP)<citet data-cite="subramanian_next_2017"><sup>subramanian</sup></citet>, or analyses predicting gene expression changes<citet data-cite="so_analysis_2017"><sup>so</sup></citet>, which could indicate the likelihood of the drug response being therapeutic. These, along with investigations of drug side effects, further target validation and drug delivery can aid in further discarding unsuccessful drugs before the clinical trial stage. However, the multi-target (polypharmacological) nature of the drug evidence presented here for the polygenic diseases studied means that such work is much more complex than the single-target validation of traditional drug discovery.
#
# Some specific questions that could be addressed that were not answered by this study are:
#
# * Does the STITCH drug-gene interaction database contain a greater proportion of well understood gene targets that align better with the targets of treatment drugs than the extra targets present in DUGGIE? Thus explaining the bias in identifying treatment drugs when using STITCH in the gene-set analyses. This could be investigated with the assumption that the target gene set of a treatment drug is enriched for genes highly associated with the disease, then testing the hypothesis that the enrichment of these highly disease-associated targets in STITCH is greater than the enrichment in extra targets added by DUGGIE, across all drugs. A further, complimentary hypothesis could also be tested that the extra DUGGIE targets allow for a greater enrichment of disease-associated targets in non-treatment drugs, thereby increasing the chance of them being significantly associated with the disease.
#
# * How does the choice of tissue types used to generate QTL data for a gene-set analysis affect the outcome for a particular disease, given that the disease is assumed to involve only a sub-set of tissue types? The broad sets of blood, brain and all tissue types used for this project did not reveal any discernible pattern, but the choice of tissue type did vary the outcome, considerably in some cases. Further investigation could involve using tissue types that have been identified from single-cell expression studies<citet data-cite="jagadeesh_identifying_2022"><sup>jagadeesh</sup></citet> as being disease-critical, and using QTL data from this tissue subset to run further gene-set analyses to ascertain if this strategy improves results. 
#
# * What are the implications of the enzyme-targeting drugs that feature prominently in the identified Parkinson's disease analyses (Table \ref{tab:pd_drug_gene_targets})? Is this a potential therapeutic avenue, or just a weak signal highlighted by the lack of more relevant targets? Previously, there have been studies linking the role of cytochrome P450 enzymes to Parkinsons, particularly CYP2D6<citet data-cite="mann_neuroprotective_2012,viaggi_cytochrome_2006,ur_rasheed_cytochrome_2017, riedl_p450_1998"><sup>mann,viaggi,ur,reidl</sup></citet>, but no firm evidence was found for Parkinson's therapeutics targeting this family of proteins - so investigating the possibility of a potential therapeutic targeting these enzymes could be fruitful.
#
# Other questions of interest that could motivate future work are:
#
# * We have a good estimate of the total number of coding genes in the human genome, around 19,000<citet data-cite="piovesan_human_2019"><sup>piovesan</sup></citet>, but a comparable estimate of, for example, the number of gene interaction pathways or the proportion of expression QTLs that could be found was not encountered whilst undertaking this study, which would have been of great use in gauging how complete the input datasets used were. Knowing how much of the input feature space is missing could give credence to the results if the missing fraction is small, and conversely would increase skepticism in them if the proportion of missing information is high. Similarly with the use pathway and protein-protein interaction data, an estimate of how complete the coverage of the actual pathway and interaction space is would be of interest.
#
# * The ATC drug classification system is not particularly relevant to genetics, as it is compiled solely on the basis of the substance and the therapeutic effect of the drug. A similar scheme that takes into account drug gene targets and mechanisms of action would be useful in undertaking future studies similar to this one.

# ## Concluding remarks
#
# The novel repurposing candidate results obtained in this study are very encouraging, with many novel candidate drugs identified for the six diseases studied for which no previous studies investigating the drug-disease links were found. These drugs were identified using a pipeline tested for effectiveness on two well-studied diseases for which a large number of treatment drugs exist, utilising two novel input data sources engineered for this project - firstly DUGGIE, a drug-target interaction database based on experimentally obtained evidence which is larger than any other such database, and secondly a collection of functional QTL based SNP to gene annotation methods. 
#
# However, there are some drawbacks to the pipeline, most notably that the direction of effect of the drug on the disease is missing, so the drug response may worsen the disease as well as potentially alleviate symptoms. Having a larger drug-gene interaction database may also cause a drug with too many ineffectual targets, i.e. ones having a low association with a disease, to subsequently have a lower disease association - whereas in reality, only the effective targets are relevant to the effectiveness of the drug.
#
# For some of the diseases studied a very large number of drugs, over one hundred for schizophrenia and bipolar disorder, were significantly associated by the analysis. Any further improvements made to the analysis may result in even higher numbers of significant associations, suggesting that more stringent quality controls and tighter statistical thresholds may be beneficial in prioritising a manageable set of candidate repurposing opportunities. Even so, as the analysis only provides a relative difference between drugs and not an absolute measure of effect on the disease, setting a reasonable set of controls and thresholds suitable for any disease would be a challenging task.
#
# Furthermore, the pipeline results are not based on protein-protein interaction network or pathway relationships; only direct drug target gene to disease gene interactions are considered. Incorporating protein-protein interaction and pathway information may make the analysis much more complex, for example involving the use of graph databases (such as Neo4j<citet data-cite="noauthor_neo4j_nodate"><sup>neo4j</sup></citet>) to analyse and control for the number of graph edges between drugs target genes and disease relevant genes. Such data may greatly increase the confidence in drugs that are identified as potential repurposing candidates by highlighting interaction pathways that are known to be involved in the disease or protein-protein interactions that confirm the drug is unsuitable due to unwanted side effects. An alternative argument against including the second order / third order / n-order interactions is that this may duplicate information already captured by the GWAS p-value associations of the n-order interaction genes with the disease, adding undesirable statistical noise.
#
# In considering any further work investigating the suitability and effectiveness of a drug repurposing candidate for a disease, even when the gene set analysis pipeline identifies a promising highly significant association there are still many hurdles and risks to navigate before the repurposing can be deemed a success. These hurdles could be pharmacological, such as identifying which tissue types are of interest, proving that the drug can be transported in the body to these tissue types and effective dosage levels, or regulatory and legal, for example ensuring adherence to patents and intellectual property laws and any relevant ethics, permissions and certifications are obtained<citet data-cite="begley_drug_2021"><sup>begley</sup></citet>. This study did not set out to cover these factors due to the large burden of work and diverse expertise required. 
#
# It has been shown in this thesis that gene set analysis can be useful for identifying potential repurposing candidates beyond looking at single targets for complex polygenic diseases, as the method considers the whole set of known genes that can be targeted in assessing the relative potential effectiveness of a drug. It may be that this polypharmacological approach allows drugs that work on many facets of a disease aetiology at once to be chosen, as the analysis will favour drugs that interact with many disease relevant genes, not just one. It further follows logically that one drug alone may not be a particularly effective therapeutic for a disease, but several used in combination may be. In this way, the pipeline could also in future be easily extended to use groups of genes targeted by more complex combinations of drugs, to evaluate their combined effectiveness over other possible sets of drug targets.
