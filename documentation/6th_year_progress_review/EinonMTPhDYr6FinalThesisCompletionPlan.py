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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Synopsis
#
#
# The world is experiencing an increasing incidence of life-changing neuropsychiatric disorder diagnoses but with few pharmaceutical treatment options available, the outlook for sufferers can be daunting. Drug companies and research organisations are investing large sums into therapeutic solutions, but existing methods have to date achieved a very low rate of success.
#
# One recent promising direction for identifying possible therapeutic drugs ar a cheaper and quicker rate is to mine the expansive data available, genetic or otherwise, covering drugs already approved to treat other indications to gauge their suitability for repurposing to treat neuropsychiatric conditions.
#
# To date, several drugs have been successfully repurposed in this way such as thalidomide (repurposed for myeloma and leprosy) and sildenafil (to treat erectile dysfunction). Genetic evidence has been found suggesting that many other opportunities exist to repurpose other drugs in this way.
#
# This study uses gene set enrichment and disease module based network analyses to seek robust genetic evidence which identifies drug repurposing opportunities for several neuropsychiatric diseases - schizophrenia, Alzheimer’s disease and Parkinson’s disease. It analyses better understood medical conditions, hypertension and hypercholesterolaemia to initially validate the data and approaches used.

# # Thesis Outline - thesis 80% complete
# * **Chapter 1: Introduction** - *40% complete*
#     * 1.1 Why is gene set analysis useful?
#     * 1.2 Gene set analysis
#     * 1.3 Drug-target gene interaction data set
#     * 1.4 QTL gene annotation
#     * 1.5 Disease-gene analysis
#     * 1.6 In-silico development environment
#     * Drug repurposing, genetic evidence, analysis techniques

# * **Chapter 2: An analysis pipeline to identify related sets of drug targets associated with disease** - *100% complete*
#
#     * An assessment of the drug target databases that are available (uniqueness and overlap) leading to an objective decision on which should be included in order to capture most information
#
#     * Design and construction of the database of drug targets, linking this to genes, and the analytical pipeline to interrogate GWAS.
#
#     * Creating an analysis pipeline using only the positional approach (as it is conventional) using diseases for which we have a known pathology and known drugs that work (call these Known Diseases)

# * **Chapter 3: Appraisal of the pipeline performance with hypertension and hypercholesterolemia** - *100% complete*
#
#     * If the more comprehensive curation of drug targets present in DUGGIE provides an improvement over other the individual DTI sources available
#     *  the pipeline, chosen annotation scheme and target data set reliably only select drugs that are used to treat the disease

# * **Chapter 4: Using functional QTL gene annotation to better identify genes associated with disease** - *100% complete*
#
#     * Which gene annotation scheme is best (proximity, funtional/QTL or a hybrid)?
#     * Running the analysis pipeline using each of these to analyse which scheme performs better

# * **Chapter 5: Application of the analysis pipeline to neuropsychiatric diseases** - *50% complete*
#
#     * Running the analysis pipeline on several Unknown Diseases
#     * Analyse and interpret the results

# * **Chapter 6: Discussion & conclusions** - *0% complete*
#
#     * Look at creating a scoring system for combining pathway/network and gene set evidence.

# * **Appendix A: Drug-gene interaction datasets** - *100% complete*
#     * A synopsis of the drug-gene interaction datasets chosen, and a brief outline of the data wrangling used.

# * **Appendix B: MAGMA permutation issues** - *100% complete*
#     * An outline of the issues found, and how better parameters were chosen.

# * **Appendix C: Pipeline scripts** - *70% complete*
#     * diagrams and listings of the scripts and jupyter notebooks used for the project.


