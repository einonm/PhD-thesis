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

# # Synopsis
#
#
# The world is experiencing an increasing incidence of life-changing neuropsychiatric disorder diagnoses but with few pharmaceutical treatment options available, the outlook for sufferers can be daunting. Drug companies and research organisations are investing large sums into therapeutic solutions, but existing methods have to date achieved a very low rate of success.
#
# One recent promising direction for identifying possible therapeutic drugs ar a cheaper and quicker rate is to mine the expansive data available, genetic or otherwise, covering drugs already approved to treat other indications to gauge their suitability for repurposing to treat neuropsychiatric conditions.
#
# To date, several drugs have been successfully repurposed in this way such as thalidomide (repurposed for myeloma and leprosy) and sildenafil (to treat erectile dysfunction). Genetic evidence has been found suggesting that many other opportunities exist to repurpose other drugs in this way.
#
# This thesis aims to build upon current gene-set analysis methods, collecting the best available drug interaction data and explore possibilities for utilising new forms of evidence to identify the most promising approved drug repurposing opportunities.
#
# The methods will be utilised to create a set of best-practice analytical pipelines whose output will be validated against known therapeutic drugs for a few well understood diseases - including hypertension, hypterchloresterolaemia and diabetes.
#
# The validated pipelines will be applied to lesser understood neuropsychiatric diseases - Parkinson's, Alzheimer's and Schizophrenia, to identify approved drugs most likely to have a causal therapeutic effect.

# # Thesis Outline - thesis 40% complete
# * **Chapter 1: Introduction** - *30% complete*
#     * 1.1 Gene set analysis
#     * 1.2 Drug-target gene interaction data set
#     * 1.3 QTL gene annotation
#     * 1.4 Disease-gene analysis
#     * 1.5 Network analysis
#     * 1.6 In-silico development environment
#     * Drug repurposing, genetic evidence, analysis techniques

# * **Chapter 2: An analysis pipeline to identify related sets of drug targets associated with disease** - *90% complete*
#
#     * An assessment of the drug target databases that are available (uniqueness and overlap) leading to an objective decision on which should be included in order to capture most information
#
#     * Design and construction of the database of drug targets, linking this to genes, and the analytical pipeline to interrogate GWAS.
#
#     * Creating an analysis pipeline using only the positional approach (as it is conventional) using diseases for which we have a known pathology and known drugs that work (call these Known Diseases)

# * **Chapter 3: Appraisal of the pipeline performance with hypertension and hypercholesterolemia** - *50% complete*
#
#     * If the more comprehensive curation of drug targets present in DUGGIE provides an improvement over other the individual DTI sources available
#     *  the pipeline, chosen annotation scheme and target data set reliably only select drugs that are used to treat the disease

# * **Chapter 4: Using functional QTL gene annotation to better identify genes associated with disease** - *20% complete*
#
#     * Which gene annotation scheme is best (proximity, funtional/QTL or a hybrid)?
#     * Running the analysis pipeline using each of these to analyse which scheme performs better

# * **Chapter 5: Application of the analysis pipeline to neuropsychiatric diseases** - *0% complete*
#
#     * Running the analysis pipeline on several Unknown Diseases
#     * Analyse and interpret the results

# * **Chapter 6: Using biological pathway and network data for drug repurposing** - *0% complete*
#     * Developing the Gene Ontology approach i.e. identifying drugs that interact with a biological Gene Ontology pathway that has been associated with a disease via GWAS (and possibly TWAS)
#     * Assess using Known Diseases and conventional statistical tests
#     * Consider more novel network analyses of the same data and compare to the results obtained using conventional statistical tests
#     * Apply the method to Unknown Diseases

# * **Chapter 7: Discussion & conclusions** - *0% complete*
#
#     * Look at creating a scoring system for combining pathway/network and gene set evidence.

# * **Appendix A: Drug-gene interaction datasets** - *100% complete*
#     * A synopsis of the drug-gene interaction datasets chosen, and a brief outline of the data wrangling used.

# * **Appendix B: MAGMA permutation issues** - *100% complete*
#     * An outline of the issues found, and how better parameters were chosen.

# * **Appendix C: Pipeline scripts** - *0% complete*
#     * diagrams and listings of the scripts and jupyter notebooks used for the project.
