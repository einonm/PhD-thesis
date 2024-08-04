# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Project background
#
# Studies have demonstrated that for genes which are genetically associated with disease, the proteins for which they encode ate likely to be enriched for existing targets of medications. It follows that these proteins are clear targets for drug development studies and offer opportunities for re-purposing of existing approved drugs. It has recently been demonstrated that drugs supported by genetic associations are more likely to successfully advance through the latter stages of clinical trials.
# Genome Wide Association Studies (GWAS) have now collectively identified common susceptibility variants at a large number of loci that increase risk for neuropsychiatric and neurological disease. Large scale exome and whole genome sequencing studies are also starting to identify low frequency genetic mutations that modify the functional properties of the encoded proteins. To evaluate the clinical impact of these susceptibility alleles, it is now essential to understand their biological significance. Such insight might improve patient diagnosis and prognosis, and with specific relevance to this project, offer the potential to develop novel or improved therapies.
#
# # Research question
#
# This project will develop a novel approach to investigate the intersection between disease risk loci and a curated database of known gene targets of all available therapeutic agents. This gene enrichment analysis will determine pharmaceutical groups of drugs (independently defined by the Anatomical Therapeutic Chemical (ATC) classification system) whose targets are significantly enriched for genetic risk factors to disease (defined by well powered GWAS). This analytical pipeline will be applied to large neuropsychiatric and neurological cohorts that are available at the MRC Centre, for which genomic data is available. It is anticipated that this approach will highlight novel opportunities for therapeutic interventions for neuropsychiatric and neurological diseases: facilitating the re-purposing of available medications that are currently not used to treat the disease or the development of improved drugs.
#
# # hypothesis
#
# Identifying drugs targeting proteins encoded by genes that are associated with increased risk to neuropsychiatric and neurological diseases will reveal novel therapeutic opportunities that can be used to inform and improve the treatment of the disease.
#
# # Aims of research
#
# Following the identification of genetic variants that increase risk to disease, a key goal is to develop therapeutics that modify disease risk in a predictable and beneficial way. Disease associated proteins are clear targets for drug development studies and offer opportunities for drug re-purposing. Moreover, any functional change that can be inferred from the genetic risk variant (e.g. loss of function risk allele) can potentially be exploited as a foundation for drug design. Investigating the intersection between known drug targets and known genetic risk factors can therefore potentially inform future pharmacological interventions.
#
# # Plan of investigation
#
# As is typical with modern software development, an agile, iterative approach will be used to develop the analysis pipeline further, as opposed to a traditional ’waterfall’ approach. It is estimated that each iteration will take 4 months (indicating ~6 more iterations will take place over the life-time of the project). Each iteration will seek to incorporate one or more of the following features; whilst culminating in a working deliverable that has been debugged, optimised and tested to a high degree of confidence:
#
# ### Update the drug targets/ontology database
# To add to the already existing data set obtained from the DrugBank resource, other sources of data and methods of data transformation will be investigated such as utilising the Similarity Ensemble Approach (SEA) and PharmGKB resources, combining them into a single set of drugs and their known or predicted protein targets, organised according to the ATC ontology.
#
# ### Increase the proportion of functionally mapped variants in gene annotations
# Not all gene annotations in the hybrid data set currently used are functionally mapped, so obtaining functionally derived data from other sources would be beneficial to increasing the proportion of functionally derived annotations.
#
# ### Identify further pharmacologically related sets of drug targets that are associated with other diseases
# Now an complete and functioning analysis pipeline is available, this can be quickly run using other annotation sets (such as GO pathway terms) and GWAS results (Alzheimers, Schizophrenia, Huntington's etc...), to both test the pipeline further and also investigate / confirm drug associations for those diseases.
#
# ### Using proximity analysis to identify drugs with therapeutic potential
# Network analysis will be used to further investigate the sets of drugs and drug targets identified by our gene enrichment approach. To achieve this, for each gene present in a gene-set that is associated with disease, we will establish its corresponding protein drug-target within our ontology database and then identify it within the previously developed protein interactome network. We hypothesize that the most effective drugs will target proteins that are, or in the immediate vicinity of, proteins known to be relevant to the pathology of each disease studied. We will apply a network-based proximity analysis to establish the drug-disease proximity within the interactome for each drug-set. This will allow us to infer the relative proximity between each drug linked to the associated gene-set and each disease. Drug-sets whose targets have the shortest path will have the highest proximity score and will be considered as the strongest candidates to have a therapeutic effect.
#
