# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + [markdown] slideshow={"slide_type": "slide"}
# # A genomic medicine approach to identifying novel drug targets
#
# ## Mark Einon - PhD project
#
# ### Overview
#
# * A gene set analysis method is used to calculate the strength of association between a disease and a drug, for many approved drugs
# * Hence this approach can be used to rank the potential of the approved drugs for repurposing to treat a disease
# * A drug's association with a disease is governed by the set of genes that are experimentally known targets for the drug
# * The genetic mechanisms behind the significant associations can be studied further by looking at the target genes involved
# * Further mechanisms involving multiple intermediate gene interactions can be explored using networks of gene interactions.
#
# <small> A genomic medicine approach to identifying novel drug targets - Mark Einon &lt;einonm@cardiff.ac.uk&gt;</small>

# + [markdown] slideshow={"slide_type": "slide"}
# * The gene set analysis method can be broadly described as:
#     * First, functional SNP-gene annotation data is used to associate genes with disease associated SNPs from GWAS results in a 'gene analysis' step
#     * Then open drug target data is used to create a set of genes targeted by each drug. These gene sets for all approved drugs are compiled into a drug-target interaction (DTI) database
#     * Finally the DTI database and the disease associated gene set are used as input to a gene set enrichment analysis, which gives a p-value association of each drug with the disease
#
# <small> A genomic medicine approach to identifying novel drug targets - Mark Einon &lt;einonm@cardiff.ac.uk&gt;</small>

# + [markdown] slideshow={"slide_type": "slide"}
# # Method pipeline overview
# <img src="../thesis/img/summary-analysis-ud.png">
#
# <small> A genomic medicine approach to identifying novel drug targets - Mark Einon &lt;einonm@cardiff.ac.uk&gt;</small>

# + [markdown] slideshow={"slide_type": "slide"}
# # Strategy and aims
#
# * The aim is to apply this pipeline to neuropsychiatric diseases; schizophrenia, Alzheimer's disease and Parkinson's disease, providing evidence for repurposing opportunities
# * To validate the pipeline, the better-understood diseases hypercholesterolaemia and hypertension are analysed and the set of significantly associated drugs, if any, compared with those prescribed to treat the condition.
#
# <small> A genomic medicine approach to identifying novel drug targets - Mark Einon &lt;einonm@cardiff.ac.uk&gt;</small>

# + [markdown] slideshow={"slide_type": "slide"}
# # How to define a drug target?
#
# * The drug-target gene interaction (DTI) database is compiled from a number of different sources, with each source having it's own unique set of supporting evidence for interactions and subsequent explicit or implied definition of a drug target
# * These DTI sources may include snapshots of other sources, which again may have differing target definitions...
# * Based on project requirements and a broad view of available DTI sources, several criteria were chosen to help define a drug target, resulting in a definition for the project of:
#
#     * *A drug target is a protein encoded by a human gene having an experimentally measured, directly associated interaction with the drug*
#
#
# * **Is this a valid definition and useful in further drug repurposing work?**
#
# <small> A genomic medicine approach to identifying novel drug targets - Mark Einon &lt;einonm@cardiff.ac.uk&gt;</small>

# + [markdown] slideshow={"slide_type": "slide"}
# # Criteria used to qualify a DTI source for inclusion
#
# * That any drug target definition defined and used to qualify inclusion into the database is compatible with the definition used for this study.
# * that the database contains substances which are approved drugs
# * the drug-gene interactions are direct and they do not happen through an intermediary, for example by a drug perturbing a gene involved in the same biological pathway as the listed target gene
# * the interactions are not predicted
# * the interaction data is manually curated by experts, not scraped directly by an algorithm into the database
# * the interactions are validated directly by experimental evidence
# * the targets in the drug-gene interactions are human proteins/genes
#
# <small> A genomic medicine approach to identifying novel drug targets - Mark Einon &lt;einonm@cardiff.ac.uk&gt;</small>
