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
from IPython.display import display, Latex, display_latex

display(Latex(r"\doublespacing"))
display(Latex(r"\setlength{\parskip}{5mm}"))

display(Latex(r"\Large \vspace{5mm}"))
# -

# **Summary**

display(Latex(r"\normalsize \vspace{2mm}"))

# Gene set enrichment analysis is an established method of identifying sub-sets of related genes that are significantly associated with a disease whose association with human genetic variation has been explored in a genome-wide association study (GWAS). This method, alongside mining the ever increasing wealth of experimentally obtained drug-target gene sets and functional genomic datasets available, gives rise to the possibility of identifying new repurposing opportunities for approved drugs. Such repurposing opportunities are particularly sought after as therapeutic solutions for debilitating neuropsychiatric and neurodegenerative diseases, for which available treatments are scarce. 
#
# This thesis describes the creation, evaluation and use of an in-silico gene set enrichment pipeline to discover novel evidence for repurposing already approved drugs to treat six neuropsychiatric and neurodegenerative diseases - schizophrenia, major depressive disorder, bipolar disorder, Alzheimer’s disease, Parkinson’s disease and Huntington’s disease. Several strategies were investigated with the aim of improving the yield and accuracy of pipeline results - curating the most comprehensive dataset of drug-gene interactions possible and exploring the use of more relevant functional QTL data with which to annotate genes using available disease GWAS results. The application of these strategies were evaluated on better understood diseases, hypercholesterolemia and hypertension, for which drug treatments are more plentiful.
#
# The strategies investigated were found to be beneficial in identifying a greater number of significantly disease-associated drugs, with nuanced results leading to the adoption of a further strategy of running a battery of analyses and ranking identified candidate repurposing drugs by the number of analyses in which they appear. Executing this battery of analyses gave extensive results, including a number of novel repurposing opportunities, for all six neuropsychiatric and neurodegenerative diseases studied.

display(Latex(r"\singlespacing"))
