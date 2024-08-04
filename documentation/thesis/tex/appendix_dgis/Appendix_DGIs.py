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
    display(Latex(r"\appendix"))
    display(Latex(r"\doublespacing"))
    display(Latex(r"\setlength{\parskip}{5mm}"))
    display(Latex(r"\newpage"))

# +
#gitspec = subprocess.check_output(
#    ["git", "describe", "--all", "--dirty", "--long"]
#).strip()
#print("This notebook was generated from git revision " + str(gitspec, "utf-8"))
# -


# # Appendix: Drug-gene interaction datasets
#
# The drug-gene interaction data sources were obtained from the literature search, either directly as the subject of a publication or from references within a publication. To enable the datasets that are extracted from these sources to be combined, a common data format was needed. To this end, each individual extracted dataset uses PubChem Compound IDs to identify drugs and/or an ATC code. Entrez Gene IDs are used to identify targets. Evidence sources for each drug-gene interaction is a PubMed identifier (PMID) by default.
#
# Two data sets for each source were produced - firstly the interaction data lists PubChem CIDs / ATC codes for each drug with an Entrez Gene ID list. Second, the evidence table having one entry per drug-gene interaction, per piece of evidence. Below is a description of each database and an outline of how these tables were extracted.

# ## IUPHAR/BPS Guide to Pharmacology Database (GtoPdb)
#
# The IUPHAR/BPS Guide to Pharmacology database (GtoPdb) is an expert-curated database of direct mappings (where the literature data permits) between chemical structures and their primary molecular targets, concentrating on those targets and ligands most relevant to therapeutics and drug discovery. It uses manual expert curation of literature, including posters and patents, as a source of evidence. PubMed identifiers for the evidence and affinity measurements are available directly from the database.
#
# From the GtoPdb website, a plain-text tab-separated file of raw interaction data was downloaded. This listed 19,775 interactions, each having up to 37 attributes such as target IDs (e.g. HUGO gene symbol, UniProt ID and Ensembl ID) along with the interacting species, affinity measurements and ligand IDs (PubChem substance ID, HUGO gene symbol). This list was filtered to only include Human targets and all columns removed except for the target Uniprot ID, ligand PubChem substance ID and evidence PubMed ID, leaving 13,000 entries. PubChem substance IDs (Which describe different forms of a compound) were normalised to PubChem compound IDs using the PubChem ID exchange service<citet data-cite="noauthor_pubchem_2019"><sup>idx</sup></citet>. With the PubChem substance IDs not converted removed, 10,648 interactions remained. On converting the target UniProt IDs to Entrez gene IDs using the 'idmapping' table available from UniProt<citet data-cite="noauthor_uniprot_2019"><sup>idmapping</sup></citet>, 10,601 interactions were left to be recorded as DGI evidence. This table was further transformed to be used in further analysis into a table indexed by drug PubChem compound ID of target Entrez IDs for each of the final remaining list of 4,918 drugs, with 345 of those listed as having more than 5 targets.

# ## Therapeutic Target Database (TTD)
#
# The Therapeutic Target Database (TTD) is a dataset providing information about therapeutic targets of approved and under-trial drugs, including other database cross-references, disease mappings, sequencing, structure, and pathway data as well as pharmacodynamic/pharmacokinetic information. Drug and target data are available from the TTD website, or via a set of downloadable tables in plain text, tab-separated files. These files containing target interaction data consist of key-value entries indexed by target ID, with some example key properties being UniProt ID, drug, type of target, name, agonist, binder, Reactome pathway and synonyms.
#
# The raw TTD target data file was downloaded from the TTD website. This was transformed into a table of 44 features indexed by the TTD target ID, of which there were 4,456. Targets having a 'Type of target' feature with a 'Successful' value were filtered, thus excluding targets marked as 'Clinical Trial' or 'Research'. This left 629 targets. Next all data columns except for 'UniProt ID' and 'Drug(s)' were dropped from the table, leaving just drug-target interaction data. Removing entries with no UniProt ID specified gave a table of 464 targets. This table was further transformed from a list of targets against multiple drugs to one listing 14,624 drugs with multiple targets. The UniProt IDs were converted to Entrez gene IDs using the 'idmapping' table available from UniProt<citet data-cite="noauthor_uniprot_2019"><sup>idmapping</sup></citet>.
#
# The TTD also provides a drug ID cross-matching table, where TTD drug names are listed with cross-matching identifiers from other public databases. The file containing this table, recorded in tab separated plain text key-value pairs indexed by drug ID, was downloaded and transformed into a table indexed by TTD drug ID with one entry per row and columns listing the drug's name, formula and its ID as identified by other databases such as PubChem, ChEBI and the ATC system. The table columns were reduced to remove all except for the PubChem Compound ID (CID) and the ATC code columns. Any entries missing a PubChem CID or at least one ATC code were also removed. This cross-matching drug table was merged with the drug-target interaction data table resulting in a final set of 1,396 ATC-coded approved drugs with their corresponding gene targets recorded as Entrez IDs.

# ## DrugBank
#
# DrugBank contains detailed drug, drug-target, drug action and drug interaction information about FDA-approved drugs as well as experimental drugs going through the FDA approval process. Interaction evidence is provided by expert curation of peer-reviewed papers, with PMIDs available directly from the database.
#
# The full DrugBank database was downloaded from the DrugBank website in XML format. This XML file was parsed to extract a table listing 4,002 drugs by PubMed compound ID against the drug's Anatomical Therapeutic Chemical (ATC) classification codes and protein targets as a list of UniProt IDs. Cleaning drugs with no targets, no ATC codes or no PubChem ID available reduces the number of drugs in the list to 1,496. The target protein UniProt IDs were converted to Entrez gene IDs using the 'idmapping' table available from UniProt<citet data-cite="noauthor_uniprot_2019"><sup>idmapping</sup></citet> leaving a final list of 1,298 drugs with a PubChem compound ID and at least one ATC code and target.
#
# The initial downloaded XML file was parsed again to extract a table of PubChem compound ID, ATC codes, target and references to evidence for the interaction. This table was saved to a plain-text tab delimited file for later use.

# ## Drug-Gene Interaction Database (DGIdb)
#
# The Drug-Gene Interaction database (DGIdb) mines existing DGI resources, providing an interface for searching lists of genes against a compendium of drug-gene interactions and potentially druggable genes, with a leaning towards interactions that are targets of interest to the study of cancers.
#
# A plain-text tab separated file of 42,727 raw interactions was downloaded from the DGIdb website, which includes interactions from all sources listed above except DrugBank and PharmGKB, which are only available in raw format from their respective websites. The interaction data downloaded lists for each interaction the gene name, Entrez ID, evidence source, PubMed IDs, interaction types, drug name and drug ChEMBL ID. After reducing these columns to gene name, Entrez ID, drug name, drug ChEMBL ID, interaction source and PMIDs then removing entries missing an Entrez ID or ChEMBL ID, 32,480 interactions remained for 6,277 drugs. This table was then transformed into a list indexed by drug ChEMBL ID of multiple drug target Entrez IDs, with 1,216 drugs having 5 or more targets. The set of ChEMBL drug IDs was converted into PubChem compound IDs by submitting the list as a job to the PubChem ID exchange service<citet data-cite="noauthor_pubchem_2019"><sup>idx</sup></citet>, resulting in a final list of 5,235 drugs listed with their targets (1,090 with 5 or more targets).
#
# The initial raw interaction table was also parsed to extract the interaction source and PubMed ID data, and filtered to contain evidence only for the interactions present in the final drug-target interaction dataset.

# ## STITCH
#
# The STITCH database records interaction networks of chemicals and proteins and provides a rich interactive network display of the data, which can be downloaded as a set of text based tables and Postgres database SQL dumps. STITCH collates a large number of existing databases and experimental results whilst also augmenting this data with the outputs of an automated text mining pipeline and an interaction prediction algorithm based on chemical similarities. The STITCH database records a confidence score of between 0 and 1 for each drug-gene interaction based on four categories of supporting evidence - 'experimental', 'database', 'text-mining' and 'prediction'. These scores are combined in a naïve bayseian fashion into a confidence score also within the range 0-1. A higher score means less qualifying interactions, but a lower false-positive rate. A threshold of 0.7 is given as a 'high confidence' limit in the STITCH FAQ<citet data-cite="noauthor_helpstitch_2019"><sup>help</sup></citet>, and this value was used as a threshold on the 'experimental' or 'database' evidence categories in deciding a STITCH drug-gene interaction's inclusion in this study.
#
# Interaction data was downloaded from the STITCH website as a subset of protein-chemical links for the Human species only, in the form of a tab-separated text table identifying chemicals by their PubChem Compound ID and proteins by their Ensembl Protein ID alongside the STITCH assigned scores in the experimental, database, prediction, text-mining categories and the combined score of these four. Only the experimental and database scores were considered, to remove any predicted results. From the initial list of 15,473,939 interactions, removing predicted interactions and rejecting any interaction where the high confidence limit of 0.7 was not reached by either the experimental or database scores resulted in 346,909 qualifying interactions of 132,934 chemical compounds. Next the Ensembl protein IDs of the chemical targets were converted to Entrez gene IDs using the Biomart service and the data was transformed into a final table of 78,478 chemicals listing 241,719 interactions, of which 3,599 chemicals have 5 or more interactions.
#
# Whilst the drug-gene interaction data exists as simple text tables, the interaction evidence can only be obtained individually from the STITCH website or by downloading several Postgres database dump files and importing them into a local Postgres database for bulk querying using Structured Query Language (SQL). Postgres database dump files were downloaded from the STITCH website (following postal acceptance of a free academic license) for the 'items' and 'evidence' database schemas. A local Postgres was created and the dump files imported to it. After studying the tables present in both schemas, an SQL query was created to produce a table of all combinations of PubChem Compound IDs for both drugs and their interacting gene proteins, with each listing a piece of evidence for the interaction. This table was limited to only those drugs having an Anatomical Therapeutic Chemical (ATC) code and excluding non-human species to reduce the size of the resulting evidence to be considered. Although the corresponding STITCH drug-gene interaction data and scores obtained were not filtered for ATC codes (as this data was not readily available), this filtering occurs later on combining drug-gene interaction data from multiple sources.

# ## DrugCentral
#
# DrugCentral integrates from a range of online public sources the pharmacologic actions and indications for active pharmaceutical ingredients approved by FDA and other regulatory agencies.
#
# A set of 17,390 drug-target interaction records were downloaded from the DrugCentral website. These records were filtered to extract for each the drug name and interacting HUGO gene ID for homo sapiens only, leaving 12,707 entries. After using a text-matching algorithm to link the drug name to an ATC code, 902 drugs were left. Subsequently converting the HUGO gene IDs to Ensembl IDs using the Biomart service provides a final set of 812 drugs with at least 1 target. Evidence for the interactions is taken from one of ChEMBL, GtoPdb, KEGG or the drug label, and a link to each website is provided via the DrugCentral website.

# ## The Toxin and Toxin-Target Database (T3DB)
#
# The T3DB captures information about clinically relevant toxins. This includes information about drugs and their targets, manually extracted from other databases, government documents, textbooks and journals.
#
# The T3DB toxins database was downloaded from the T3DB website in plain-text XML format. This file was parsed to construct a table of PubChem compound IDs against Uniprot IDs with corresponding pubmed references, containing 3,673 items. The UniProt IDs were converted to Entrez gene IDs using the 'idmapping' table available from UniProt<citet data-cite="noauthor_uniprot_2019"><sup>idmapping</sup></citet>. Any items missing pubmed or Entrez IDs were removed, leaving a total of 3,070 drugs with their interacting genes.

# ## Drug Signatures Database (DsigDB)
#
# DsigDB is a collection of drug and small molecule related gene sets based on quantitative inhibition and/or gene expression changes data. It is split into four distinct parts, labelled D1 to D4. D1 lists approved drugs taken from FDA data. D2 given kinase inhibitor data from pharmaceutical databases. D3 contains perturbagen signatures whilst D4 is a collection of drug signatures from other databases (TTD and CTD), along with additional data from automated text mining. Because of the conflict with the target definition employed for this study, only data from the D1 and D2 collections are valid for use.
#
# The entire detailed database was downloaded in plain text format from the DsigDB website and filtered to obtain the D1 and D2 collections, consisting of 28,771 rows of drug name and single HUGO ID interacting genes. The drug names were merged with any corresponding ATC codes using a text-matching algorithm, resulting in 417 unique ATC-coded drugs with at least one target. The HUGO gene IDs were converted to Ensembl IDs using the Biomart service.
