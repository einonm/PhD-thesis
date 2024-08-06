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
from IPython.display import display, Latex, display_latex

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

import csv

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

# # Curation of a drug-target interaction database

# ## Introduction
#
# This part of the study describes the curation of a drug-target interaction database incorporating the most comprehensive interaction coverage available, used as an input to the gene-set analysis described in chapter \ref{general-introduction}, and explores the improvement and usefulness of the curated database over those already existing. The resulting database has been named DUGGIE, an abbreviation of 'DrUG-Gene IntEractions'. 
#
# As shown in figure \ref{fig:summary-analysis}, along with the disease-gene association Z-scores, the other main input of the drug-disease gene set analysis is a comprehensive set of drug-target gene interactions. Data on such interactions can be found in freely available drug and target databases, collated for a variety of reasons but having drug-gene interactions as one element. These interactions can be listed with varying types and levels of supporting evidence, for example the Toxin and Toxin-Target Database (T3DB)<citet data-cite="wishart_t3db_2015"><sup>wishart</sup></citet> combines target information with toxin data, whilst DrugBank<citet data-cite="law_drugbank_2014"><sup>law</sup></citet> specialises in combining drug data with their chemical properties.
#
# For a drug-gene interaction to appear in DUGGIE, it should be based on expert-curated experimental evidence, not predicted, and involve compounds that have an associated and direct interaction with human proteins. However, not all drug-gene databases define a drug target similarly, if at all. The majority of databases studied involve a degree of evidence from literature searches but these vary greatly between simple searches for co-occurring keywords in abstracts, through complex machine learning algorithms (both considered to be predicted evidence for the sake of this study) to curation of the literature evidence by experts, which is considered to be permissible as experimental evidence. The type of information captured by the included source databases - high throughput screening (HTS), PubMed publications, FDA approvals, pharmaceutical companies and patents are shown in figure \ref{fig:dgidb_network}.
#
# According to the criteria above, the definition of a valid drug target consistent with targets suitable for use in this study can be stated as:
#
#  *A drug target is a protein encoded by a human gene having an experimentally measured, directly associated interaction with the drug.*
#
#
# ## Aims of the chapter
#
# This chapter details activities to gather a database of drug-gene interactions from as many independent, freely-available sources as possible, along with any provided supporting evidence for the interaction. A literature search was undertaken to gather as many potential sources as possible, with a check made for only interactions involving drugs that directly interact with target proteins and are backed by experimental evidence, by recording any drug target definitions applied by contributing data sources and checking them for consistency with the chosen definition of a drug target above. 
#
# Many freely available drug-gene interaction databases were also found to take as a starting point other similar databases available at the time of their curation, which could lead to circular references and databases based on data sets that are no longer obtainable, which is not desirable as all evidence for a proposed interaction should be available to validate. A check of cross-inclusion between the qualifying databases was undertaken to manage these effects. It is also assumed that if a database states or infers a target definition, then all drug-gene interactions included by the database conform to this definition.
#
# After the drug-gene interaction database search and quality checks, a selection of the qualifying databases was made to maximise the number of drug-gene interactions in DUGGIE. Following the non-trivial merging of interaction data from this database selection to form DUGGIE, a brief analysis of the quality and quantity of interactions in DUGGIE over other databases is shown.

# ## Methods and materials

# ### Literature search
#
# A search of the literature was conducted to discover peer-reviewed publications describing currently available drug-gene interaction databases. Each viable source was further researched to ascertain the data availability, target definition, contributing sources (which may be other databases, so as to avoid circular references) and data formats used. From a subset of sources selected to permit the largest coverage, the drug-gene interactions and supporting evidence were extracted and merged into one consistent set of interactions suitable for a gene-set analysis with a separate supporting evidence dataset.
#
# The search function available within the PubMed web site<citet data-cite="noauthor_pubmed_nodate"><sup>pubmed</sup></citet> was used to search the PubMed citation database using the boolean search string '(((drug) AND (gene)) AND (interaction)) AND (database)', constraining the results to studies involving the human species and whose publication text is marked as being freely available. The results were limited to the past 10 years (from March 2019). The returned search results list consisting of title, data, journal, Digital Object Identifier (DOI), authors and abstract was downloaded and manually curated to remove publications not reporting on databases containing drug-gene interactions. Of those remaining publications, any where the advertised source of data was not available at the time of checking, such as a non-functioning website or download link, were also removed. Any references to other qualifying drug-gene databases present in this filtered list of publications were also considered for inclusion.
#
# This final list of publications was read in detail and references to the sources of drug-gene interaction data, such as publication text mining or high throughput screening, recorded along with any cross-database integration - where one database includes data extracted from another.

# ### Checking the source databases match the criteria
#
# As part of the detailed reading of the publications, each source drug-gene database was checked against a number of criteria and any database not meeting them discarded. The criteria are:
#
# * that the database contains substances that are approved drugs
# * the drug-gene interactions are direct and they do not happen through an intermediary, for example by a drug perturbing a gene involved in the same biological pathway as the listed target gene
# * the interactions are not predicted
# * the interaction data is manually curated by experts, not scraped directly by an algorithm into the database
# * the interactions are validated directly by experimental evidence
# * the targets in the drug-gene interactions are human proteins/genes
# * any drug target definition used to qualify inclusion into the database is compatible with the definition used for this study.
#
# To ascertain if each database meets the criteria of including only targets that have a direct interaction with a drug, the relevant database literature and online documentation was searched with any description contributing to a target definition recorded.

# ### The Anatomical Therapeutic Chemical (ATC) classification system
#
# The ATC classification system<citet data-cite="noauthor_anatomical_2017"><sup>atc</sup></citet>, along with its associated Defined Daily Dose (DDD) is designed to exchange and compare data on drug utilisation worldwide, and managed by the World Health Organization (WHO). It catalogues a non-exhaustive set of approved, licensed drugs by a 7 character code which is split into 5 levels:
#
# 1) Anatomical/pharmacological main group
# 2) Pharmacological/therapeutic main group
# 3) Chemical/therapeutic/pharmacological subgroup
# 4) Chemical/therapeutic/pharmacological subgroup
# 5) Chemical substance
#
# The 2nd, 3rd and 4th levels are often used to identify pharmacological subgroups where appropriate, instead of therapeutic or chemical subgroups. For example, the level 5 ATC code (chemical substance) for *metformin* is **A10BA02**. This substance also belongs to 4 other groups, one at each level, which can also be identified from its level 5 code: **A** is its 1st level, anatomical main group code - *Alimentary tract and metabolism*. **A10** is the 2nd level code, therapeutic subgroup - *Drugs used in diabetes*. **A10B** is the 3rd level, pharmacological subgroup - *Blood glucose lowering drugs, excl. insulin* and finally **A10BA** is the 4th level code, chemical subgroup - *Biguanides*.
#
# In effect, if all level 5 coded substances are taken as a set, this is also the complete set of substances described by all level 4 codes, or level 3 codes etc. but grouped differently. It follows that the set of gene targets of drugs at level 5 is also the same set of genes, grouped differently, at each of the other 4 levels. Using available drug-target databases it is possible to populate a table of level 5 chemical substances with known protein / gene targets of each drug.
#
# As drugs can have several different therapeutic or pharmacological applications (e.g. ethanol can be used topically or orally), the same substance can appear in different groups at the same level, and hence there are multiple ATC codes for some substances, in a one drug to many ATC codes relationship. Also, an ATC code may refer to combinations of two or more drugs - there are also relationships of one ATC code to many substances.

# ### Gene identification
#
# The NCBI gene database<citet data-cite="brown_gene_2015"><sup>brown</sup></citet> contains comprehensive gene-specific records relating to a particular nucleotide reference sequence specific to a species, including pathways, maps, variations, phenotypes and other resources. Each entry in the database is given a unique numeric identifier and as the tool used to search the database is named 'Entrez', these IDs are colloquially referred to as 'Entrez Gene IDs' and are chosen as the default gene identification system for collating drug-gene interaction data sets.
#
# In comparison, the widely used HGNC/HUGO gene names<citet data-cite="gray_review_2016"><sup>gray</sup></citet> are specific to the human species only so were not used to allow the possibility of analysis of non-human genes. Ensembl gene IDs, whilst comparable and directly translatable to Entrez IDs and having a broader coverage of data categories including transcription variations, were also discounted as a default choice with the aim of minimising the number of lossy conversions between Ensembl and Entrez IDs. However, the final merged drug-gene interaction data sets were converted to Ensembl gene IDs for use in the gene set enrichment analysis, as the opposing disease associated gene set was generated from QTL data that is mainly recorded using Ensembl gene IDs.
#
# UniProt identifiers, which often are used to record the targets of drugs also have an official protocol for mapping to Entrez Gene IDs<citet data-cite="noauthor_mapping_2019"><sup>mapping</sup></citet>. Whilst multiple isoforms of a protein, for example, can also be mapped individually to the same Ensembl ID from their unique UniProt IDs in a one-to-many relationship, their use in multiple conversions between datasets could lead to data being discarded unnecessarily as this mapping can also involve missing cross-references in one conversion direction only. Hence care was taken to reduce the number of conversions in order to minimise this risk.

# ### PubChem compound IDs
#
# PubChem<citet data-cite="kim_pubchem_2019"><sup>pubchem</sup></citet> compound IDs (CIDs) are by far the most prevalent means available to identify drugs in the drug-gene databases used, hence they are taken as the default identification system for drugs for this study. PubChem compound IDs are unique chemical structure identifiers, providing a canonical representation of a compound whilst dealing with ambiguities such as tautomers and salt forms which are otherwise differentiated by PubChem using substance IDs. Unfortunately not all PubChem substance IDs map to a compound ID, one notable example of these are monoclonal antibody drugs which are only identified by substance IDs, their structure being too complex for a compound ID and hence do not appear in this study - removing around 30 drugs from the ATC group L01XC. There is also some indication that PubChem compound IDs have complications with not dealing with all the different structural variations of a drug (approved drugs being especially vulnerable to this) within one ID, but having to use several identifiers<citet data-cite="muresan_mapping_2012"><sup>muresan</sup></citet> - this may translate into some targets unknowingly being dropped from a drug's gene target list.

# ### Drug-gene interaction datasets
#
# The drug-gene interaction data sources were obtained from the literature search, either directly as the subject of a publication or from references within a publication. To enable the datasets that are extracted from these sources to be combined, a common data format was needed. To this end, each individual extracted dataset uses PubChem Compound IDs to identify drugs and/or an ATC code. Entrez Gene IDs are used to identify targets. Evidence sources for each drug-gene interaction is a PubMed identifier (PMID) by default.
#
# Two data sets for each source were produced - firstly the interaction data lists PubChem CIDs / ATC codes for each drug with an Entrez Gene ID list. Second, the evidence table having one entry per drug-gene interaction, per piece of evidence. A description of each database and an outline of how these tables were extracted can be found in Appendix \ref{appendix-drug-gene-interaction-datasets}.

# ## Results

# ### Literature search results
#
# The initial PubMed search resulted in a list of 646 publications (\hyperref[supp:sf1]{Supplementary File 1}) which after manual curation reduced to 24 papers directly describing drug-gene interaction databases (\hyperref[supp:sf2]{Supplementary File 2}). After removing those which no longer have data available and adding other drug-gene interaction database publications referenced from the set of curated PubMed papers provided a list of 75 publications. These publications describe a total of 19 drug-gene interaction databases, listed in table \ref{tab:dgi_dbs}.
#
# ### Criteria check results
#
# The results of the criteria check are show in table \ref{tab:dgi_dbs}. Several databases were discarded either for having their drug-gene interactions wholly coming from other included databases or a failure to meet all criteria - these are shown in red. Table \ref{tab:targ_defs} shows the results of the target definition search with discarded databases again marked in red.

# \newpage

# +
columns = "|DGI DB Name | Abbr | Approved drugs | Direct interactions | Not predicted | Expert curated | Exp. evidence | Human targets | Reference |"
dbcriteria_df = pd.DataFrame(columns=re.sub(r" ?\| ?", "|", columns).split("|"))
db_criteria_content = [
    '| Binding Database\cite{gilson_bindingdb_2016} | BindingDB | Yes | Yes | Yes | Yes | Yes | Yes | www.bindingdb.org<citet data-cite="gilson_bindingdb_2016"><sup>gilson</sup></citet> |',
    '| ChEMBL\cite{gaulton_chembl_2012} | ChEMBL | Filter | Yes | Yes | Yes | Yes | Yes | www.ebi.ac.uk/chembl<citet data-cite="gaulton_chembl:_2012"><sup>gaulton</sup></citet> |',
    '| Comparative Toxico-genomic Database\cite{davis_comparative_2015} | CTD  | Filter | Yes  | Filter | Yes | Filter | Filter | ctdbase.org<citet data-cite="davis_comparative_2015"><sup>davis</sup></citet> |',
    '| \color{red} Connectivity Map\cite{lamb_connectivity_2006} | CMAP | Yes | \color{red} No | Yes | Yes | N/A | Yes | portals.broadinstitute.org/cmap<citet data-cite="lamb_connectivity_2006"><sup>lamb</sup></citet> |',
    '| DrugBank\cite{law_drugbank_2014} | DB | All | Yes | Yes | Yes | Yes | Yes | www.drugbank.ca<citet data-cite="law_drugbank_2014"><sup>law</sup></citet> |',
    '| Drug Central\cite{ursu_drugcentral_2017} | DC | Filter | Yes | Yes | Yes | Yes | Yes | http://drugcentral.org |',
    '| Drug Gene Interaction Database\cite{wagner_dgidb_2016} | DGIdb | Filter | Yes | Yes | Check contributors | Check contributors | Yes | www.dgidb.org<citet data-cite="wagner_dgidb_2016"><sup>wagner</sup></citet> |',
    '| Drug SIGnatures Database\cite{yoo_dsigdb_2015} | DSigDB | Filter | Yes | Filter | Filter | Yes | Yes | tanlab.ucdenver.edu/DSigDB<citet data-cite="yoo_dsigdb:_2015"><sup>yoo</sup></citet> |',
    '| International Union of Basic and Clinical Pharmacology db\cite{sharman_iuphar-db_2013}| IUPHARdb | Filter | Yes | Yes | Yes | Filter | Filter | www.iuphar-db.org<citet data-cite="sharman_iuphar-db_2013"><sup>sharman</sup></citet> |',
    '| IUPHAR/BPS Guide to PHARMACOLOGY\cite{pawson_iupharbps_2014} | GtoPdb | Filter | Yes | Yes | Yes | Filter | Filter | www.guidetopharmacology.org<citet data-cite="pawson_iupharbps_2014"><sup>pawson</sup></citet> |',
    '| \color{red} Manually Annotated Targets And Drugs Online Resource\cite{gunther_supertarget_2008} (broken link)| MATADOR  | Filter | Filter | Filter | Filter | Filter | Yes | matador.embl.de<citet data-cite="gunther_supertarget_2008"><sup>gunther</sup></citet> |',
    '| \color{red} PharmGKB\cite{thorn_pharmgkb_2013} | PGKB | Filter | From DrugBank/DGIdb | - | - |- | -| www.pharmgkb.org<citet data-cite="thorn_pharmgkb:_2013"><sup>thorn</sup></citet> |',
    '| Psychoactive Drug Screening Program\cite{roth_multiplicity_2000} | PDSP(Ki) | Filter | Yes | Yes | Yes | Yes | Yes | pdsp.unc.edu/pdspweb/<citet data-cite="roth_multiplicity_2000"><sup>roth</sup></citet> |',
    '| \color{red} PubChem\cite{kim_pubchem_2016} | PubChem | Filter|  \color{red} From DGIdb | - | -  | - | - | pubchem.ncbi.nlm.nih.gov<citet data-cite="kim_pubchem_2016"><sup>kim</sup></citet> |',
    '| Search Tool for InteracTions of CHemicals\cite{szklarczyk_stitch_2016} | STITCH | Filter | Filter | Filter | Filter | Filter | Yes | stitch.embl.de<citet data-cite="szklarczyk_stitch_2016"><sup>szklarczyk</sup></citet> |',
    '| \color{red} SuperLigands\cite{michalsky_superligands_2005} | SL | Yes | \color{red} No | No | N/A | No | Yes | bioinf.charite.de/superligands<citet data-cite="michalsky_superligands_2005"><sup>michalsky</sup></citet> |',
    '| \color{red} SuperTarget\cite{hecker_supertarget_2012} | ST  | Filter | Filter | \color{red} No | Filter | Filter | Yes | insilico.charite.de/supertarget<citet data-cite="hecker_supertarget_2012"><sup>hecker</sup></citet> |',
    '| Therapeutics Target Database\cite{yang_therapeutic_2016}  | TTD | Filter | Yes | Yes | Yes  | Yes | Yes | idd.nus.edu.sg/group/cjttd<citet data-cite="yang_therapeutic_2016"><sup>yang</sup></citet> |',
    '| Toxin and Toxin Target Database\cite{wishart_t3db_2015} | T3DB | Filter | Yes | Yes | Yes | Filter | Yes | www.t3db.org<citet data-cite="wishart_t3db:_2015"><sup></sup>wishart</citet> |',
]

for item in db_criteria_content:
    dbcriteria_df.loc[len(dbcriteria_df)] = re.sub(r" ?\| ?", "|", item).split("|")

dbcriteria_df.drop(["Reference"], axis=1, inplace=True)

if formatting == "HTML":
    display(HTML(dbcriteria_df.to_html()))

if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            dbcriteria_df.to_latex(
                column_format="lR{3.9cm}R{1.3cm}p{1.28cm}R{1.5cm}p{1cm}p{1cm}p{1cm}p{1cm}l",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=( r'Drug-gene interaction databases found from a literature search. \normalfont{"Filter" and"Check contributors" implies results '
                        r'need manual filtering before inclusion, otherwise "Yes" or "No". DGI DB Names in red fail to meet all criteria '
                        r'or lack original data.}' ),
                label="tab:dgi_dbs",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
columns = "| DGI DB | Inclusion criteria |"
inclusion_df = pd.DataFrame(columns=re.sub(r" ?\| ?", "|", columns).split("|"))
inclusion_content = [
    '| BindingDB | "..database of measured protein–ligand affinity data" - implied that a target is a protein binding to a drug with affinity measurement evidence\cite{gilson_bindingdb_2016} | ',
    '| ChEMBL | "..targets (the proteins or systems being monitored by an assay)"\cite{gaulton_chembl_2012} | ',
    '| CTD  | "CTD curators developed a hierarchical vocabulary of 50 diverse terms (e.g. binding, phosphorylation, activity) to describe specific molecular interactions between a chemical and gene. A complete list of action codes with their definitions is available via the Help link for interactions on the CTD glossary page or query pages."\cite{davis_comparative_2009}|',
    "| \color{red} CMAP | No direct targets - gene expression signatures only.|",
    "| DB | 'Mechanism of action' listed in database entries - but with many 'unknown' values. \"Drug targets are identified and confirmed using multiple sources (PubMed, TTD, FDA labels, RxList, PharmGKB, textbooks)\"\cite{wishart_drugbank_2006} ; \"565 non-redundant protein/DNA targets being identified for FDA-approved drugs compared to 524 non-redundant targets  dentified in release 1.0. The identification of so many more targets was aided by PolySearch... a text-mining tool developed in our laboratory to facilitate these kinds of searches\"\cite{wishart_drugbank_2008}|",
    '| DC | "Mechanism of action annotations in DrugCentral are aggregated from several expert curated resources: ChEMBL database, Guide to Pharmacology, KEGG Drug, as well as manual curations extracted from drug labels and literature."\cite{ursu_drugcentral_2017}| ',
    '| DGIdb | "Interactions in DGIdb are defined as a relationship between a gene and a drug with an associated interaction type (e.g., inhibitor) from a specified source."\cite{griffith_dgidb_2013}|',
    '| DSigDB | Collects targets in 4 collections, based on different criteria: "D1: approved drugs...we obtained all the approved drugs from US FDA website, and retrieved bioactivity data for these drugs from PubChem and ChEMBL. Genes with ‘active’ bioassay results recorded in these databases were compiled as the drug target genes...D2: large-scale in vitro kinase profiling assays from literature and two databases (Medical Research Council Kinase Inhibitor database and Harvard Medical School Library of Integrated Network-based Cellular Signatures database). We considered the kinase a target of a kinase inhibitor if the IC 50 /K d /K i x 1 mM or the Percent of inhibition over Control x 15\% from the assays."\cite{yoo_dsigdb_2015}|',
    "| IUPHARdb | (An older version of GtoPdb) | ",
    '| GtoPdb | "Our generic use of the term ‘target’ refers to a record in the database that has been resolved to a UniProtKB/Swiss-Prot ID as our primary identifier."\cite{southan_iupharbps_2016} |',
    '| \color{red} MATADOR  | "For a subset of the drug-target relations, namely those where our text-mining approach indicated a wealth of additional information, the type of binding was further analyzed and direct and indirect interactions were manually distinguished. Indirect interactions can, for example, be caused by active metabolites of the drug or by changes in the expression of a protein."\cite{hecker_supertarget_2012} |',
    '| \color{red} PGKB | "Drug information including pharmacological effects, mechanisms of action, and structures was obtained from Drugbank [4] and Pubchem [9]."\cite{thorn_pharmgkb_2013} |',
    '| PDSP(Ki) | "...provides information on the abilities of drugs to interact with an expanding number of molecular targets. The Ki database serves as a data warehouse for published and internally-derived Ki, or affinity, values for a large number of drugs and drug candidates"\cite{noauthor_pdsp_nodate}| ',
    "| \color{red} PubChem | Uses DGIdb as a source of drug-gene interaction data |",
    '| STITCH | "...chemical–protein associations are integrated from pathway and experimental databases, as well as from the literature. Lastly, as many associations as possible are annotated with interaction types", "Experimental evidence of direct chemical–protein binding is derived from the PDSP K i Database (11) and the protein data bank (PDB) (23). Additional interactions between metabolites and proteins are extracted from the manually annotated pathway databases such as KEGG (12), Reactome (14) and the NCI-Nature Pathway Interaction Database...and drug–target relations are imported from DrugBank (8) and MATADOR (10)"\cite{kuhn_stitch_2008}|',
    "| \color{red} SL | No direct targets - predicted by tanimoto similarity only.  |",
    '| \color{red} ST  | "We consider a drug-target relation as a specific interaction of a small chemical compound administered to treat or diagnose a disease and a macromolecule, namely protein, DNA or RNA."\cite{gunther_supertarget_2008} ; "Drug-target relationships described in SuperTarget were obtained in three different ways. Starting with 2400 drugs and their synonyms from the SuperDrug Database (5), the text mining tool EbiMed (6) was used to extract relevant text passages containing potential drug-target relations from about 15 millions public abstracts listed in PubMed. Many thousands of false positive or irrelevant relations were eliminated by manual curation. In parallel, potential drug-target relations were automatically extracted from Medline by searching for synonyms of drugs, proteins and Medical Subject Headings (MeSH terms) describing groups of proteins (7). MeSH terms were used to capture and down-weight interactions that are not explicitly described in the abstracts e.g. for protein families or protein complexes. In the case of families, the specific interacting family member might not be known yet, whereas in the case of complexes, the drug might interact with more than one subunit. Proteins associated to MeSH terms were assigned by a semi-automated procedure relying on mappings provided by MeSH and synonyms of proteins that are aggregated in the STRING resource (8)...The most probable candidates were identified using a benchmarking scheme (8) and manually curated. In a last step, relations from other databases, namely DrugBank (3), KEGG (9), PDB (10), SuperLigands (11) and TTD (4), were checked for drug-target interactions not identified with the preceding steps. If those interactions could be confirmed by literature listed in PubMed, the references were included in SuperTarget otherwise the describing database is referenced."\cite{hecker_supertarget_2012}|',
    '| TTD | Note: Does not link a drug to the target-"The criteria for identifying the primary target of a drug or targets of a multi-target drug is based on the developer or literature reported cell-based or in vivo evidence that links the target to the therapeutic effect of the drug" \cite{zhu_update_2010} |',
    '| T3DB | "Target - Substances to which the toxin chemically binds"\cite{noauthor_t3db_nodate}|',
]

for item in inclusion_content:
    inclusion_df.loc[len(inclusion_df)] = re.sub(r" ?\| ?", "|", item).split("|")

if formatting == "HTML":
    display(HTML(inclusion_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            inclusion_df.to_latex(
                column_format="lp{1.5cm}p{13cm}l",
                multirow=True,
                multicolumn=False,
                index=False,
                caption="Definition of a drug target used for inclusion in each drug-gene interaction database.",
                label="tab:targ_defs",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# ### Cross-inclusion of databases
#
# Further reading was carried out of the collated literature and online content for each database to discover which databases import others in their entirety. The results of this search is listed in table \ref{tab:db_includes}, and displayed in network form along with and primary evidence sources incorporated in each database in figure \ref{fig:dgidb_network}. The last update dates are also shown, giving helpful information for choosing databases to include - if a database includes another which has been updated substantially in the meantime, it may be useful to include both to capture more data. This does not capture corrections to mistakes, however - which may be propagated using this method. It is assumed that this error propagation is of negligible effect.

# +
columns = "| DGI DB | Date of last update | Included Databases |"
inclusions_df = pd.DataFrame(columns=re.sub(r" ?\| ?", "|", columns).split("|"))
inclusions_content = [
    "|  TTD | July 2019 | | ",
    "| CTD | May 2020 | |",
    "| DB | 22-04-2020 | |",
    "| STITCH | November 2015 (v5) | CTD, DB, MATADOR, PDSP, PDB, BinidingDB, PGKB, GLIDA |",
    "| T3DB | Dec 2014 (Estimated) | |",
    "| ChEMBL | March 2020 (release 26) | |",
    "| DGIdb | Jan 2018 | PGKB, DB, TTD, TEND, TALC, MyCancerGenome |",
    "| GtoPdb | Apr 2020 | |",
    "| DSigDB | May 2015 | ChEMBL(D1), PubChem(D1) and all D2 only |",
    "| BindingDB | March 2020 |  |",
    "| DC | August 2018 | ChEMBL, GtoPdb, KEGG Drug  |",
    "| PDSP(Ki) | May 2020 | |",
]

for item in inclusions_content:
    inclusions_df.loc[len(inclusions_df)] = re.sub(r" ?\| ?", "|", item).split("|")

if formatting == "HTML":
    display(HTML(inclusions_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            inclusions_df.to_latex(
                column_format="lp{1.5cm}p{3cm}p{8cm}l",
                multirow=True,
                multicolumn=True,
                index=False,
                caption="Drug-gene interaction database cross-inclusions.",
                label="tab:db_includes",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# The relationships between DGI databases, their primary sources and each other resulting from the literature search is captured in the network diagram in figure \ref{fig:dgidb_network}.

if formatting == "LaTeX":
    display(Latex(r"\newpage"))
    display(Latex(r"\begin{figure}[htpb]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Network diagram of relationships between the drug-target interaction databases "
            r"investigated. \normalfont{Blue nodes are primary data sources, orange nodes are drug-target "
            r"interaction databases, DUGGIE is displayed in red.}}"
        )
    )
    display(Latex(r"\includegraphics{../../img/dgidb_network5}"))
    display(Latex(r"\label{fig:dgidb_network}"))
    display(Latex(r"\end{figure}"))

# <img src="../../img/dgidb_network5.png">

# From this information, a set of drug-gene interaction databases was selected to provide coverage of the largest number of nodes, taking inclusion of other databases into account. The drug-gene interaction sources chosen for further consideration were: STITCH, TTD, DrugBank, DGIdb, DtoGdb, T3DB, DrugCentral and DsigDB. This gives complete coverage of the graph with all database nodes within two edges of the DUGGIE database. Figure \ref{fig:dgidb_network} also shows DUGGIE and this final set of contributing databases with their inclusion relationships and sources.

# ### Merging drug-gene interaction datasets.
#
# The eight drug-gene interaction tables obtained from STITCH, TTD, DrugBank, GtoPdb, DrugCentral, DsigDB, T3DB and DGIdb were merged together in turn by a series of outer table joins based on PubChem CID or ATC code, as detailed in Figure \ref{fig:chapter2-pipeline_2}. From this table of 86,885 entries the eight target gene list columns were combined into one list for each drug and any ATC codes collated, resulting in a two-column table of ATC code lists against target protein lists in a many-to-many relationship.
#
# Separately, the relationships between ATC codes were explored and text processing used to generate three lists covering all assigned ATC codes:
#
#  1. Individual, single drug ATC codes. e.g. prednisolone - A07EA01
#  2. Mixture drug ATC codes, where two or more individual drugs are used together. The ATC codes of the individual drugs used were noted against each mixture ATC code, e.g. prednisolone and promethazine - V03AB05; A07EA01 and D04AA10
#  3. Combination drugs, which are either derivatives of individual drugs or mixture drugs with a broad category and cannot be assigned to individual ATC codes (e.g. psycholeptics). It is assumed this category would either have a target list duplicated with an individual drug or not have an assignable target list unless one is assigned explicitly by a source interaction database. e.g. prednisolone, combinations - A01AC54
#
# This categorisation of the ATC code relationships explains the many-to-many relationship between ATC codes and their corresponding drug targets. Because of the one-to-many directional relationship linking one drug to many ATC codes, there are multiple ATC codes listed for some target gene lists due to combination drug ATC codes. There are also ATC codes where a one-to-many relationship between one mixture ATC code and many drugs exist.
#
# To perform a gene-set analysis against each ATC code, a set of single unique ATC codes is required with each ATC code having a group of target genes. To obtain this, the table was 'exploded' such that each ATC code exists in a separate row. For resulting rows having the same ATC code but differing target lists, as is the case for mixture drugs, these were grouped and their target gene sets combined.
#
# Next the 'exploded' table was split into three according to the individual, mixture and combination ATC lists. The mixture ATC drugs were augmented by targets belonging to their constituent individual ATC drugs. Combination ATC drugs that were identified as combinations of a single ATC drug and have matching target lists were removed as duplicates.
#
# Following the conversion of the Entrez gene IDs to Ensembl gene IDs using the Biomart service, further filtering on the list was performed to remove any drug having less then 5 gene targets to account for any targets not converted. Drugs were also removed where their gene target list was duplicated by another drug, to cut down on unnecessary multiple testing.

# #### Drugs and targets
#
# A drug-target gene interaction database, DUGGIE, was created from 8 source databases that fitted the outlined criteria. The DUGGIE database, as a plain-text tab separated file can be found as \hyperref[supp:sf3]{Supplementary File 3}. The STITCH sub-set of the DUGGIE database, which is the largest contributing database prior to the creation of DUGGIE and used for comparison analyses in subsequent chapters, can be found as \hyperref[supp:sf4]{Supplementary File 4}.
#
# There are a total of 4,970 approved drugs with codes in the current ATC system, of which 4,096 are individual drugs, 202 mixture drugs and 672 combination drugs. From these, 2,977 DUGGIE candidate drugs were annotated with at least one target gene and of these, following the extra augmentation of mixture ATC coded drug targets, 1,984 have 5 or more targets. These ATC coded drugs target 5,600 unique proteins between them, in 64,312 unique interactions.
#
# Along with these figures, table \ref{tab:atc_stats} displays for each contributing drug-gene interaction source the number of candidate drugs and targets for drugs having more than one target, alongside the number of interactions, ATC drugs and targets contributed for drugs having 5 or more targets.

# \newpage

# +
columns = "|DGI source| Drugs >=1 target|targets(drugs >=1 target) | drugs >=5 targets| targets(drugs >=5 targets)|DGIs(drugs >=5 targets)|"
tallies_df = pd.DataFrame(columns=re.sub(r" ?\| ?", "|", columns).split("|"))
tallies_content = [
    "|STITCH      |2111|4799|1128|4784|33673|",
    "|DGIdb       |2080|1622|1053|1592|14990|",
    "|DrugCentral |2262|1296|956|1182|14414|",
    "|DrugBank    |2216|1659|756|1608|10663|",
    "|DsigDB      |1233|1194|857|1179|16019|",
    "|GtoPdb      |1234|464|220|439|2884|",
    "|TTD         |1532|277|5|230|1331|",
    "|T3DB        |814|180|2|174|730|",
    "|DUGGIE    |2977|5637|1984|5600|64312|",
]

for item in tallies_content:
    tallies_df.loc[len(tallies_df)] = re.sub(r" ?\| ?", "|", item).split("|")

if formatting == "HTML":
    display(HTML(tallies_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            tallies_df.to_latex(
                column_format="lp{2cm}p{1.4cm}p{1.8cm}p{1.4cm}p{2cm}p{2cm}l",
                multirow=True,
                multicolumn=False,
                index=False,
                caption=(r"ATC drugs and targets contributed by each source drug-gene interaction database to DUGGIE. "
                         r"\normalfont{For each database, showing the number of drugs and targets of drugs with one or "
                         r"more targets, the number of drugs and targets of drugs with five or more targets, and finally "
                         r"the number of drug-target interactions of drugs with more than 5 targets.}"
                        ),
                label="tab:atc_stats",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# The nature of the contributions to DUGGIE from each constituent database was explored further by creating two correlation triangles - in Figure \ref{fig:drug_jaccard_dbs} showing the Jaccard index, a measure of the set overlap, of the drugs contributed to DUGGIE between each pair of databases and in Figure \ref{fig:target_jaccard_dbs} showing the Jaccard index of the targets contributed to DUGGIE between each pair. Bar graphs of the number of drugs and targets contributed by each database in Figures \ref{fig:drug_db_sizes} and \ref{fig:gene_db_sizes} respectively.
#
# From these, it can be seen that each database has at least some unique contribution to DUGGIE, with some databases  being mainly made up of data extracted from others, such as DrugBank data appearing in DGIdb and STITCH, a relationship corroborated by Figure \ref{fig:dgidb_network}.

if formatting == "LaTeX":
    display(Latex(r"\newpage"))
    display(Latex(r"\begin{figure}[h]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Correlation triangle between database drug contributions using the Jaccard Index as a "
            r"correlation measure. \normalfont{The colour scale represents the Jaccard index value, which also "
            r"appears in the centre of each square.}}"
        )
    )

# +
db_overlap = pd.read_csv(
    "../../../../data/target-dbs/all_dgi_drug_jaccard.csv", index_col=0
)
mask_ut = np.triu(np.ones(db_overlap.shape)).astype(np.bool)
np.fill_diagonal(mask_ut, 0)

plt.figure(figsize=(12, 10))
# sns.set(font_scale=1.5)

sns.heatmap(db_overlap, mask=mask_ut, cmap="viridis", annot=True, square=True)
plt.xlabel("Drug target database")
plt.ylabel("Drug target database")
plt.xticks(rotation=45)
plt.show()
# -

# Note: caption could be here for display under the figure, or in
# previous latex code cell for display above the figure.
if formatting == "LaTeX":
    display(Latex(r"\label{fig:drug_jaccard_dbs}"))
    display(Latex(r"\end{figure}"))

if formatting == "LaTeX":
    display(Latex(r"\newpage"))
    display(Latex(r"\begin{figure}[h]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Bar graph of tallies of drug contribution sizes by database. \normalfont{The same  drug "
            r"can appear in multiple tallies.}}"
        )
    )

# +
with open("../../../../data/target-dbs/all_dgi_drugset_sizes.csv") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        drugset_sizes = row
        drugset_sizes.update((k, int(v)) for k, v in drugset_sizes.items())

plt.figure(figsize=(10, 10))
plt.bar(list(drugset_sizes.keys()), drugset_sizes.values())
plt.xlabel("Drug target database")
plt.ylabel("Number of drugs")
plt.show()
# -

# Note: caption could be here for display under the figure, or in
# previous latex code cell for display above the figure.
if formatting == "LaTeX":
    display(Latex(r"\label{fig:drug_db_sizes}"))
    display(Latex(r"\end{figure}"))

if formatting == "LaTeX":
    display(Latex(r"\newpage"))
    display(Latex(r"\begin{figure}[h]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Correlation triangle between database target contributions using the Jaccard Index as a "
            r"correlation measure. \normalfont{The colour scale represents the Jaccard index value, which also "
            r"appears in the centre of each square.}}"
        )
    )

# +
db_overlap = pd.read_csv(
    "../../../../data/target-dbs/all_dgi_targets_jaccard.csv", index_col=0
)
mask_ut = np.triu(np.ones(db_overlap.shape)).astype(np.bool)
np.fill_diagonal(mask_ut, 0)

plt.figure(figsize=(12, 10))
# sns.set(font_scale=1.5)

sns.heatmap(db_overlap, mask=mask_ut, cmap="viridis", annot=True, square=True)
plt.xlabel("Drug target database")
plt.ylabel("Drug target database")
plt.xticks(rotation=45)
plt.show()
# -

# Note: caption could be here for display under the figure, or in
# previous latex code cell for display above the figure.
if formatting == "LaTeX":
    display(Latex(r"\label{fig:target_jaccard_dbs}"))
    display(Latex(r"\end{figure}"))

if formatting == "LaTeX":
    display(Latex(r"\newpage"))
    display(Latex(r"\begin{figure}[h]"))
    display(Latex(r"\centering"))
    display(
        Latex(
            r"\caption{Bar graph of tallies of gene contribution sizes by database. \normalfont{The same gene "
            r"can appear in multiple tallies.}}"
        )
    )

# +
with open("../../../../data/target-dbs/all_dgi_geneset_sizes.csv") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        geneset_sizes = row
        geneset_sizes.update((k, int(v)) for k, v in geneset_sizes.items())

plt.figure(figsize=(10, 10))
plt.bar(list(geneset_sizes.keys()), geneset_sizes.values())
plt.xlabel("Drug target database")
plt.ylabel("Number of targets")
plt.show()
# -

# Note: caption could be here for display under the figure, or in
# previous latex code cell for display above the figure.
if formatting == "LaTeX":
    display(Latex(r"\label{fig:gene_db_sizes}"))
    display(Latex(r"\end{figure}"))


# create an all_targets_entrez/len_all_entrez columns of de-duped targets for an ATC drug table
def collate_all_targets(in_df):
    df = in_df.copy()
    # Convert each string of entrez IDs and ATC codes to a set list
    for name in ["targets_entrez_" + dgi for dgi in DGIs]:
        df[name] = [x.split() for x in df[name]]
        df[name] = [set(x) for x in df[name]]
        [x.remove("nan") if "nan" in x else x for x in df[name]]
        df[name] = [list(x) for x in df[name]]

    df["all_targets_entrez"] = (
        df["targets_entrez_TTD"]
        + df["targets_entrez_DrugBank"]
        + df["targets_entrez_DGIdb"]
        + df["targets_entrez_STITCH"]
        + df["targets_entrez_GtoPdb"]
        + df["targets_entrez_T3DB"]
        + df["targets_entrez_DrugCentral"]
        + df["targets_entrez_DsigDB"]
    )

    df["all_targets_entrez"] = [set(x) for x in df["all_targets_entrez"]]
    df["len_all_entrez"] = [len(x) for x in df["all_targets_entrez"]]

    return df


DGIs = ["STITCH", "DGIdb", "DrugCentral", "DrugBank", "DsigDB", "GtoPdb", "TTD", "T3DB"]
intermediate_atc_df = pd.read_csv(
    "../../../../data/target-dbs/intermediate_atc_targets.tsv", sep="\t"
)
intermediate_atc_df.fillna(value="nan", inplace=True)
intermediate_atc_df = collate_all_targets(intermediate_atc_df)

# #### Targets
#
# Following the augmentation of the mixture ATC drug targets by individual constituent ATC code target lists, the 1,984 DUGGIE ATC coded drugs were split into tables of 1,502 individual drugs, 309 combination and 173 mixture drugs. 22 other ATC codes were discarded at this point, due to them not being valid codes in the current ATC system and hence having no available description with which to be grouped - these ignored codes could exist because of contributing database errors or historical alterations to the ATC system.
#
# The combination drugs were analysed further to ascertain which have target lists that are duplicates of their corresponding individual ATC coded drugs. 242 of the combination drugs were found to have duplicates, and were discarded leaving a set of 67 combination ATC coded drugs. The three drug groups were then combined back into one list of 1,742 ATC coded drugs with 5 or more targets.
#
# For this list of drugs now in DUGGIE, what proportion of targets does each constituent drug-target interaction source contribute? Figure \ref{fig:boxplot_atc} shows the contribution to identified targets for each drug per database source. One data point on this plot corresponds to the percentage of one drug's gene targets contributed by the drug-gene interaction source. As can be seen STITCH is by far the most prolific contributor, with a sliding scale of average contributions down to the Toxin and Toxin Target Database with the lowest level of contribution.

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[h]"))
    display(Latex(r"\centering"))
    # display(Latex(r'\vspace{-10mm}'))
    display(
        Latex(
            r"\caption{Boxplot diagram of target gene coverage of DUGGIE from each drug-gene interaction "
            r"source. \normalfont{Each data point is the percentage coverage of all targets of a drug "
            r"appearing in DUGGIE by that drug's targets given by the listed DGI source.}}"
        )
    )

# +
# How much of the final gene set is each DB capturing?
for dgi in DGIs:
    intermediate_atc_df[dgi] = [
        int(
            len(intermediate_atc_df.at[x, "targets_entrez_" + dgi])
            / intermediate_atc_df.at[x, "len_all_entrez"]
            * 100
        )
        for x in intermediate_atc_df.index
    ]

fig, ax = plt.subplots(figsize=(10, 6))
graph = sns.boxplot(data=intermediate_atc_df[DGIs], ax=ax).set(
    xlabel="Contributing database", ylabel="Target gene coverage %"
)
# -

if formatting == "LaTeX":
    display(Latex(r"\label{fig:boxplot_atc}"))
    display(Latex(r"\end{figure}"))

# Exploring unique target contributions per contributing database in the same way, figure \ref{fig:boxplot2_atc} displays how successful each drug-gene interaction source is in providing targets unique to that source. Similarly, one data point on this plot corresponds to the percentage of one drug's gene targets contributed uniquely by a source. STITCH has a high contribution of this type, but DGIdb has the largest mean contribution of unique targets. DrugBank's unique contribution is much lower, possibly due to its inclusion in many other databases.

if formatting == "LaTeX":
    display(Latex(r"\begin{figure}[H]"))
    display(Latex(r"\centering"))
    # display(Latex(r'\vspace{-10mm}'))
    display(
        Latex(
            r"\caption{Boxplot diagram of target gene coverage of DUGGIE uniquely contributed by each "
            r"drug-gene interaction source. \normalfont{Each data point is the percentage unique coverage "
            r"of all targets of a drug appearing in DUGGIE by that drug's targets given by in the listed "
            r"DGI source.}}"
        )
    )

# +
unique_df = intermediate_atc_df

unique_df["targets_entrez_no_TTD"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_DrugBank"] = (
    unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_DGIdb"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_STITCH"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_GtoPdb"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_T3DB"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_TTD"]
)

unique_df["targets_entrez_no_DsigDB"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_T3DB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_TTD"]
)

unique_df["targets_entrez_no_DrugCentral"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_T3DB"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_TTD"]
)

for dgi in DGIs:
    unique_df["targets_entrez_no_" + dgi] = [
        set(x) for x in unique_df["targets_entrez_no_" + dgi]
    ]


# (#Targets unique to DB) = #DBtargets - (#DBtargets also in other DBs)
# # %age is ((#Targets unique to DB) / total targets) * 100
for index, row in unique_df.iterrows():
    for DGI in DGIs:
        unique_df.at[index, DGI] = (
            len(
                set(unique_df.at[index, "targets_entrez_" + DGI])
                - set(unique_df.at[index, "targets_entrez_" + DGI]).intersection(
                    set(unique_df.at[index, "targets_entrez_no_" + DGI])
                )
            )
            / unique_df.at[index, "len_all_entrez"]
        ) * 100

fig, ax = plt.subplots(figsize=(10, 6))
graph = sns.boxplot(data=unique_df[DGIs], ax=ax).set(
    xlabel="Contributing database", ylabel="Coverage unique to DB %"
)
# -

if formatting == "LaTeX":
    display(Latex(r"\label{fig:boxplot2_atc}"))
    display(Latex(r"\end{figure}"))

# In the penultimate target gene IDs conversion from Entrez to Ensembl, the vast majority of targets underwent a successful conversion - 24 Entrez IDs failed to have Ensembl ID alternatives, affecting the target lists of 92 drugs but not dropping any below the 5 target threshold for inclusion, so no drugs were discarded and 1,742 drugs remain included in DUGGIE.
#
# Finally, following the removal of drugs which duplicate the target list of other drugs, which would burden the analysis with unnecessary statistical tests, 1,323 approved drugs qualify for inclusion in DUGGIE and hence to be used for this study.

# ### Major Histocompatibility Complex (MHC) considerations
#
# The Major Histocompatibility Complex (MHC) is a complex region of chromosome 6 in humans that contains many immune related genes, exhibiting a high degree of polymorphism and strong linkage disequilibrium (LD). As it can prove problematic in interpreting genomic results in complex human diseases, particularly schizophrenia, it is often excluded from analyses<citet data-cite="pardinas_common_2018"><sup>pardinas</sup></citet>. 
#
# To ascertain the extent to which the MHC region is involved in gene set analyses involving DUGGIE, a list of Ensembl genes overlapping the MHC region on chromosome 6 was compiled (between 26,000,000 and 33,000,000 million base pairs / Mb). Gene lists for each drug present in DUGGIE were searched for the 113 MHC genes, and it was found that only 4 drugs targeted MHC genes - hydrogen peroxide, glycerol, zinc preparations and glutathione, with MHC genes comprising no more than 1% of genes targeted by any of these drugs. This implies that the MHC region is inconsequential in gene set analyses involving DUGGIE and no further consideration is required.

# ## Discussion
#
# This chapter outlined the details of the curation of the most comprehensive drug-gene interaction dataset possible, named DUGGIE, resulting in a larger number of ATC-coded drugs than any of its constituent datasets, with an increase of 76% drugs over the next largest (STITCH). The coverage of gene targets is also greater with 17% more targets involved in interactions for those drugs, and DUGGIE targets are involved in almost 2.2 times more drug-target interactions than STITCH targets whilst still having searchable evidence for these interactions. This provides a more robust basis than any previously created drug-gene interaction database on which to run gene set enrichment analysis pipelines, and will enhance the foundation on which the appraisal of a gene set analysis pipeline with hypertension and hypercholesterolemia will take place, as described in the next chapter, chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}.
