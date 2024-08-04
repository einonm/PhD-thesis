# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
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

import qvalue
from IPython.display import display, HTML, Latex, display_latex

import warnings

warnings.filterwarnings("ignore")

import subprocess

# pip install seaborn
import seaborn as sns
import pylab
import scipy.stats as stats
from scipy.stats import chi2, norm

# %matplotlib inline
import venn

# Also generate tables in latex format, with formatting opts
pd.set_option("display.latex.repr", True)
pd.set_option("display.latex.longtable", True)
pd.set_option("display.latex.multirow", True)
pd.set_option("display.latex.multicolumn", True)
pd.set_option("display.max_colwidth", 180)

# For plots
fig_unit = 10
fig_width = 3
fig_height = 2
# -

display(Latex(r"\doublespacing"))
display(Latex(r"\setlength{\parskip}{5mm}"))
display(Latex(r"\newpage"))

# # Abstract
#
# Drug discovery, the process by which new candidate medicines are found, is an expensive, difficult and lengthy endeavour with a low rate of success. Drug repositioning (or repurposing) involves investigating the effectiveness of drugs to treat diseases other than those for which they are already approved and is seen as a cheaper, quicker way of producing novel therapeutic interventions.
#
# Along with the premise that drugs supported by genetic associations to disease are more likely to successfully advance through clinical trials, the use of modern in-silico screening techniques and rapidly growing databases covering the pharmacological, genomic and transcriptomic domains have uncovered several novel drug repurposing opportunities to date and are considered to be a fertile area for yielding more.
#
# One technique used for repurposing investigations is gene set enrichment analysis, which can identify over-represented sets of genes within a larger pool of genes. This  allows sets of genes which encode the protein targets of a drug to be analysed against a larger group of genes associated with a disease and a significant over-representation of drug target genes indicate that a drug may be a candidate for repurposing. Additional evidence can be uncovered through network-based techniques examining protein-protein and drug-target interaction relationships with the cluster of genes associated with a disease, known as the disease module.
#
# This study will use gene set enrichment and disease module based network analyses to seek robust genetic evidence which identifies drug repurposing opportunities for several neuropsychiatric diseases - schizophrenia, Alzheimer's disease and Parkinson's disease. It will use more well understood diseases, hypertension and hypercholesterolaemia, to initially validate the data and approaches used.
#
# Gene-disease associations will be generated from disease genome-wide association study (GWAS) data, where Single Nucleotide Polymorphisms (SNPs) are mapped to genes using evidence provided by functional quantitative trait loci (QTLs). Drug-gene associations will be collected from existing freely available sources and combined, with the aim to provide increased coverage of both drugs and their gene targets that exceeds any currently available single source.
#
# Managing the complexity and scale of integrating these data sets along with applying a detailed reproducible analysis is understood to be a difficult challenge. Another aim of this study is to describe and share the methods and techniques created to overcome these difficulties.

# ## Gene set enrichment analysis
#
# A gene set enrichment analysis will perform a competitive test for each drug, to ascertain if the set of gene targets for the drug are more strongly associated with the disease gene set than a set of disease genes chosen at random. the MAGMA tool performs this test whilst attempting to correct for known confounding effects.
#
# Two gene set curation methods will be described, one for disease-gene data and another drug-target gene data. For each drug, these will be used as input to the enrichment analysis performed using the MAGMA software tool to eventually produce a list of drugs and their association p-values with the disease. A summary of the process can be seen in figure \ref{fig:summary-analysis}.
#
# <img src="./summary-analysis.png">

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\includegraphics{summary-analysis}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:summary-analysis} - Summary of the main parts of a drug-disease association study. \normalfont{Blue squares represent data sets, green diamonds analysis steps.}}"
    )
)
display(Latex(r"\label{fig:summary-analysis}"))
display(Latex(r"\end{figure}"))

# ## Drug-target gene interaction data set
#
# A drug-target gene interaction database provides the collection of gene sets needed to perform a gene set enrichment analysis against the set of genes associated with a disease. There are many freely available sources of drug-target interactions, all collated for varying reasons with the observation that each may define a drug target and provide evidence for inclusion differently. In order to gain the broadest, most complete set of interactions data from these individual sources will be combined, and in combining them a compatible definition of a drug target should be determined along with consistent schemes for identifying drugs and targets (such as Entrez gene IDs or UniProt protein IDs).
#
# To compile the most comprehensive set of drug gene interactions available, a search of the literature will be conducted for existing published drug-gene interaction databases. A drug target definition will be sought that is both consistent with any definitions used to compile the constituent drug-gene interaction data sets, and with the idea of a drug target that is based on ligand binding supported by experimental evidence, as is aimed for in this study.
#
# Participating drugs will be required to be identified by an Anatomical Therapeutic Chemical (ATC) code, a method of grouping drugs which will provide both a useful ontology whilst allowing only approved drugs to take part in the analysis.
#
# The supporting evidence used to identify each drug-gene interaction will also be recorded to enable to automatic presentation of this evidence if a repurposing opportunity is identified, allowing the drug-gene relationship to be studied and verified further.
#
# The drug-target database thus created will list for each drug its interacting set of genes, and for each drug-gene pair, supporting evidence for the relationship.

# ## Functional gene annotation
#
# Gene annotation in this case is the assignment of Single Nucleotide Polymorphisms (SNPs) present in the genome to genes, based on the SNP's effect or function. A large number of SNPs involved in a Genome Wide Association Study (GWAS) are non-coding but functionally relevant and assignment of these compared to assignment of coding SNPs to genes is less straightforward.
#
# The most accessible strategy for annotation is a positional or proximity based method where a SNP is assigned to the nearest gene on the same chromosome. However, this approach doesn't take into account any further known biological evidence and has been often found to assign SNPs to genes that differ from other sources of evidence such as Quantitative Trait Loci (QTL) expression measurements that can, in many cases, implicate the greatest effect on a gene close but not adjacent to the SNP being considered.
#
# In an attempt to improve the biological basis for SNP annotation, QTL data can be used. A QTL is a region of DNA whose variation correlates with the variation of a quantitative trait such as an mRNA expression level (eQTL) or CpG methylation of a DNA region (mQTL). These DNA loci are identified by a direct association test between SNPs and a quantitative measurement of the phenotype (e.g. methylation or gene expression) from tens or hundreds of individuals, and provide more direct biological evidence than associating SNPs positionally. As these measured traits can vary by tissue type, QTL data sets are often stratified thus and this provides an opportunity to select a subset of the most relevant tissue types if enough aetiology is known of the disease being investigated - reducing the multiple testing required and preserving statistical power. For Parkinson's, brain specific QTLs are an obvious choice and are used here, but more recent research also indicates that the gut and appendix tissues may be of interest to include.
#
# As QTL data may not include results significant enough to identify SNPs for all known genes, a hybrid approach where genes not annotated by any QTL derived SNPs are annotated using positionally derived SNP lists instead - this increases the coverage of available genes, giving more opportunity to find novel associations.
#
# The resultant gene annotation data set can be used as a basis to assign the degree of association with a disease for each gene, given a set of associations of SNPs with the disease such as is produced by a GWAS.

# ## Disease-gene analysis
#
# This step takes a diseases' GWAS summary statistics and the functionally annotated SNP to gene assignments and transforms them into a set of disease to gene association Z-scores suitable for use in a gene set enrichment analysis.
# The MAGMA tool will be used to provide this functionality, which has various statistical options on how to combine SNP based p-values into gene based p-values, as well as estimating correlations that reflect the LD between genes used to compensate for inter-gene dependencies when performing a gene set analysis.
#
#
# ## Network analysis
#
# Gene network graphs containing a known disease module can be used to measure the relative proximity of drug target genes to genes in the module, taking into account known biases such as network incompleteness due to missing interactions and datasets focusing on well studied disease genes. This method will be used to further investigate and refine the evidence for drug repurposing opportunities.

# \newpage
# # Introduction: Drug-target gene database creation
#
# One major component of a drug-disease gene set enrichment analysis is a comprehensive set of drug-target gene interactions (DGIs), as shown in figure \ref{fig:summary-analysis}. DGIs can be found in freely available drug and target databases, collated for a variety of reasons but having drug-gene interactions as one element, with varying types and levels of contributing evidence for the interaction.
#
# This study will seek to gather from as many independent sources a database of drug-gene interactions, along with any supporting evidence provided. It will also attempt to create a drug target definition that is compatible with that applied by contributing data sources and a broad definition chosen for this study of interactions involving drugs that bind to proteins and are backed by experimental evidence.
#
# For a drug-gene interaction to appear in the curated dataset described here, it should be based on experimental evidence (i.e. not predicted) and involve compounds that bind to human proteins. However, not all drug-gene databases define a drug target similarly, if at all. Almost all databases involve a degree of evidence from literature searches, but these vary greatly between simple searches for co-occurring keywords in abstracts, through complex machine learning algorithms (both considered to be predicted evidence for the sake of this study) to curation of the literature evidence by experts, which is considered to be permissible as experimental evidence.
#
# Many also take as a starting point any other similar database that exist at the time, which could lead to circular references and databases based on data sets that are no longer available - which is relevant in the case here where evidence for a proposed interaction is no longer available to validate.
#
# A literature search will be conducted to discover available drug-gene interaction sources. Each viable source will be further researched to ascertain the data availability, target definition, contributing sources (which may be other databases, so as to avoid circular references), and data formats used. From selected sources the drug-gene interactions and supporting evidence will be extracted and merged into one consistent set of interactions suitable for a gene-set analysis, with separate supporting evidence.

# ## Literature search
#
# In order to compile the most comprehensive set of drug gene interactions available, a search of the literature was conducted for existing published drug-gene interaction databases. The search function available within PubMed was used to search the PubMed citation database for items containing the set of terms 'drug', 'gene', 'interaction' and 'database', constraining the results to studies involving the human species and whose publication text is marked as being freely available. The returned search results list consisting of title, data, journal, Digital object identifier (DOI), authors and abstract was downloaded and manually curated to remove publications not actually reporting on database(s) containing drug-gene interactions. Of those remaining publications, any where the advertised source of data was not available at the time of checking, such as a non-functioning website or download link, were also removed.
#
# Any references to other qualifying drug-gene databases present in this filtered list of publications were also considered for inclusion.
#
# This final list of publications was read and references to the source(s) of DGI data (such as PubMed publication text mining or high throughput screening) recorded along with cross-database integrations, where one database includes data extracted from another.
#
# ## Literature search results
#
# The initial PubMed search resulted in a list of 646 publications which after manual curation reduced to 24 papers to study. After removing those which no longer had data available and adding other drug-gene interaction database publications referenced from the set of curated PubMed papers provided a list of 43 publications. These publications described a total of 18 drug-gene interaction databases, listed in table \ref{tab:dgi_dbs}. Note that if a drug-gene database was published or last updated prior to 2010, but contributed data to another database that was published after 2010 then it is included as discarding it may also discard useful meta-data.

#  |DGI DB name | Abbreviation | Reference |
#  | --- | --- | --- |
#  | Therapeutics Target Database | TTD | idd.nus.edu.sg/group/cjttd<citet data-cite="yang_therapeutic_2016"><sup>yang</sup></citet> |
#  | Comparative Toxico-genomic Database | CTD  |  ctdbase.org<citet data-cite="davis_comparative_2015"><sup>davis</sup></citet> |
#  | PubChem | PubChem |  pubchem.ncbi.nlm.nih.gov<citet data-cite="kim_pubchem_2016"><sup>kim</sup></citet> |
#  | DrugBank | DB | www.drugbank.ca<citet data-cite="law_drugbank_2014"><sup>law</sup></citet> |
#  | Search Tool for InteracTions of CHemicals | STITCH | stitch.embl.de<citet data-cite="szklarczyk_stitch_2016"><sup>szklarczyk</sup></citet> |
#  | Target-toxin Database | T3DB  | www.t3db.org<citet data-cite="wishart_t3db:_2015"><sup></sup>wishart</citet> |
#  | SuperTarget | ST  | insilico.charite.de/supertarget<citet data-cite="hecker_supertarget_2012"><sup>hecker</sup></citet> |
#  | PharmacoGenetics Data, PharmGKB | PGKB  | www.pharmgkb.org<citet data-cite="thorn_pharmgkb:_2013"><sup>thorn</sup></citet> |
#  | ChEMBL | ChEMBL |  www.ebi.ac.uk/chembl<citet data-cite="gaulton_chembl:_2012"><sup>gaulton</sup></citet> |
#  | Drug Gene Interaction Database | DGIdb | www.dgidb.org<citet data-cite="wagner_dgidb_2016"><sup>wagner</sup></citet> |
#  | IUPHAR/BPS Guide to PHARMACOLOGY | GtoPdb | www.guidetopharmacology.org<citet data-cite="pawson_iuphar/bps_2014"><sup>pawson</sup></citet> |
#  | International Union of Basic and Clinical Pharmacology db| IUPHARdb | www.iuphar-db.org<citet data-cite="sharman_iuphar-db:_2013"><sup>sharman</sup></citet> |
#  | Drug SIGnatures Database | DSigDB | tanlab.ucdenver.edu/DSigDB<citet data-cite="yoo_dsigdb:_2015"><sup>yoo</sup></citet> |
#  | Binding Database | BindingDB |  www.bindingdb.org<citet data-cite="gilson_bindingdb_2016"><sup>gilson</sup></citet> |
#  | Psychoactive Drug Screening Program | PDSP |  pdspdb.unc.edu/pdspWeb<citet data-cite="roth_multiplicity_2000"><sup>roth</sup></citet> |
#  | Connectivity Map | CMAP | portals.broadinstitute.org/cmap<citet data-cite="lamb_connectivity_2006"><sup>lamb</sup></citet> |
#  | Manually Annotated Targets And Drugs Online Resource | MATADOR  | matador.embl.de<citet data-cite="gunther_supertarget_2008"><sup>gunther</sup></citet> |
#  | SuperLigands | SL | bioinf.charite.de/superligands<citet data-cite="michalsky_superligands_2005"><sup>michalsky</sup></citet> |

display(Latex(r"\caption{Table \ref{tab:dgi_dbs} - Drug-gene interaction databases}"))
display(Latex(r"\label{tab:dgi_dbs}"))
display(Latex(r"\end{longtable}%"))

# The relationships between DGI databases, their primary sources and each other resulting from the literature search is captured in the network diagram in figure \ref{fig:dgidb_network}.
#
# <img src="./dgidb_network.png">

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\includegraphics{dgidb_network}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:dgidb_network} - Network diagram of relationships between the drug-target interaction databases investigated. \normalfont{Pink nodes are primary data sources, green nodes are DTI databases.}}"
    )
)
display(Latex(r"\label{fig:dgidb_network}"))
display(Latex(r"\end{figure}"))

#  From this information, a set of DGI databases were selected to provide coverage of a large part of the graph, taking inclusion of other DGI databases into account. Some data sources were excluded due to issues with extracting or querying bulk data, for example SuperTarget does not have bulk download capability, and the Java application used to run queries was found to be non-functioning. The drug-gene interaction sources chosen were: STITCH, TTD, DrugBank, DGIdb, and DtoGdb. The nature and details of these databases are paraphrased from the literature in the methods and materials section below.

# ## Updates since original search
#
# During the lifetime of this study, many new drug-gene interaction data sources were published. These were considered for inclusion or exclusion from the curated data set. These considerations were:
#
# * Open Targets<citet data-cite="koscielny_open_2017"><sup>koscielny</sup></citet> - excluded as Open Targets links diseases to targets and diseases to drugs but not drugs to targets.
# * Drug Central<citet data-cite="ursu_drugcentral:_2017"><sup>ursu</sup></citet> - Scrapes data from ChEMBL and FDA drug labels, excluded as some study participating data sources already contain data from these primary sources.
# * Drug Repurposing Hub<citet data-cite="corsello_drug_2017"><sup>corsello</sup></citet> - Excluded as main focus is to capture drug chemical properties, but does include drug target data from the literature and other databases which are already include by other participating data sources.
# * KEGG<citet data-cite="kanehisa_kegg:_2017"><sup>kanehisa</sup></citet> -  New DRUG database added, containing targets from Japanese and FDA drug labels. Excluded as these are mostly available from other sources already curated and are unlikely to add new targets.

# \newpage
# # Methods and materials
#
# ## The Anatomical Therapeutic Chemical (ATC) classification system
#
# The ATC classification system<citet data-cite="noauthor_anatomical_nodate"><sup>atc</sup></citet>, along with it's associated Defined Daily Dose (DDD) is designed to exchange and compare data on drug utilisation worldwide, and managed by the World Health Organization (WHO). It catalogs a non-exhaustive set of approved, licensed drugs by a 7 character code which is split into 5 levels:
#
# 1) Anatomical/pharmacological main group
# 2) Pharmacological/therapeutic main group
# 3) Chemical/therapeutic/pharmacological subgroup
# 4) Chemical/therapeutic/pharmacological subgroup
# 5) Chemical substance
#
# The 2nd, 3rd and 4th levels are often used to identify pharmacological subgroups where appropriate, instead of therapeutic or chemical subgroups. For example, the level 5 ATC code (chemical substance) for *metaformin* is **A10BA02**. This substance also belongs to 4 other groups, one at each level, which can also be identified from its level 5 code: **A** is its 1st level, anatomical main group code - *Alimentary tract and metabolism*. **A10** is the 2nd level code, therapeutic subgroup - *Drugs used in diabetes*. **A10B** is the 3rd level, pharmacological subgroup - *Blood glucose lowering drugs, excl. insulin* and finally **A10BA** is the 4th level code, chemical subgroup - *Biguanides*.
#
# In effect, if all level 5 coded substances are taken as a set, this is also the complete set of substances described by all level 4 codes, or level 3 codes etc. but grouped differently. It follows that the set of gene targets of drugs at level 5 is also the same set of genes, grouped differently, at each of the other 4 levels.
#
# Using available drug-target databases it is possible to populate a table of level 5 chemical substances with known protein / gene targets of each drug.
#
# As drugs can have several different therapeutic or pharmacological applications (e.g. ethanol can be used topically or orally), the same substance can appear in different groups at the same level. Also, an ATC code may refer to combinations of two or more drugs - meaning that a direct conversion between ATC codes and chemical identifiers may not always be possible.

# ## Gene identification
#
# The NCBI gene database<citet data-cite="brown_gene:_2015"><sup>brown</sup></citet> contains comprehensive gene-specific records relating to a particular nucleotide reference sequence specific to a species, including pathways, maps, variations, phenotypes and other resources. Each entry in the database is given a unique numeric identifier and as the tool used to search the database is named 'Entrez', these IDs are colloquially referred to as 'Entrez Gene IDs' and are chosen as the default gene identification system for collating drug-gene interaction data sets.
#
# In comparison, the widely used HGNC/HUGO gene names<citet data-cite="gray_review_2016"><sup>gray</sup></citet> are specific to the human species only so were discounted to allow the possibility of analysis of non-human genes. Ensembl gene IDs, whilst comparable and directly translatable to Entrez IDs and having a broader coverage of data categories including transcription variations, were also discounted as a default choice with the aim of minimising the number of lossy conversions between Ensembl and Entrez IDs. However, the final merged drug-gene interaction data sets will be converted to Ensembl gene IDs to use in the gene set enrichment analysis, as the opposing disease associated gene set will be generated from QTL data that is mainly recorded using Ensembl gene IDs.
#
# UniProt identifiers, which often are used to record the targets of drugs also have an official protocol for mapping to Entrez Gene IDs<citet data-cite="noauthor_mapping_nodate"><sup>mapping</sup></citet>. While specific isoforms of a protein, for example, can also be mapped individually to Ensembl IDs from UniProt IDs, their use in conversions between datasets can lead to data being discarded unnecessarily unless the relationship between Ensembl IDs is factored into the process.

# ## PubChem compound IDs
#
# PubChem compound IDs (CIDs) were by far the most prevalent means available to identify drugs in the drug-gene databases used, and hence are taken as the default identification system for drugs for this study. PubChem compound IDs are unique chemical structure identifiers, providing a canonical representation of a compound whilst dealing with ambiguities such as tautomers and salt forms which are otherwise differentiated by PubChem using substance IDs. Unfortunately not all PubChem substance IDs map to a compound ID, one notable example of these are monoclonal antibody drugs which are only identified by substance IDs, their structure being too complex for a compound ID and hence do not appear in this study - removing around 30 drugs from the ATC group L01XC. There is also some indication that PubChem compound IDs have complications with not dealing with all the different structural variations of a drug (approved drugs being especially vulnerable to this) within one ID, but having to use several identifiers<citet data-cite="muresan_mapping_2012"><sup>muresan</sup></citet> - this may translate into some targets unknowingly being dropped from a drug's gene target list.

# ## In-silico development environment
#
# Data manipulation and analysis were carried out in a Linux environment, with big data computation taking place on Linux clustered HPC hardware running the Slurm scheduler and data wrangling activities performed on a personal Linux machine. Orchestration of the computation was achieved using Bash shell scripts with data manipulation and analysis activities coded in the Python language within the Jupyter Notebook environment. Notable Python libraries used for the analysis were Numpy, SciPy, Seaborn and Pandas. All scripts and notebook versions were recorded and stored by the Git version control system, allowing any part of the analysis pipeline to be repeated given identical data inputs.
#
# The source text for generation of this document, along with the in-line python code for creating the tables and graphs from raw results data is also available as a Jupyter notebook.

# ## Drug-gene interaction datasets
#
# The drug-gene interaction data sources described below were obtained from the literature search described previously, either directly as the subject of a publication or from references within a publication. To enable the datasets that are extracted from these sources to be combined, they need to be in a common data format. To this end, each individual extracted dataset will use PubChem Compound IDs to identify drugs with an ATC code if present (at least one data source will need this, or an external PubChem-ATC conversion table will be required), Entrez Gene IDs to identify targets, and an evidence source for each drug-gene interaction - the default for this will be a PubMed identifier (PMID). Two data sets for each source should be produced - the interaction data will list PubChem CIDs (optionally with ATC codes) with each having an Entrez Gene ID list. The evidence table will have one entry per drug-gene interaction, per piece of evidence.

# ### IUPHAR/BPS Guide to Pharmacology Database (GtoPdb)
#
# The IUPHAR/BPS Guide to Pharmacology database (GtoPdb) is an expert-curated database of direct mappings (where the literature data permits) between chemical structures and their primary molecular targets, concentrating on those targets and ligands most relevant to therapeutics and drug discovery. It uses manual expert curation of literature, including posters and patents, as a source of evidence. PubMed identifiers for the evidence and affinity measurements are available directly from the database.
#
# From the GtoPdb website, a plain-text tab-separated file of raw interaction data was downloaded. This listed 19,078 interactions, each having up to 37 attributes such as target IDs (e.g. HUGO gene symbol, UniProt ID and Ensembl ID) along with the interacting species, affinity measurements and ligand IDs (PubChem substance ID, HUGO gene symbol). This list was filtered to only include Human targets and all columns removed except for the target Uniprot ID, ligand PubChem substance ID and evidence PubMed ID, leaving 12,523 entries. PubChem substance IDs (Which describe different forms of a compound) were normalised to PubChem compound IDs using the PubChem ID exchange service<citet data-cite="noauthor_pubchem_nodate"><sup>idx</sup></citet>. With the PubChem substance IDs not converted removed, 10,611 interactions remained. On converting the target UniProt IDs to Entrez gene IDs using the 'idmapping' table available from UniProt<citet data-cite="noauthor_uniprot_nodate"><sup>idmapping</sup></citet>, 10,596 interactions were left to be recorded as DGI evidence. This table was further transformed to be used in further analysis into a table indexed by drug PubChem compound ID of target Entrez IDs for each of the final remaining list of 4,931 drugs (343 with more than 5 targets).
#
# GtoPdb identifies binding ligands and their protein targets, defining a drug target as such.

# ### Therapeutic Target Database (TTD)
#
# The Therapeutic Target Database (TTD) is a dataset providing information about therapeutic targets of approved and under-trial drugs, including other database cross-references, disease mappings, sequencing, structure, and pathway data as well as pharmacodynamic/pharmacokinetic information. Drug and target data are available from the TTD website, or via a set of downloadable tables in plain text, tab-separated files. These files containing target interaction data consist of key-value entries indexed by target ID, with some example key properties being UniProt ID, drug, type of target, name, agonist, binder, Reactome pathway and synonyms.
#
# The raw TTD target data file was downloaded from the TTD website. This was transformed programmatically into a table indexed by the TTD target ID (of which there were 4,456), with each key (numbering 44 in total) present in the database forming a column in the table and their values extracted from the downloaded data appropriately. Targets having a 'Type of target' of 'Successful' were extracted, excluding targets marked as 'Clinical Trial' or 'Research' targets, leaving 629 targets listed. Next all data columns except for 'UniProt ID' and 'Drug(s)' were dropped from the table, leaving just the bare drug-target interaction data. Removing entries with no UniProt IDs or drugs listed gave a table of 464 targets. This table was further transformed from a list of targets against multiple drugs to one listing 14,624 drugs with multiple targets. The UniProt IDs were converted to Entrez gene IDs using the 'idmapping' table available from UniProt<citet data-cite="noauthor_uniprot_nodate"><sup>idmapping</sup></citet>.
#
# The TTD also provides a drug ID cross-matching table, where TTD drug names are listed with cross-matching identifiers from other public databases. The file containing this table, recorded in tab separated plain text key-value pairs indexed by drug ID, was downloaded and transformed into a table indexed by TTD drug ID with one entry per row and 8 columns listing the drug's name, formula and its ID as identified by other databases such as PubChem, ChEBI and the ATC system. The table columns were reduced to remove all except for the PubChem Compound ID (CID) and the ATC code columns. Any entries missing a PubChem CID or at least one ATC code were also removed. This cross-matching drug table was merged with the drug-target interaction data table resulting in a final set of 1,010 ATC-coded approved drugs with their corresponding gene targets recorded as Entrez IDs.
#
# References to the TTD's sources of DGI data are not available in any downloadable tables, but are from the TTD website via each drug's information page, so web scraping tools or a manual check will be utilised to record this information. TTD criteria for a target is described as<citet data-cite="zhu_update_2010"><sup>zhu</sup></citet>: 'The criteria for identifying the primary target of a drug or targets of a multi-target drug is based on the developer or literature reported cell-based or in vivo evidence that links the target to the therapeutic effect of the drug.'

# ### DrugBank
#
# DrugBank contains detailed drug, drug-target, drug action and drug interaction information about FDA-approved drugs as well as experimental drugs going through the FDA approval process. Interaction evidence is provided by expert curation of peer-reviewed papers, with PMIDs available directly from the database.
#
# The full DrugBank database was downloaded from the DrugBank website in XML format. This XML file was parsed to extract a table listing 3,883 drugs by PubMed compound ID against the drug's Anatomical Therapeutic Chemical (ATC) classification codes and protein targets as a list of UniProt IDs. Cleaning drugs with no targets, no ATC codes or no PubChem ID available reduces the number of drugs in the list to 1,481. The target protein UniProt IDs were converted to Entrez gene IDs using the 'idmapping' table available from UniProt<citet data-cite="noauthor_uniprot_nodate"><sup>idmapping</sup></citet> leaving a final list of 1,297 drugs with a PubChem compound ID and at least one ATC code and target.
#
# The initial downloaded XML file was parsed again to extract a table of PubChem compound ID, ATC codes, target and references to evidence for the interaction. This table was saved to a plain-text tab delimited file for later use.
#
# A target, as defined by DrugBank<citet data-cite="noauthor_drugbank_nodate"><sup>db_target</sup></citet> is, to quote: "A protein, macromolecule, small molecule, etc. to which a given drug binds or otherwise interacts with, resulting in an alteration of the normal function of the bound molecule and desirable therapeutic effects or unwanted adverse effects. Some common examples of drug targets include macromolecules like nucleic acids and proteins such as enzymes, ion channels, and receptors. Drugs can be either naturally or synthetically designed to specifically interact with such entities as targets."

# ### DGIdb
#
# The Drug-Gene Interaction database (DGIdb) mines existing DGI resources, providing an interface for searching lists of genes against a compendium of drug-gene interactions and potentially druggable genes, with a leaning towards interactions that are targets of interest to the study of cancers. The DGIdb website<citet data-cite="noauthor_dgidb_nodate"><sup>dgidb</sup></citet> lists it's interation sources as Cancer Commons,  Guide to Pharmacology Interactions, CKB, NCI Cancer Gene Index, Chembl Interactions, My Cancer Genome, TTD, TEND, Tgd Clinical Trial, CIViC, TALC, DoCM, Cancer Genome Interpreter (CGI), Clearity Foundation Biomarkers, Clarity Foundation Clinical Trial, OncoKB, FDA, DrugBank and PharmGKB.
#
# A plain-text tab separated file of 42,727 raw interactions was downloaded from the DGIdb website, which includes interactions from all sources listed above except DrugBank and PharmGKB, which are only available in raw format from their respective websites. The interaction data downloaded lists for each interaction the gene name, Entrez ID, evidence source, PubMed IDs, interaction types, drug name and drug ChEMBL ID. After reducing these columns to gene name, Entrez ID, drug name, drug ChEMBL ID, interaction source and PMIDs then removing entries missing an Entrez ID or ChEMBL ID, 32,480 interactions remained for 6,277 drugs. This table was then transformed into a list indexed by drug ChEMBL ID of multiple drug target Entrez IDs, with 1,216 drugs having 5 or more targets. The set of ChEMBL drug IDs was converted into PubChem compound IDs by submitting the list as a job to the PubChem ID exchange service<citet data-cite="noauthor_pubchem_nodate"><sup>idx</sup></citet>, resulting in a final list of 5,235 drugs listed with their targets (1,090 with 5 or more targets).
#
# The initial raw interaction table was also parsed to extract the interaction source and PubMed ID data, and filtered to contain evidence only for the interactions present in the final drug-target interaction dataset.
#
# Interactions in DGIdb are defined as a relationship between a gene and a drug with an associated interaction type (e.g., inhibitor) from a specified source.

# ### STITCH
#
# The STITCH database records interaction networks of chemicals and proteins and provides a rich interactive network display of the data, which can be downloaded as a set of text based tables and Postgres database SQL dumps. STITCH collates a large number of existing databases and experimental results whilst also augmenting this data with the outputs of an automated text mining pipeline and an interaction prediction algorithm based on chemical similarities. The STITCH database records a confidence score of between 0 and 1 for each drug-gene interaction based on these four categories of supporting evidence - 'experimental', 'database', 'text-mining' and 'prediction'. These scores are combined in a naïve bayseian fashion into a confidence score also within the range 0-1. A higher score means less qualifying interactions, but a lower false-positive rate. A threshold of 0.7 is given as a 'high confidence' limit in the STITCH FAQ<citet data-cite="noauthor_help/stitch:_nodate"><sup>help</sup></citet>, and this value will be used as a threshold in deciding a STITCH drug-gene interaction's inclusion in this study.
#
# Interaction data was downloaded from the STITCH website as a subset of protein-chemical links for the Human species only, in the form of a tab-separated text table identifying chemicals by their PubChem Compound ID and proteins by their Ensembl Protein ID alongside the STITCH assigned scores in the experimental, database, prediction, text-mining categories and the combined score of these four. Only the experimental and database scores were considered, to remove any predicted results. From the initial list of 15,473,939 interactions, removing predicted interactions and rejecting any interaction where the high confidence limit of 0.7 was not reached by either the experimental or database scores resulted in 346,909 qualifying interactions of 132,934 chemical compounds. Next the Ensembl protein IDs of the chemical targets were converted to Entrez gene IDs using the Biomart service and the data was transformed into a final table of 78,478 chemicals listing 241,719 interactions, of which 3,599 chemicals have 5 or more interactions.
#
# Whilst the drug-gene interaction data exists as simple text tables, the interaction evidence can only be obtained individually from the STITCH website or by downloading several Postgres database dump files and importing them into a local Postgres database for bulk querying using Structured Query Language (SQL). Postgres database dump files were downloaded from the STITCH website (following postal acceptance of a free academic license) for the 'items' and 'evidence' database schemas. A local Postgres was created and the dump files imported to it. After studying the tables present in both schemas, an SQL query was created to produce a table of all combinations of PubChem Compound IDs for both drugs and their interacting gene proteins, with each listing a piece of evidence for the interaction. This table was limited to only those drugs having an Anatomical Therapeutic Chemical (ATC) code and excluding non-human species to reduce the size of the resulting evidence to be considered. Although the corresponding STITCH drug-gene interaction data and scores obtained were not filtered for ATC codes (as this data was not readily available), this filtering will occur later on combining drug-gene interaction data from multiple sources.
#
# Evidence for a STITCH interaction is categorised directionally into one of activation, phenotype (phenotypic effects, or: predicted to have the same phenotype), binding, pred_bind (predicted to bind), catalysis, inhibition or reaction<citet data-cite="noauthor_stitch_nodate"><sup>stitch_readme</sup></citet>. Existing databases from which STITCH incorporates evidence are DrugBank, GPCR-ligand database (GILDA), Matador, the Therapeutic Targets Database (TTD) and the Comparative Toxicogenomics Database (CTD). It also uses data from pathway databases such as the Kyoto Encyclopedia of Genes and Genomes (KEGG), NCI/Nature Pathway Interaction Database, Reactome and BioCyc.
#
# As the STITCH database does not provide the definition of a drug target which incorporates targets found in contributing datasets, one will have to be inferred.

# ## Merging drug-gene interaction datasets
#
# The five drug-gene interaction tables obtained from STITCH, TTD, DrugBank, GtoPdb and DGIdb were merged together in turn by joining them based on PubChem CID in an 'outer' join to obtain a table indexed by the union of all PubChem CIDs present in the contributing tables with data columns of Entrez gene IDs, one from each constituent table. Two ATC code columns were also present, one from TTD and another from Databank.
#
# From this table of 85,248 entries, for each drug the five target gene lists were combined into one de-duplicated list and the two ATC code column values were combined together with any entry removed not having an ATC code thus assigned, reducing the table to 1,651 approved drug entries in total. Following the conversion of the Entrez gene IDs to Ensembl gene IDs using the Biomart service, further filtering on the list was performed to remove any drug having less then 5 gene targets. Drugs were also removed where it's gene target list was duplicated by another drug, to cut down on unnecessary multiple testing and resulting in a final yield of 850 approved drug-target gene sets suitable for use in a gene set analysis.

# \newpage
# # Results
#
# #### Defining a drug target
#
# All drug-gene interaction sources studied have differing definitions of a drug target, which are listed below:
#
# * STITCH - Target interactions can only be inferred from the drug-gene interaction sources imported into STITCH. These are DrugBank, GILDA, Matador, TTD and the CTD.
# * DrugBank - defines a target as 'a protein, macromolecule, small molecule, etc. to which a given drug binds or otherwise interacts with, resulting in an alteration of the normal function of the bound molecule and desirable therapeutic effects or unwanted adverse effects. Some common examples of drug targets include macromolecules like nucleic acids and proteins such as enzymes, ion channels, and receptors'.
# * GLIDA - defines a target as a GPCR that binds to a ligand.
# * Matador - Includes 'direct' targets (e.g binding) and 'indirect' (e.g expression changes) targets, but differentiates between the two. Assuming direct targets only are permissible, this is defined concisely by Matador as a 'binding' target.
# * TTD - The criteria for identifying the primary target of a drug or targets of a multi-target drug is based on the developer or literature reported cell-based or in vivo evidence that links the target to the therapeutic effect of the drug.
# * CTD - Has many comprehensive classes of interactions that define a target: abundance, activity, binding, co-treatment, expression, folding, localisation, metabolic processing, mutagenesis, reaction, response to substance, splicing, stability and transport. This data source is incorporated into the STITCH database.
# * DGIdb - a relationship between a gene and a drug with an associated interaction type (e.g., inhibitor) from a specified data source. Also acknowledges that contributing data sources have different definitions of genes and drugs, and provides references to all sources for this reason. As there are 19 known drug-gene interaction data sources used for DGIdb, the definitions of each and compatibility between them is not investigated further here, but left to the investigation of the much smaller set drug-gene interactions appearing in any significant drug-disease association analyses.
# * Guide to Pharmacology - Defines targets as proteins that bind to ligands.
#
# #### Drugs and targets
#
# For each contributing drug-gene interaction source, the number of ATC drugs present, number of ATC drugs with five or more targets and total number of unique targets contributed were calculated and displayed in table \ref{tab:atc_stats}
#
# |DGI source|ATC drugs|ATC drugs >= 5 targets|Total targets|
# | --- | --- | --- | --- |
#  |DrugBank|1245|358|1610|
#  |STITCH  |1239|545|4804|
#  |DGIdb   |1196|509|1622|
#  |GtoPdb  |667|108|452|
#  |TTD     |923|4|279|
#  |Curated|1749|878|5460|

display(
    Latex(
        r"\caption{Table \ref{tab:atc_stats} - ATC drugs and targets contributed by each drug-gene interaction source}"
    )
)
display(Latex(r"\label{tab:atc_stats}"))
display(Latex(r"\end{longtable}%"))

# From a total of 4,823 approved drugs with ATC codes available, 1,749 were annotated with target genes. Following conversion of gene IDs from Entrez to Ensembl, removing drugs with duplicate target lists and filtering of drugs with 5 or more targets listed, 878 approved drugs were eligible to be used in the gene set analysis stage.
#
# Using the set of 1,749 ATC-coded approved drugs, the data can be studied further to ensure it is consistent and as expected. Figure \ref{fig:venn_atc} is a Venn diagram showing the overlap of drugs contributed from each source, for approved drugs only.

# +
DGIs = ["GtoPdb", "TTD", "DrugBank", "DGIdb", "STITCH"]
jointall_df = pd.read_csv("../../data/target-dbs/jointall_targets.tsv", sep="\t")
jointall_df.fillna(value="nan", inplace=True)

# Convert each string of entrez IDs and ATC codes to a set list
for name in ["targets_entrez_" + dgi for dgi in DGIs]:
    jointall_df[name] = [
        jointall_df.at[x, name].split() if jointall_df.at[x, name] != "nan" else []
        for x in jointall_df[name].index
    ]
    jointall_df[name] = [set(x) for x in jointall_df[name]]
    [x.remove("nan") if "nan" in x else x for x in jointall_df[name]]
    jointall_df[name] = [list(x) for x in jointall_df[name]]

jointall_df["all_ATC"] = (
    jointall_df["ATC_TTD"].str.split() + jointall_df["ATC_DrugBank"].str.split()
)
# -

# Concatenate all ATC values per pubChemID - note may have duplicates for different pubchem IDs!
jointall_df["all_ATC"] = [set(x) for x in jointall_df["all_ATC"]]
[x.remove("nan") if "nan" in x else x for x in jointall_df["all_ATC"]]
jointall_df["all_ATC"] = [list(x) for x in jointall_df["all_ATC"]]
jointall_df["ATC"] = [x[0] if x else "nan" for x in jointall_df["all_ATC"]]

# +
jointall_df["all_targets_entrez"] = (
    jointall_df["targets_entrez_TTD"]
    + jointall_df["targets_entrez_DrugBank"]
    + jointall_df["targets_entrez_DGIdb"]
    + jointall_df["targets_entrez_STITCH"]
    + jointall_df["targets_entrez_GtoPdb"]
)

jointall_df["all_targets_entrez"] = [set(x) for x in jointall_df["all_targets_entrez"]]

# +
for index, row in jointall_df.iterrows():
    jointall_df.at[index, "entrez"] = " ".join(
        [str(x) for x in jointall_df.at[index, "all_targets_entrez"]]
    )

# Remove drugs with no targets
jointall_df = jointall_df[jointall_df["entrez"] != ""]

# +
atc_sets = []
set_names = []

for dgi in DGIs:
    new_set = set(
        [
            jointall_df.at[x, "ATC"]
            for x in jointall_df["targets_entrez_" + dgi].index
            if len(jointall_df.at[x, "targets_entrez_" + dgi]) != 0
        ]
    )
    if "nan" in new_set:
        new_set.remove("nan")
    atc_sets.append(new_set)
    set_names.append(dgi + ": " + str(len(new_set)))

# !! Need to sort out duplicates here!
result = jointall_df[jointall_df["ATC"].isin(atc_sets[0])]
dups = jointall_df[jointall_df.duplicated(subset=["ATC"], keep=False)]
dups.sort_values(["ATC"], inplace=True)
# display(dups)

labels = venn.get_labels(atc_sets, fill=["number"])
fig, ax = venn.venn5(labels, names=set_names)
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-15mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:venn_atc} - Venn diagram of approved drug contributions from each drug-gene interaction source. \normalfont{Showing shared and individually contributed drugs by database. Key shows each DGI source with total drugs contributed.}}"
    )
)
display(Latex(r"\label{fig:venn_atc}"))
display(Latex(r"\end{figure}"))

# How much of the final gene set is each DB capturing?
for dgidb_name in DGIs:
    jointall_df[dgidb_name] = [
        int(
            len(jointall_df.at[x, "targets_entrez_" + dgidb_name])
            / len(jointall_df.at[x, "all_targets_entrez"])
            * 100
        )
        for x in jointall_df.index
    ]

atc_only_df = jointall_df[jointall_df["ATC"] != "nan"][DGIs]

# #### Targets
#
# For these approved drugs, what proportion of targets does each drug-target interaction source contribute? Figure \ref{fig:boxplot_atc} is a box plot showing the average and spread of percentage contributions to the targets gathered for each drug, by each source. One data point on this plot corresponds to the percentage of a drug-gene interaction source's contribution to one drug's gene targets. As can be seen STITCH is by far the most prolific contributor, with a sliding scale of average contributions down to the Guide to Pharmacology database with the lowest level of contribution.

graph = sns.boxplot(data=atc_only_df[DGIs]).set(
    xlabel="Contributing database",
    ylabel="Target gene coverage %",
    title="Target gene coverage by contributing database for each ATC drug",
)

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:boxplot_atc} - Boxplot diagram of target gene coverage from each drug-gene interaction source. \normalfont{Each data point is the percentage coverage of all discovered drug targets for a drug by that DGI source.}}"
    )
)
display(Latex(r"\label{fig:boxplot_atc}"))
display(Latex(r"\end{figure}"))

# # Discussion
#
# A drug-gene interaction dataset has been curated with a larger number of ATC-coded drugs than any constituent dataset, an increase of 61% over the largest. The coverage of gene targets is also greater, with 14% more targets involved in interactions for those drugs which are supported by searchable evidence for those interactions.
#
# Constructing a solid set of validated drug-target gene interactions from multiple sources was found to not be a simple task; the approach taken here may not prevent all non-experimentally proven interactions from entering the analysis but does provide for a second quality control step where the gathered supporting evidence of an interaction can be analysed manually for targets involved in a significant drug-disease association of interest.
#
# Attempting to define what a 'target' is that conforms with all the data sources incorporated into the curated drug-target interaction data set is equally difficult - as some other importers of multiple datasets were found to do, one solution is to create a weaker broad definition of a target and defer the scrutiny of the evidence to the end-user - in this case to the drug-gene targets that make up a significant association of a drug with a disease, if any are to be found. This is the approach that will also be taken here, with a drug target defined as a protein that binds to a drug, backed by experimental evidence.

# #### Document provenance

gitspec = subprocess.check_output(
    ["git", "describe", "--all", "--dirty", "--long"]
).strip()
print("This notebook was generated from git revision " + str(gitspec, "utf-8"))
