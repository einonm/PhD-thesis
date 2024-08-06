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
from IPython.display import display, HTML, Latex

import warnings

warnings.filterwarnings("ignore")

import subprocess

# #!pip install seaborn
import seaborn as sns
import pylab
import scipy.stats as stats
from scipy.stats import chi2, norm

import statsmodels.api as sm

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score

# %matplotlib inline
# #!pip install venn
import venn

sys.path.append("../../../../lib/")

import results_display as rdisp

import re

# #!pip install upsetplot
import upsetplot as usp

# Also generate tables in latex format, with formatting opts
pd.set_option("display.max_colwidth", 2000)
pd.set_option('display.precision', 2)

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


# # Using functional QTL gene annotation to better identify genes associated with disease

# ## Introduction
#
# Proximity annotation is a popular and effective method of assigning a SNP to a gene in order to subsequently transform SNP-disease associations generated from GWAS studies into gene-disease associations. Proximity annotation appoints a SNP to its nearest gene(s), when a SNP's closest gene is defined as either the gene where the SNP is located between the gene start and stop codons or where the SNP occurs within a pragmatic set distance either side of the gene. This is a simple but crude approach to take, based on the assumption that a SNP near a gene will most likely affect the function of that gene - however, this is not always the case, with a previous estimate that at least 10% of SNPs found to influence the expression level of a gene are not located within 15 kilobases (kb) of the gene<citet data-cite="pickrell_understanding_2010"><sup>pickrell</sup></citet>. Functional QTL annotation is an attempt to add greater biological relevance to the assignment of SNPs to genes, given that a large proportion of risk loci are located in non-coding regions of the genome<citet data-cite="hindorff_potential_2009"><sup>hindorff</sup></citet>, and that trait-assoicated SNPs are more likely to be eQTLs<citet data-cite="nicolae_trait-associated_2010,hormozdiari_leveraging_2018"><sup>nicolae,hormozdiari</sup></citet>.
#
# In the current study expression and methylation trait evidence were used to functionally link SNPs to the gene most associated with the trait. However, this scheme has a few disadvantages. The first is that not all SNPs found to have an association with the disease have a functional link to a gene captured by known QTL annotation data, meaning that some genes will have fewer SNPs functionally annotated than proximally annotated and hence a proportion of SNP-disease associations will be discarded. The second is that as QTL trait data from more than one tissue type is used, it is more likely that a SNP may be functionally assigned to multiple genes, adding to statistical noise.
#
# To combat the first disadvantage, a hybrid approach was taken where QTL trait data was combined with proximity annotations in two ways; once by creating a 'missing gene' hybrid annotation by only including proximity annotations for genes which had no resulting functional annotation, and secondly by a 'combined' hybrid annotation method where all SNPs from the proximity annotation were taken and added to the functional annotations for every gene. Both methods maximise the number of genes able to take part in a competitive gene set analysis, with each annotation method producing a different SNP list for each gene. The missing-gene hybrid approach prioritises the functional evidence over the proximal evidence, whilst the combined hybrid approach prioritises the larger proximity data set but is even more likely to assign a SNP to multiple genes.
#
# Any method assigning p-value associations to genes using SNP p-values obtained from GWAS summary results must also deal with the effects of linkage disequilibrium (LD). LD arises from the nature of meiosis, a type of cell division that produces gametes, occurring when sections of DNA containing many SNPs are routinely combined as one entity into each gamete, the distribution of the sections is not random and hence the SNPs within a block are not independently inherited. Because of the random cutting into sections, over many meioses these independent blocks are divided up into smaller and smaller regions, but never disappear. LD means that alleles at a causal locus may be correlated with those at nearby marker loci. This may result in several non-causal SNPs showing significant association with the trait, hence adding noise to the analysis. The true causal SNPs may not even appear in the SNPs that are studied due to incomplete SNP coverage, especially if the causal variants are rare. The MAGMA gene analysis attempts to account for LD by using a reference set of genotype data matching the ancestry and SNP positions used for the GWAS summary results providing the SNP p-values.

# ## Aims of the chapter
#
# This chapter describes and discusses gene set enrichment analyses that were performed with the drug-gene interaction databases derived in chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}, namely DUGGIE and STITCH, and utilising a range of methods for annotating GWAS derived SNPs to genes. A comparison of the results of these analyses was also carried out indicating the most effective annotation methods for identifying drugs whose targets are enriched for GWAS signal, and thus potential candidates for repurposing. The most successful annotation methods will subsequently be used to identify potential repurposable drugs for neuropsychiatric diseases, described in the following results chapter.
#
# In addition to providing further evidence to test the hypotheses that treatment drugs are enriched in gene set analysis pipeline results and that the DUGGIE drug-gene interaction database can identify more significant drugs and treatment drugs than previously collated databases, this chapter tests and discusses the following hypotheses:
#
# * Functionally annotating genes increases the number of significant drugs and treatment drugs in analysis pipeline results
# * Augmenting functional annotations with positional data increases the number of significant drugs and treatment drugs in analysis pipeline results
# * Using QTL data obtained from larger numbers of tissue samples increases the number of significant drugs and treatment drugs in analysis pipeline results
# * Using the specific brain tissue mQTL data alongside eQTL data increases the number of significant drugs and treatment drugs in analysis pipeline results
#
# Several schemes were employed to annotate SNPs to genes to enable the assignment of gene-disease association p-values by the pipeline and the subsequent results compared. The annotation methods explored were proximity, functional (QTL based), and two hybrid annotations of the proximity and functional methods. For each disease considered, each of these methods (excluding proximity) was applied to functional eQTL data derived from blood, brain tissue, brain tissue including frontal cortex mQTL data, all available tissues and all available tissues including mQTL frontal cortex data. This gave a collection of sixteen pipeline execution results differing only in the annotation method and data used to translate disease GWAS derived SNP association p-values into gene p-values. These schemes are shown in in table \ref{tab:db_annotations}.
#
# All sixteen analyses are executed using the DUGGIE drug-gene interaction dataset. To provide additional data points and further evidence exploring the differences of DUGGIE and STITCH, these sixteen annotations were also run with STITCH as the input drug-gene interaction dataset, resulting in a multi-dimensional matrix of 32 results to analyse.
#
# Hypertension and hypercholesterolemia summary GWAS results were again used as baseline diseases for the comparison due to their complex genetic aetiology and plethora of available treatment drugs. A description of these diseases can be found in the introduction to chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}, with 124 treatment drugs and 1194 non-treatment drugs for hypertension involved in the analysis pipeline and 40 treatment drugs and 1278 non-treatment drugs for hypercholesterolemia. The same metrics of interest as used in chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia} were used, with a key metric used of the count of significant drugs supported by an analysis of drugs, broken into treatment and non-treatment sets for each disease. The analysis once again used confusion matrices to measure the capability of each annotation method to identify treatment drugs and reject non-treatment drugs by calculating the Mann-Whitney test statistic, area under the curve of the receiver-operator characteristic, sensitivity, specificity and false discovery rates. The nature of the overlap of significant identified drugs between results sets was also examined.
#
# Alongside blood tissue and all tissue data groups, brain tissue was also analysed as it will be of importance for investigating brain diseases in further chapters, despite it having little known specific relevance to hypertension and hypercholesterolemia. It also provides an intermediate number of functionally annotated genes and tissue sample sizes between that of blood tissue, which has a lower number of annotated genes but higher sample numbers and all tissues, which has the highest number of annotated genes but with eQTLs generated from smaller sample sizes. This will allow the effects of varying numbers of annotated genes and sample sizes on enrichments to be explored.

# ## Methods and Materials
#
# ### Annotation methods
#
# Genome annotation is the process of assigning a property to a genetic location. In this chapter, SNP locations are assigned to genes where the SNP variation is found to have an association with the function or expression of the protein encoded by the gene.
#
# #### Proximity annotation method
#
# This method simply assigns a SNP to a gene based on their relative location, if the SNP is within the gene or a set genomic distance before the transcription start site or after the end of the final exon, taken as within 10 kilobases in this study.
#
# #### Functional QTL annotation method
#
# Functional QTL annotation links a region of DNA to a quantitative trait. Mostly for this study, the region of DNA is a SNP and the trait of interest is gene expression - SNPs are assigned to one or more genes based on association strength between measured gene expression and SNP variation. In addition to the eQTL datasets from various tissues, an mQTL dataset is also used, where SNPs are associated with genes using a methylation quantitative trait. Here SNP variation is firstly associated with cytosine/guanine nucleotide pairs, known as CpG sites (CpGs), then CpG sites that are known to cause silencing or activation of particular genes allow the CpGs-SNP associations to be translated into SNP-gene assignments.
#
# ####  Hybrid annotation methods
#
# Unfortunately, as functional QTL data is based on experimental evidence with incomplete coverage of all genes, the number of functionally annotated genes available is far less than that for proximity annotation. This alone reduces the power of a gene analysis result and makes comparing functional and proximal annotations difficult. To mitigate this effect and allow such comparisons, two hybrid annotation methods were employed, combining functional and proximal annotation data.
#
# The first is a 'missing gene' hybrid approach, where proximity annotation data is only used for genes where there is insufficient functional annotation data. So genes have either SNPs assigned according to functional evidence, or in the absence of sufficient functional evidence, SNPs assigned according to their proximity to the gene.
#
# The second method used is a 'combination' hybrid approach. For genes with only proximity data or only functional data, that data is used to annotate the gene. Where a gene has both proximity and functional SNP annotation data, these two sets of SNPs are combined by collating all SNPs in one list and removing any duplicates.

# ### Datasets
#
# #### Source QTL datasets
#
# ##### Genotype-Tissue Expression (GTEx) eQTLs
#
# Tissue specific expression QTL data from GTEx analysis version 8, obtained from the GTEx portal <citet data-cite="noauthor_gtex_nodate"><sup>gtex</sup></citet> with dbGap accession number phs000424.v8.p2 were used as a source of functional SNP-gene associations. Data from 48 non-diseased tissue types were used. The GTEx data was generated<citet data-cite="noauthor_gtex_nodate-1"><sup>gtex_method</sup></citet> using genotype data based on whole genome sequencing of 838 donor samples that also have RNA-seq data available. With a cis-eQTL defined as an eQTL within 1Mb of a gene transcription start site (TSS), cis-eQTL mapping was performed by the FastQTL package using adaptive permutations and a false discovery rate threshold of 0.05 or less. 
#
# ##### Methylation QTL
#
# Genome-wide methylation profiles from the substantia nigra and the frontal cortex of 134 individuals with Parkinson's disease were obtained<citet data-cite="kia_identification_2021"><sup>pd_mqtl</sup></citet>. cis-methylation QTLs, defined as correlations between a SNP and DNA methylation levels of a CpG site within a 500kb window of the SNP. Five covariates were included and the strongest SNP-CpG pairs were retained at a false discovery threshold of 0.05 or less. The CpG sites were then mapped to genes if they were within 10kb of a gene TSS or end base position. 37,460 CpG sites were eligible for inclusion in this way. Only the frontal cortex data was used in this study, as substantia nigra tissue is likely to be modified from Parkinson's disease.
#
# ##### PsychENCODE QTL
#
# The PsychENCODE project has identified over 2.5 million cis-eQTLs from 1387 samples of brain tissue from the prefrontal cortex, cerebellar cortex and thalamus cortex, gathered from related projects such as ENCODE, the CommonMind Consortium, GTEx and Roadmap.<citet data-cite="noauthor_comprehensive_nodate"><sup>psychencode</sup></citet>. The eQTLs were generated using the standard GTEx pipeline<citet data-cite="noauthor_gtex_nodate-1"><sup>gtex_method</sup></citet>, and the QTLtools package for eQTL identification, with a multiple testing correction made by limiting FDR values to less than 0.05.
#
# ##### eQTLGen
#
# The eQTLGen project has discovered cis-eQTLs from a large-scale meta analysis on 31,684 blood samples from 37 eQTGen Consortium cohorts, identifying cis-eQTLs for 16,987 genes. QC is performed using the published eQTLMappingPipeline<citet data-cite="vosa_large-scale_2021"><sup>eqtlgen</sup></citet> method.
#
# #### Drug-gene interaction datasets DUGGIE and STITCH
#
# The curated drug-target gene interaction database described in chapter \ref{curation-of-a-drug-target-interaction-database}, DUGGIE, is used as the principle source of drug-gene interactions for the gene set analyses in this chapter. As a source of more evidence to support the conclusions reached in chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia} and as it is easily supported by the analysis pipeline, all executions using DUGGIE were duplicated exactly, but with the STITCH drug-gene interaction database in place of DUGGIE.
#
# #### Proximity mapping data
#
# A proximity mapping reference file provided with MAGMA and obtained from the MAGMA website<citet data-cite="noauthor_magma_nodate-1"><sup>magmaweb</sup></citet> was used, created from the 1000 Genomes project European sub-population data and used as one of the pair of input gene sets in the gene-set analyses.
#

# #### Disease treatment drug sets
#
# The DrugBank website<citet data-cite="law_drugbank_2014"><sup>db</sup></citet> was used to create a list of drugs used to treat a disease by using the website's search mechanism to search for drugs marked as treatment for an indication of 'High blood pressure (hypertension)' in the case of hypertension and 'Heterozygous Familial Hypercholesterolaemia' in the case of hypercholesterolemia. Drugs in these treatment sets were removed where evidence was found indicating that they are not currently used in practice or withdrawn.

# #### Annotation data sets
#
# Several annotation data sets were compiled to explore the effect of using functionally annotated GWAS results by comparing against a standard positional annotation scheme. Source QTL datasets were collected together in three distinct groups differentiated by tissue type - brain, blood and all tissues, with the brain and all tissues groups also being split in two, one where the methylation QTL data was added, and one without - with each of these groups further split according to the 'missing gene' and 'combination' hybrid methods utilising positional annotation data, plus with a no combination set, resulted in fifteen different analyses considered in total against a positionally annotated analysis. The all tissue annotation does not include QTLs sourced from brain tissue GTEx data sets as these are already incorporated into PsychENCODE. Table \ref{tab:db_annotations} details all sixteen datasets, along with the simple statistics for each of total number of genes, the mean number of SNPs per gene and standard deviation of SNPs per gene. 

# \newpage

# +
columns = "| Annotation | QTL sources | Hybrid | Total genes | Mean SNPs/gene | Std SNPs/gene |"
annotations_df = pd.DataFrame(columns=re.sub(r" ?\| ?", "|", columns).split("|"))
annotations_content = [
    "| Proximity | N/A | N/A | 19863 | 197 | 371 |",
    "| Brain |  psychEncode  | None |  17863 | 105 | 177 |",
    "| Brain-missing  | psychEncode | Missing gene | 22970 | 112 | 188 |" ,
    "| Brain-combined | psychEncode | Combined | 22970 | 236 | 378 |",
    "| Brain-mQTL |  psychEncode/mQTL  | None | 19113 | 131 | 202 |",
    "| Brain-mQTL-missing  | psychEncode/mQTL | Missing gene | 22974 | 127 | 191 |" ,
    "| Brain-mQTL-combined | psychEncode/mQTL | Combined | 22974  | 254 | 391 |",
    "| Blood  | eQTLgen/GTEx(Blood)  | None | 14375 | 661 | 888 |",
    "| Blood-missing | eQTLgen/GTEx(Blood)  | Missing gene | 21247 | 510 | 799 |",
    "| Blood-combined | eQTLgen/GTEx(Blood)  | Combined | 21247 | 575 | 835 |",
    "| AlleQTLs | GTEx(No brain)/psychEncode/eQTLgen | None | 20837 | 612 | 918 |",
    "| AlleQTLs-missing | GTEx(No brain)/psychEncode/eQTLgen | Missing gene | 23255 | 562 | 883 |",
    "| AlleQTLs-combined | GTEx(No brain)/psychEncode/eQTLgen | Combined | 23255 | 645 | 929 |",
    "| AlleQTLs-mQTL | GTEx(No brain)/psychEncode/eQTLgen/mQTL | None | 21194  | 617 | 920 |",
    "| AlleQTLs-mQTL-missing | GTEx(No brain)/psychEncode/eQTLgen/mQTL | Missing gene | 23257 | 572 | 891 |",
    "| AlleQTLs-mQTL-combined | GTEx(No brain)/psychEncode/eQTLgen/mQTL | Combined | 23257 | 655 | 936 |",
]

for item in annotations_content:
    annotations_df.loc[len(annotations_df)] = re.sub(r" ?\| ?", "|", item).split("|")

drop_cols = [0,7]
annotations_df.drop(annotations_df.columns[drop_cols], axis=1, inplace=True)
    

if formatting == "HTML":
    display(HTML(annotations_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            annotations_df.to_latex(
                column_format="p{3.5cm}p{5.9cm}p{1.7cm}p{0.8cm}p{1.1cm}p{1.1cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=(r"List of compiled annotation data sets, describing the functional data sources and "
                         r"basic statistics for each. \normalfont{Contributing QTL sources are "
                         r"given for each annotation along with the method used to combine this functional "
                         r"data with the proximity data, if applicable.}"
                        ),
                label="tab:db_annotations",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# ### Statistical analysis
#
# #### Comparing gene set analysis results
#
# A hierarchy of statistics was used to compare results from executions of the gene set analysis pipeline, as listed in table \ref{tab:metrics_hierarchy}. The principal metric employed was the number of drugs found to be significantly associated with the disease, as an analysis providing the greatest number of significantly associated drugs will be the most likely to identify novel repurposing opportunities. A significant association in this case was defined to be one where following the transformation of each drug's disease association p-value to a Storey q-value<citet data-cite="storey_direct_2002"><sup>storey</sup></citet>, the resulting drug-disease association q-value was less than 0.05.
#
# The second metric used to compare gene set analysis results was the number of treatment drugs that were found to be significantly associated with the disease, again selected by having a Storey q-value of less than 0.05 and appearing in the list of treatment drugs collated from the DrugBank website. This metric can be seen as a proxy measure of the capability of the executed pipeline to find novel drug repurposing opportunities, were the disease not to have any known treatment drugs.
#
# The final two metrics in the hierarchy are the Mann-Whitney and Fisher p-values. The Mann-Whitney non-parametric test as used here produces a p-value that measures the capability of the analysis to differentiate between the treatment drug and non-treatment drug groups without imposing a false positive rate threshold, with a lower p-value indicating a greater capability. The Fisher p-value similarly measures the difference in rates of significant enrichment between treatment and non-treatment drugs at a specified significance value, chosen to be 0.05. Both rely on the number of significant treatment drugs found, but give a broader view by incorporating the complementary measurements of non-significant treatment, significant non-treatment and non-significant non-treatment drug numbers.
#
# However, metrics involving the treatment/non-treatment classification of drugs are not as robust as the number of significant drugs found, as many drugs that could be used for treatment of the disease do not appear in the treatment drug list due to there being other more effective drugs or known side effects. Such drugs may well be correctly significantly associated with a disease by the analysis pipeline but marked as non-treatment due to factors not considered by the analysis. Alongside this, the gene-set analysis pipeline is not capable of differentiating between non-treatment drugs with side-effects that exacerbate the disease and treatment drugs that alleviate causes of the disease. Hence comparisons using treatment drug classifications are not an ideal indicator of therapeutic effect, but can still prove insightful when comparing the relative difference between results.

# +
metrics_df = pd.DataFrame(
    [
        [1, "Number of significant drugs"], [2, "Number of significant treatment drugs"],
        [3, "Mann-Whitney p-value - treatment vs. non-treatment"], [4, "Fisher p-value - treatment vs. non-treatment"],
    ],
    columns=["Order of comparison", "Metric"],
)
metrics_df.set_index("Order of comparison", inplace=True)

if formatting == "HTML":
    display(HTML(metrics_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            metrics_df.to_latex(
                column_format="p{3cm}p{8cm}",
                multirow=True,
                multicolumn=True,
                caption=r"Hierarchy of metrics obtained from gene-set analysis pipeline executions, used "
                        r"to compare and rank execution results.",
                label="tab:metrics_hierarchy",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# #### UpSet plots and significant drug enrichment across annotation results
#
# An UpSet plot<citet data-cite="lex_upset_2014"><sup>upset</sup></citet> shows the relationships between large numbers of sets in a matrix-like diagram, incorporating set sizes, set intersections and aggregates alongside summary statistics. For each disease, two UpSet plots were created to show the start and end of the list of annotation result drug set intersections, ordered by the number of annotation result sets that partake in the intersection, otherwise referred to as the degree of overlap. From these UpSet plots, the overlap of significant drugs can be studied to discover the scale of agreement between annotation results. 
#
# Further to this, the nature of the drugs most commonly found in annotation results was investigated by creating a table of drugs most frequently enriched in results from the set of 32 annotations, along with their treatment/non-treatment status and if they are known to affect blood pressure in the case of hypertension or cholesterol levels in the case of hypercholesterolemia. The drugs less frequently enriched in results were also studied to ascertain if significant treatment drugs are enriched in the set of annotation results.

# ## Hypertension results

# +
data_path = "../../../../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"

ht_gwas = "GA_I10"
run_id = "3qtl"
MAGMA_VER = 109

N_GENE_THRESH = 4
SIGNIF_THRESH = 0.05

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

# All GO genesets
go_genesets_df = pd.read_csv(
    data_path + "magma_analysis/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt",
    header=None,
    sep="\t",
    index_col=0,
)
# All GO + ATC genesets
genesets_df = pd.concat([atc_genesets_df, go_genesets_df])

# Emsembl to HGNC translation table
gene_hgnc_df = pd.read_csv(
    data_path + "wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    header=None,
    names=["Ensembl", "1", "chrom", "chromStart", "chromEnd", "HGNC"],
)[["Ensembl", "HGNC"]]

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

ht_dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            ht_gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        ht_dti_df_sets[dtis.index(dti)].append(newlist)
# -

# For the results involving DUGGIE as the source of drug-gene interaction data, the highest number of significant results, 52, was achieved with functional annotations composed of QTLs from all tissues, including mQTL data from the frontal cortex brain region, with positional data included using the 'combined' hybrid approach, compared to 33 drugs identified from the proximity analysis. However, the proximity analysis identified the greater number of treatment drugs, 18, in comparison to 14 significant treatment drugs from the combined hybrid annotation involving both brain tissue eQTLs and brain tissue eQTLs with frontal cortex mQTLs. The confusion matrix summary for all sixteen DUGGIE analyses can be found in table \ref{tab:ht_cmatrix_drugs_c4}. 
#
# \newpage

# +
ht_stats_df = pd.DataFrame()
dti = "duggie"
for annot in annots:
    (
        ht_dti_df_sets[dtis.index(dti)][annots.index(annot)],
        group_labels,
        colours,
    ) = rdisp.group_drugs(
        ht_dti_df_sets[dtis.index(dti)][annots.index(annot)], ht_gwas
    )
    ht_stats_df = ht_stats_df.append(
        rdisp.show_treatment_drug_stats(
            list([ht_dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, SIGNIF_THRESH
        )
    )
    
# Some formatting
ht_stats_df["Sens."] = ht_stats_df["Sens."].map("{:,.3f}".format)
ht_stats_df["Spec."] = ht_stats_df["Spec."].map("{:,.3f}".format)
ht_stats_df["FDR"] = ht_stats_df["FDR"].map("{:,.3f}".format)
ht_stats_df.drop(["DTI DB"], axis=1, inplace=True)
ht_stats_df['Annotation'] = display_annots

ht_cmatrix_stats_df = ht_stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(ht_cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_cmatrix_stats_df.to_latex(
                column_format="p{3.2cm}p{1.1cm}p{2cm}R{1.5cm}p{1.6cm}R{2cm}R{2cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=r"Confusion matrix results from multiple runs of a gene-set analysis pipeline differing "
                        r"in the gene annotation used, analysing hypertension using the DUGGIE drug-gene "
                        r"interaction dataset from q-value results at a 5\% significance level.",
                label="tab:ht_cmatrix_drugs_c4",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# This trend is repeated in table \ref{tab:ht_deriv_drugs_c4}, consisting of derived statistics from the data in table \ref{tab:ht_cmatrix_drugs_c4}, with both proximity and brain tissue-combined annotations showing the greatest ability to identify a difference between significant enrichment rates or treatment and non-treatment drugs as evidenced by the very low Fisher p-values.
#
# \newpage

# +
ht_deriv_stats_df = ht_stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

if formatting == "HTML":
    display(HTML(ht_deriv_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_deriv_stats_df.to_latex(
                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=r"Derived confusion matrix results from multiple runs of a gene-set analysis pipeline differing "
                        r"in the gene annotation used, analysing hypertension using the DUGGIE drug-gene "
                        r"interaction dataset from q-value results at a 5\% significance level.",
                label="tab:ht_deriv_drugs_c4",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# The Mann-Whitney and Fisher p-values also indicate that the gene set analysis pipeline identifies treatment drugs at a significantly higher rate than non-treatment drugs, giving an indication of the validity of the pipeline, with only three of the sixteen annotations not achieving Fisher p-values of less than 0.05. Two of the Mann-Whitney p-values, both for blood tissue, also did not produce a p-value less then 0.05, with the proximity annotation generating the smallest and highly significant Fisher and Mann-Whitney p-values of 3.17e-11 and 5.17e-16 respectively.

# ### Hypertension - drug-gene interaction dataset comparison
#
# To briefly check the conclusions of chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}, that DUGGIE has an advantage over STITCH in identifying larger numbers of significant drugs, a comparison between the number of significant drugs, significant treatment drugs, derived Fisher p-values and Mann-Whitney p-values for each drug-gene interaction dataset can be found in table \ref{tab:ht_dti_comparison}. This table arranges the results of the 32 gene set analysis pipeline executions, differing in annotation and drug-gene interaction dataset, such that this comparison can be made.
#
# As can be seen from this table, the maximum number of significant drugs discovered for each annotation (highlighted in green) is predominantly achieved using DUGGIE. DUGGIE also finds the highest number of significant drugs, 52, using the 'all eQTL with mQTL and combined hybrid approach' annotation, highlighted in red. DUGGIE also finds the greatest number of treatment drugs in more cases, 8, than STITCH, with 5 - but STITCH appears to have an advantage where a larger number of treatment drugs (i.e. over 10) are found to be significantly associated with hypertension, including for the proximity annotation where the greatest number of treatment drugs, 21, are discovered.

# Generate DUGGIE/STITCH comparison table
ht_stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            ht_dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            ht_dti_df_sets[dtis.index(dti)][annots.index(annot)], ht_gwas
        )
    
        ht_stats_df = ht_stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([ht_dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, SIGNIF_THRESH
            )
        )

# +
ht_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'DUGGIE Sig.', 'STITCH Sig.', 
                                        'DUGGIE Sig. Treatment' , 'STITCH Sig. Treatment',
                                        'DUGGIE Mann-Whitney p-value' , 'STITCH Mann-Whitney p-value',
                                        'DUGGIE Fisher p-value' , 'STITCH Fisher p-value',
                                       ])

ht_dgi_stats_df['DUGGIE Sig.'] = ht_stats_df[ht_stats_df['DTI DB'] == 'duggie']['Sig.']
ht_dgi_stats_df['DUGGIE Sig. Treatment'] = ht_stats_df[ht_stats_df['DTI DB'] == 'duggie']['Sig. Treatment']
ht_dgi_stats_df['DUGGIE Fisher p-value'] = ht_stats_df[ht_stats_df['DTI DB'] == 'duggie']['Fisher p-value'].astype(float)
ht_dgi_stats_df['DUGGIE Mann-Whitney p-value'] = ht_stats_df[ht_stats_df['DTI DB'] == 'duggie']['Mann-Whitney p-value'].astype(float)
ht_dgi_stats_df['STITCH Sig.'] = ht_stats_df[ht_stats_df['DTI DB'] == 'stitch']['Sig.']
ht_dgi_stats_df['STITCH Sig. Treatment'] = ht_stats_df[ht_stats_df['DTI DB'] == 'stitch']['Sig. Treatment']
ht_dgi_stats_df['STITCH Fisher p-value'] = ht_stats_df[ht_stats_df['DTI DB'] == 'stitch']['Fisher p-value'].astype(float)
ht_dgi_stats_df['STITCH Mann-Whitney p-value'] = ht_stats_df[ht_stats_df['DTI DB'] == 'stitch']['Mann-Whitney p-value'].astype(float)

ht_dgi_stats_df['Annotation'] = display_annots
ht_dgi_stats_df.set_index('Annotation', inplace=True)

# +
# Highlight max /min values from each pair of values
ht_dgi_stats_style = ht_dgi_stats_df.style.highlight_max(axis=1, subset=['DUGGIE Sig.', 'STITCH Sig.'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_max(axis=1, subset=['DUGGIE Sig. Treatment', 'STITCH Sig. Treatment'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_min(axis=1, subset=['DUGGIE Fisher p-value', 'STITCH Fisher p-value'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_min(axis=1, subset=['DUGGIE Mann-Whitney p-value', 'STITCH Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
ht_dgi_stats_style.highlight_max(axis=None, subset=['DUGGIE Sig.', 'STITCH Sig.'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_max(axis=None, subset=['DUGGIE Sig. Treatment', 'STITCH Sig. Treatment'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_min(axis=None, subset=['DUGGIE Fisher p-value', 'STITCH Fisher p-value'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_min(axis=None, subset=['DUGGIE Mann-Whitney p-value', 'STITCH Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = ht_dgi_stats_style.format("{:,.2e}", subset=['DUGGIE Fisher p-value', 
                                                     'STITCH Fisher p-value', 
                                                     'DUGGIE Mann-Whitney p-value', 
                                                     'STITCH Mann-Whitney p-value'])
# -

if formatting == "HTML":
    display(HTML(ht_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format="p{3.5cm}|p{1cm}p{1cm}|p{1cm}p{1cm}|p{1.2cm}p{1.2cm}|p{1.2cm}p{1.2cm}",
                caption=r"Table of gene set analysis results analysing hypertension from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"the drug-gene interaction data set DUGGIE against STITCH."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:ht_dti_comparison",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# \newpage

# The Fisher and Mann-Whitney p-values in table \ref{tab:ht_dti_comparison} also corroborate the conclusions made in chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia} - that both STITCH and DUGGIE are comparable and effective at detecting the difference between treatment and non-treatment hypertension drugs, as can be evidenced by the set of Fisher p-values, and DUGGIE is consistently but marginally better at differentiating between the treatment and non-treatment drug groups, as shown by the set of Mann-Whitney p-values.

# ### Hypertension - functional, proximal and hybrid annotation comparison
#
# To view the evidence of the effect of functionally annotating genes and augmenting functional annotations with positional data, the matrix of 32 annotation execution results were tabled such that proximity annotation results are compared against functional QTL annotation, the 'missing gene' hybrid and 'combined gene' hybrid annotation approaches. As there are only two proximity results to compare, one each for DUGGIE and STITCH, these are duplicated in a column against their respective annotation combinations involving that drug-gene interaction dataset for comparison. As this produces a larger table, the table is split into two - one of significant tallies, table \ref{tab:ht_annotation_comparison1}, and the other of derived values, table \ref{tab:ht_annotation_comparison2}.
#
# Across both tables, it can be seen that the proximity annotation is the successful scheme in most annotation categories - with 33 significantly associated drugs identified by both DUGGIE and STITCH it was the largest in 7 out of 10 annotation comparisons. However, notably the proximity annotation does not achieve the highest overall count of significant drugs (52), which is obtained using the combined hybrid approach utilising both eQTLs and mQTLs for all tissues and the DUGGIE drug-gene interaction dataset. When looking at the significant treatment drugs for each annotation category, the proximity annotation achieves the highest score in all cases, with only the combined hybrid functional approach coming close and this dominance of the proximity annotation is repeated across the Fisher and Mann-Whitney p-values with only one non-proximity top result, a Mann-Whitney p-value using the eQTL brain combined hybrid annotation and STITCH drug-gene interaction dataset.
#
# As the QTL functional annotation metrics, both absolute counts and p-values, are consistently below that of the proximity annotations, it cannot be said that this evidence supports the hypothesis that functionally annotating genes increases the enrichment of significant drugs, at least for hypertension. Despite this poorer result, the combined hybrid approach which augments the functional annotation with proximally derived data manages to achieve a large improvement in one of the cases, the combined hybrid approach utilising both eQTLs and mQTLs for all tissues and the DUGGIE drug-gene interaction dataset.
#
# This one case discovering the most significant drugs, however, does not lead on to provide compelling evidence that the hypothesis of augmenting functional annotations with positional data increases significant treatment drug enrichment, as it identifies less treatment drugs (10) than the corresponding non-functional pipeline execution with proximally annotated genes (18).

# +
ht_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'Prox Sig.', 'QTL Sig.', 'Missing Sig.', 'Combined Sig.', 
                                        'Prox Sig. Treatment' , 'QTL Sig. Treatment', 'Missing Sig. Treatment',
                                        'Combined Sig. Treatment',
                                       ])

display_annots3 = [
    "eQTL brain DUGGIE",
    "e/mQTL brain DUGGIE",
    "eQTL blood DUGGIE",
    "eQTL all DUGGIE",
    "e/mQTL all DUGGIE",
    "eQTL brain STITCH",
    "e/mQTL brain STITCH",
    "eQTL blood STITCH",
    "eQTL all STITCH",
    "e/mQTL all STITCH",
]

temp = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('prox')]['Sig.'])
ht_dgi_stats_df['Prox Sig.'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
ht_dgi_stats_df['QTL Sig.'] = list(ht_stats_df[~ht_stats_df['Annotation'].str.contains('hybrid|prox')]['Sig.'])
ht_dgi_stats_df['Missing Sig.'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('hybrid$')]['Sig.'])
ht_dgi_stats_df['Combined Sig.'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('both')]['Sig.'])

temp = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('prox')]['Sig. Treatment'])
ht_dgi_stats_df['Prox Sig. Treatment'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
ht_dgi_stats_df['QTL Sig. Treatment'] = list(ht_stats_df[~ht_stats_df['Annotation'].str.contains('hybrid|prox')]['Sig. Treatment'])
ht_dgi_stats_df['Missing Sig. Treatment'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('hybrid$')]['Sig. Treatment'])
ht_dgi_stats_df['Combined Sig. Treatment'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('both')]['Sig. Treatment'])

ht_dgi_stats_df['Annotation'] = display_annots3
ht_dgi_stats_df.set_index('Annotation', inplace=True)

# +
# Highlight max /min values from each pair of values
ht_dgi_stats_style = ht_dgi_stats_df.style.highlight_max(axis=1, subset=['Prox Sig.', 'QTL Sig.',
                                                                         'Missing Sig.', 'Combined Sig.'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_max(axis=1, subset=['Prox Sig. Treatment',  'QTL Sig. Treatment',
                                                 'Missing Sig. Treatment', 'Combined Sig. Treatment'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
ht_dgi_stats_style.highlight_max(axis=None, subset=['Prox Sig.', 'QTL Sig.',
                                                    'Missing Sig.', 'Combined Sig.'],
                                                        props='color:{Red};')
dummy=ht_dgi_stats_style.highlight_max(axis=None, subset=['Prox Sig. Treatment',  'QTL Sig. Treatment',
                                                 'Missing Sig. Treatment', 'Combined Sig. Treatment'],
                                                        props='color:{Red};')
# -

if formatting == "HTML":
    display(HTML(ht_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format=r"p{3.2cm}|p{0.6cm}p{0.6cm}p{0.8cm}p{1.5cm}|p{0.6cm}p{0.6cm}"
                              r"p{0.8cm}p{1.5cm}",
                caption=r"Table of basic gene set analysis results analysing hypertension from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"functional, proximal and hybrid annotations."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:ht_annotation_comparison1",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
ht_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'Prox Mann-Whitney p-value' , 'QTL Mann-Whitney p-value',
                                        'Missing Mann-Whitney p-value', 'Combined Mann-Whitney p-value',
                                        'Prox Fisher p-value', 'QTL Fisher p-value',
                                        'Missing Fisher p-value', 'Combined Fisher p-value', 
                                       ])

display_annots3 = [
    "eQTL brain DUGGIE",
    "e/mQTL brain DUGGIE",
    "eQTL blood DUGGIE",
    "eQTL all DUGGIE",
    "e/mQTL all DUGGIE",
    "eQTL brain STITCH",
    "e/mQTL brain STITCH",
    "eQTL blood STITCH",
    "eQTL all STITCH",
    "e/mQTL all STITCH",
]

temp = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('prox')]['Fisher p-value'].astype(float))
ht_dgi_stats_df['Prox Fisher p-value'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
ht_dgi_stats_df['QTL Fisher p-value'] = list(ht_stats_df[~ht_stats_df[
                                                'Annotation'].str.contains('hybrid|prox')]['Fisher p-value'].astype(float))
ht_dgi_stats_df['Missing Fisher p-value'] = list(ht_stats_df[ht_stats_df[
                                                'Annotation'].str.contains('hybrid$')]['Fisher p-value'].astype(float))
ht_dgi_stats_df['Combined Fisher p-value'] = list(ht_stats_df[ht_stats_df[
                                                'Annotation'].str.contains('both')]['Fisher p-value'].astype(float))

temp = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('prox')]['Mann-Whitney p-value'].astype(float))
ht_dgi_stats_df['Prox Mann-Whitney p-value'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
ht_dgi_stats_df['QTL Mann-Whitney p-value'] = list(ht_stats_df[~ht_stats_df[
                                                'Annotation'].str.contains('hybrid|prox')]['Mann-Whitney p-value'].astype(float))
ht_dgi_stats_df['Missing Mann-Whitney p-value'] = list(ht_stats_df[ht_stats_df[
                                                'Annotation'].str.contains('hybrid$')]['Mann-Whitney p-value'].astype(float))
ht_dgi_stats_df['Combined Mann-Whitney p-value'] = list(ht_stats_df[ht_stats_df[
                                                'Annotation'].str.contains('both')]['Mann-Whitney p-value'].astype(float))

ht_dgi_stats_df['Annotation'] = display_annots3
ht_dgi_stats_df.set_index('Annotation', inplace=True)

# +
ht_dgi_stats_style = ht_dgi_stats_df.style

# Highlight max /min values from each pair of values
ht_dgi_stats_style.highlight_min(axis=1, subset=['Prox Fisher p-value', 'QTL Fisher p-value', 
                                                 'Missing Fisher p-value',  'Combined Fisher p-value'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_min(axis=1, subset=['Prox Mann-Whitney p-value', 'QTL Mann-Whitney p-value', 
                                                 'Missing Mann-Whitney p-value', 'Combined Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
ht_dgi_stats_style.highlight_min(axis=None, subset=['Prox Fisher p-value', 'QTL Fisher p-value', 
                                                 'Missing Fisher p-value',  'Combined Fisher p-value'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_min(axis=None, subset=['Prox Mann-Whitney p-value', 'QTL Mann-Whitney p-value', 
                                                 'Missing Mann-Whitney p-value', 'Combined Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = ht_dgi_stats_style.format("{:,.2e}", subset=['Prox Fisher p-value', 
                                                     'QTL Fisher p-value', 
                                                     'Missing Fisher p-value', 
                                                     'Combined Fisher p-value', 
                                                     'Prox Mann-Whitney p-value', 
                                                     'QTL Mann-Whitney p-value', 
                                                     'Missing Mann-Whitney p-value',
                                                     'Combined Mann-Whitney p-value'])
# -

# TODO - just derived values here (as table is too big otherwise)
if formatting == "HTML":
    display(HTML(ht_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format=r"p{3.2cm}|p{1.15cm}p{1.15cm}p{1.15cm}p{1.5cm}|p{1.15cm}p{1.15cm}"
                              r"p{1.15cm}p{1.5cm}",
                caption=r"Table of derived gene set analysis results analysing hypertension from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"functional, proximal and hybrid annotations."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:ht_annotation_comparison2",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Hypertension - functional annotation tissue sample size comparison
#
# The results from the 32 analyses were again re-arranged in table \ref{tab:ht_tissue_comparison}, this time to compare results from functional annotation of different tissue sample sizes, varying over the hybrid approach used, mQTL addition and drug-gene interaction dataset used. This allows a comparison between blood tissue, brain tissue and all tissue results that differ in the sample sizes used to generate functional QTLs. The blood tissue annotation is mainly sourced from eQTLGen QTLs, which have the highest sample size of 31,684 samples, with a smaller contribution from GTEx blood QTLs obtained from 838 samples. All QTLs contributing to the brain tissue annotation come from PsychENCODE, obtained from 1,387 samples. The brain annotation including mQTL data also includes these PsychENCODE QTLs, plus mQTL data obtained from 134 tissue samples. Finally the majority of all tissue data is composed of GTEx QTLs from 48 tissue types, taken from 838 donors whilst also incorporating the eQTLGen and PsychENCODE datasets. Again the all tissue with mQTL data annotation has in addition mQTL data from 134 samples. To summarise, the blood tissue annotations are majority sourced using the largest sample size, followed by brain tissue annotation with a smaller sample size and finally the all tissue annotations mostly consist of QTL data sourced from the smallest sample size.
#
# No combination of blood results with mQTL data are available, as the mQTL data used are specific to brain tissue - these are marked as N/A.
#
# Of these tissue types, blood identified the largest number of significant drugs in half of the six annotation types for which blood data was available, however the differences were modest compared to the other three instances, where brain or all tissues found higher numbers of significant drugs differing sometimes by a large multiplying factor compared to blood tissue results. The all tissue annotation using both eQTL and mQTL data with the combined hybrid approach also produced the highest tally with 52 significantly associated drugs. The lack of success by blood tissue annotations could be because the blood tissue analyses have the lowest number of annotated genes identified by functional annotation (at 14,375, from table \ref{tab:db_annotations}), reducing the power of the functional analysis, but the missing/combined hybrid annotation approaches which involve comparable numbers of annotated genes between the three tissue types also show poor numbers of significant drugs for blood tissue based annotations.
#
# Brain tissue annotations were found to identify the largest number of treatment drugs used to treat hypertension, both identifying the largest number of drugs for any annotation (19) and identifying more drugs than other tissues for the largest number of annotations (7). This could be due to brain tissue containing many of the functional genes involved in hypertension, with the greater number of genes functionally annotated under the all-tissue analyses adding statistical noise, especially as many of the GTEx eQTL samples are quite small compared to PsychENCODE which could point to power issues<citet data-cite="button_power_2013"><sup>button</sup></citet>. The Fisher and Mann-Whitney p-values mirror this greater success of brain tissue annotations, with the all-tissue annotations suffering from identifying a high number of significant non-treatment drugs.
#
# This evidence does not support the hypothesis that using QTL data obtained from larger numbers of tissue samples increases the number of significant drugs and treatment drugs, but it is also likely that statistical noise brought by annotating more genes based on relatively small eQTL datasets is a confounding factor.

# +
ht_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'Blood Sig.', 'Brain Sig.', 'All Sig.', 
                                        'Blood Sig. Treatment' , 'Brain Sig. Treatment', 'All Sig. Treatment',
                                        'Blood Mann-Whitney p-value' , 'Brain Mann-Whitney p-value',
                                        'All Mann-Whitney p-value',
                                        'Blood Fisher p-value' , 'Brain Fisher p-value', 'All Fisher p-value',
                                       ])

three_na = pd.Series([np.NaN, np.NaN, np.NaN])

display_annots2 = [
    "eQTL DUGGIE",
    "eQTL-miss. DUGGIE",
    "eQTL-comb. DUGGIE",
    "e/mQTL DUGGIE",
    "e/mQTL-miss. DUGGIE",
    "e/mQTL-comb. DUGGIE",
    "eQTL STITCH",
    "eQTL-miss. STITCH",
    "eQTL-comb. STITCH",
    "e/mQTL STITCH",
    "e/mQTL-miss. STITCH",
    "e/mQTL-comb. STITCH",
]

temp = ht_stats_df[ht_stats_df['Annotation'].str.contains('blood')]['Sig.']
ht_dgi_stats_df['Blood Sig.'] = list(temp[0:3].append(three_na).append(temp[3:6]).append(three_na))
ht_dgi_stats_df['Brain Sig.'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('brain')]['Sig.'])
ht_dgi_stats_df['All Sig.'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains('tissues')]['Sig.'])

temp = ht_stats_df[ht_stats_df['Annotation'].str.contains('blood')]['Sig. Treatment']
ht_dgi_stats_df['Blood Sig. Treatment'] = list(temp[0:3].astype(int).append(three_na).append(temp[3:6]).append(three_na))
ht_dgi_stats_df['Brain Sig. Treatment'] = list(ht_stats_df[ht_stats_df[
                                            'Annotation'].str.contains('brain')]['Sig. Treatment'])
ht_dgi_stats_df['All Sig. Treatment'] = list(ht_stats_df[ht_stats_df[
                                            'Annotation'].str.contains('tissues')]['Sig. Treatment'])

temp = ht_stats_df[ht_stats_df['Annotation'].str.contains('blood')]['Fisher p-value'].astype(float)
ht_dgi_stats_df['Blood Fisher p-value'] = list(temp[0:3].append(three_na).append(temp[3:6]).append(three_na))
ht_dgi_stats_df['Brain Fisher p-value'] = list(ht_stats_df[ht_stats_df[
                                            'Annotation'].str.contains('brain')]['Fisher p-value'].astype(float))
ht_dgi_stats_df['All Fisher p-value'] = list(ht_stats_df[ht_stats_df[
                                            'Annotation'].str.contains('tissues')]['Fisher p-value'].astype(float))

temp = ht_stats_df[ht_stats_df['Annotation'].str.contains('blood')]['Mann-Whitney p-value'].astype(float)
ht_dgi_stats_df['Blood Mann-Whitney p-value'] = list(temp[0:3].append(three_na).append(temp[3:6]).append(three_na))
ht_dgi_stats_df['Brain Mann-Whitney p-value'] = list(ht_stats_df[ht_stats_df[
                                            'Annotation'].str.contains('brain')]['Mann-Whitney p-value'].astype(float))
ht_dgi_stats_df['All Mann-Whitney p-value'] = list(ht_stats_df[ht_stats_df[
                                            'Annotation'].str.contains('tissues')]['Mann-Whitney p-value'].astype(float))

ht_dgi_stats_df['Annotation'] = display_annots2
ht_dgi_stats_df.set_index('Annotation', inplace=True)

# +
ht_dgi_stats_style = ht_dgi_stats_df.style

# Highlight max /min values from each pair of values
ht_dgi_stats_style = ht_dgi_stats_df.style.highlight_max(axis=1, subset=['Blood Sig.', 'Brain Sig.', 'All Sig.'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_max(axis=1, subset=['Blood Sig. Treatment', 
                                                 'Brain Sig. Treatment', 'All Sig. Treatment'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_min(axis=1, subset=['Blood Fisher p-value',
                                                 'Brain Fisher p-value', 'All Fisher p-value'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_min(axis=1, subset=['Blood Mann-Whitney p-value', 
                                                 'Brain Mann-Whitney p-value', 'All Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
ht_dgi_stats_style.highlight_max(axis=None, subset=['Blood Sig.', 'Brain Sig.', 'All Sig.'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_max(axis=None, subset=['Blood Sig. Treatment', 
                                                    'Brain Sig. Treatment', 'All Sig. Treatment'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_min(axis=None, subset=['Blood Fisher p-value',
                                                    'Brain Fisher p-value', 'All Fisher p-value'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_min(axis=None, subset=['Blood Mann-Whitney p-value', 
                                                    'Brain Mann-Whitney p-value', 'All Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = ht_dgi_stats_style.format("{:,.2e}", subset=['Blood Fisher p-value', 
                                                     'Brain Fisher p-value', 
                                                     'All Fisher p-value', 
                                                     'Blood Mann-Whitney p-value', 
                                                     'Brain Mann-Whitney p-value', 
                                                     'All Mann-Whitney p-value'], na_rep = 'N/A')
# format ints
dummy = ht_dgi_stats_style.format("{:.0f}", subset=['Blood Sig.', 
                                                    'Blood Sig. Treatment',
                                                   ], na_rep = 'N/A')
# -

# \newpage

if formatting == "HTML":
    display(HTML(ht_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\tiny"))
    display(
        Latex(
            ht_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format=r"l|p{0.35cm}p{0.25cm}p{0.25cm}|p{0.4cm}p{0.4cm}p{0.55cm}|"
                              r"p{1.01cm}p{1.01cm}p{1.01cm}|p{0.95cm}p{0.95cm}p{0.95cm}",
                caption=r"Table of basic gene set analysis results analysing hypertension from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"tissue specific functional annotations."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:ht_tissue_comparison",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Hypertension - contribution of mQTL data to the enrichment of drugs
#
# The final hypothesis investigated, that the specific brain tissue mQTL data alongside eQTL data increases the enrichment of treatment drugs, was considered by once more displaying the matrix of results in a suitable table, table \ref{tab:ht_mqtl_comparison} - in this case by constructing two columns of 'with mQTL' and 'without mQTL' results for each of the four quantities of significant drugs, significant treatment drugs, Fisher and Mann-Whitney p-values. No blood tissue results are included as these do not involve mQTL data.
#
# By classifying the results in this way, it can be seen that including mQTL functional data does improve the number of significant drugs identified in most annotations - 8, compared to only 4 where not including mQTL data identified more. The top count of significant drugs, 52, was also achieved by including the mQTL functional dataset.
#
# This trend also extends to the number of significant treatment drugs identified as well as the Fisher p-values comparing the rates of significant treatment versus non-treatment drugs, although the highest significant treatment drug score across all the annotations, 19, is shared between with-mQTL and without-mQTL annotations. However, across all Mann-Whitney p-values the no-mQTL annotations have more desirable smaller values in most cases apart from only two mQTL cases (both occuring from combined scheme annotations), possibly because the no-mQTL analyses are ranking treatment drugs by p-value consistently higher than non-treatment drugs than the with-mQTL analyses, even among drugs that do not reach the formal cut-off for significance (FDR < 0.05).
#
# There was no clear evidence found here to support or reject the hypothesis that using mQTL data results in a greater enrichment of significant drugs and significant treatment drugs.

# +
ht_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'mQTL Sig.', 'no mQTL Sig.', 
                                        'mQTL Sig. Treatment' , 'no mQTL Sig. Treatment',
                                        'mQTL Mann-Whitney p-value' , 'no mQTL Mann-Whitney p-value',
                                        'mQTL Fisher p-value' , 'no mQTL Fisher p-value',
                                       ])

display_annots4 = [
    "Brain DUGGIE",
    "Brain-miss. DUGGIE",
    "Brain-comb. DUGGIE",
    "All DUGGIE",
    "All-miss. DUGGIE",
    "All-comb. DUGGIE",
    "Brain STITCH",
    "Brain-miss. STITCH",
    "Brain-comb. STITCH",
    "All STITCH",
    "All-miss. STITCH",
    "All-comb. STITCH",
]

ht_dgi_stats_df['mQTL Sig.'] = ht_stats_df[~ht_stats_df['Annotation'].str.contains("eqtl|blood|prox")]['Sig.']
ht_dgi_stats_df['mQTL Sig. Treatment'] = ht_stats_df[~ht_stats_df['Annotation'].str.contains(
                                                    "eqtl|blood|prox")]['Sig. Treatment']
ht_dgi_stats_df['mQTL Fisher p-value'] = ht_stats_df[~ht_stats_df['Annotation'].str.contains(
                                                    "eqtl|blood|prox")]['Fisher p-value'].astype(float)
ht_dgi_stats_df['mQTL Mann-Whitney p-value'] = ht_stats_df[~ht_stats_df['Annotation'].str.contains(
                                                    "eqtl|blood|prox")]['Mann-Whitney p-value'].astype(float)

ht_dgi_stats_df['no mQTL Sig.'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~ht_stats_df['Annotation'].str.contains("blood")]['Sig.'])
ht_dgi_stats_df['no mQTL Sig. Treatment'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~ht_stats_df['Annotation'].str.contains("blood")]['Sig. Treatment'])
ht_dgi_stats_df['no mQTL Fisher p-value'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~ht_stats_df['Annotation'].str.contains("blood")]['Fisher p-value'].astype(float))
ht_dgi_stats_df['no mQTL Mann-Whitney p-value'] = list(ht_stats_df[ht_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~ht_stats_df['Annotation'].str.contains("blood")]['Mann-Whitney p-value'].astype(float))

ht_dgi_stats_df['Annotation'] = display_annots4
ht_dgi_stats_df.set_index('Annotation', inplace=True)

# +
# Highlight max /min values from each pair of values
ht_dgi_stats_style = ht_dgi_stats_df.style.highlight_max(axis=1, subset=['mQTL Sig.', 'no mQTL Sig.'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_max(axis=1, subset=['mQTL Sig. Treatment', 'no mQTL Sig. Treatment'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_min(axis=1, subset=['mQTL Fisher p-value', 'no mQTL Fisher p-value'],
                                                        props='color:{Green};')
ht_dgi_stats_style.highlight_min(axis=1, subset=['mQTL Mann-Whitney p-value', 'no mQTL Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
ht_dgi_stats_style.highlight_max(axis=None, subset=['mQTL Sig.', 'no mQTL Sig.'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_max(axis=None, subset=['mQTL Sig. Treatment', 'no mQTL Sig. Treatment'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_min(axis=None, subset=['mQTL Fisher p-value', 'no mQTL Fisher p-value'],
                                                        props='color:{Red};')
ht_dgi_stats_style.highlight_min(axis=None, subset=['mQTL Mann-Whitney p-value', 'no mQTL Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = ht_dgi_stats_style.format("{:,.2e}", subset=['mQTL Fisher p-value', 
                                                     'no mQTL Fisher p-value', 
                                                     'mQTL Mann-Whitney p-value', 
                                                     'no mQTL Mann-Whitney p-value'])
# -

if formatting == "HTML":
    display(HTML(ht_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            ht_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format="p{3.5cm}|p{1cm}p{1cm}|p{1cm}p{1cm}|p{1.2cm}p{1.2cm}|p{1.2cm}p{1.2cm}",
                caption=r"Table of basic gene set analysis results analysing hypertension from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"annotations including, and not including mQTL data."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:ht_mqtl_comparison",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Hypertension - overlap of significant drugs found per annotation
#
# UpSet plots showing the overlap of significantly identified drugs between gene set analysis pipeline results for each of the 30 annotations that found at least one significant drug can be seen in figures \ref{fig:ht_upset_0_3} and \ref{fig:ht_upset_4_plus}, split for ease of display. Figure \ref{fig:ht_upset_0_3} shows the overlaps of 3 degrees or less, i.e. drugs identified as significant by three or fewer annotations, where only 10 analysis result drug sets have drugs unique to themselves, composed mainly of result drug sets from analyses that identify the greatest number of significant drugs overall. The largest set of 52 significantly enriched drugs, resulting from employing DUGGIE with an annotation using all eQTLs plus mQTLs and the combined proximity annotation scheme, contains 20 unique drugs not found by any other annotation, which is consistent with the difference in set size to the next largest set, illustrating that the higher powered annotations are uncovering new information.

# +
ATC_LEVELS = rdisp.get_atc_levels(data_path)

# We want all results, so we set gene count and q-value theshold to 1.
for dti in dtis:
    for annot in annots:
        rdisp.summarise_drug_results(
            ht_gwas,
            annot,
            ATC_LEVELS,
            magma_results_dir,
            summary_results_dir,
            dti,
            run_id,
            MAGMA_VER,
            1,
            1,
        )

# +
drug_sets_df = pd.DataFrame()

for dti in dtis:
    file_postfix = dti + "_qvals.tsv"
    for annot in annots:
        label = dti.upper() + " " + display_annots[annots.index(annot)]
        pval_file = [
            file
            for file in os.listdir(summary_results_dir)
            if file.endswith(file_postfix)
            and "found-" + ht_gwas + "-" + run_id + "-" + annot + "-" in file
        ]
        table = pd.read_csv(os.path.join(summary_results_dir, pval_file[0]), sep="\t")
        table2 = table[table["Q " + annot + " " + dti] < SIGNIF_THRESH]
        drug_set = table2["ATC_CODE"]
        # Add series as column to get a matrix of boolean columns
        drug_sets_df[label] = False
        for drug in drug_set:
            if drug not in drug_sets_df.index:
                drug_sets_df = drug_sets_df.append(pd.Series([False], name=drug))
            drug_sets_df.at[drug, label] = True
        
drug_sets_df = drug_sets_df.fillna(False)
drug_sets_df = drug_sets_df.drop(0, axis=1)
# -

bold_caption = (
    r"UpSet plot showing the overlap of drug sets identified to be significantly associated with hypertension "
    r" from multiple runs of a gene-set analysis pipeline differing in the gene annotation used, for up to 3 "
    r" degrees of overlap."
)
normal_caption = (
    r"Vertical bars show the gene set details per named annotation, with a histogram showing the significant"
    r" gene tallies per annotation. The horizontal histogram shows the size of overlap between the sets indicated"
    " horizontally in the dot matrix. Single dots represent the tally of genes unique to the named annotation gene"
    r" set."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=3)

rdisp.end_caption(formatting, r"fig:ht_upset_0_3")

# Looking at figure \ref{fig:ht_upset_4_plus} of higher degree overlaps, particularly the bottom left of the matrix which has a high density of dots, it can be seen that the drugs appearing as enriched under the greatest number of annotations also appear heavily in the sets with the highest tallies, with no grouping of drugs appearing to be specific to analyses with lower numbers of significantly associated drugs. This implies that generally, the annotations achieving higher rates of significant enrichment are supplementing and not replacing drugs found by annotations with lower rates of significant enrichment and that there is a degree of consensus between all annotations.

bold_caption = (
    r"UpSet plot showing the overlap of drug sets identified to be significantly associated with hypertension "
    r"from multiple runs of a gene-set analysis pipeline differing in the gene annotation used, "
    r"for 4 or more degrees of overlap."
)
normal_caption = (
    r"Vertical bars show the gene set details per named annotation, with a histogram showing the significant"
    r" gene tallies per annotation. The horizontal histogram shows the size of overlap between the sets indicated"
    " horizontally in the dot matrix. Single dots represent the tally of genes unique to the named annotation gene"
    r" set."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

usp_figure = plt.figure(dpi=300)
dummy = usp.plot(usp.from_indicators(drug_sets_df), 
                 subset_size='count', 
                 fig=usp_figure,
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=4,
                 max_degree=100)

rdisp.end_caption(formatting, r"fig:ht_upset_4_plus")
if formatting == "LaTeX":
    display(Latex(r"\newpage"))

# ### Hypertension - drug enrichment across annotation results
#
# To briefly explore the drugs that are consistently enriched across many annotations, those with the highest frequency of enrichment in analysis results across the 32 annotations are listed in table \ref{tab:ht_most_frequent_drugs}, by ATC code and drug name as well as if they are a known treatment drug for hypertension and if the drug is known to affect blood pressure. Those drugs appearing in seven or more analysis results are listed, giving 16 drugs in total with the most frequent, aliskiren, appearing in 24 of the 32 analyses. From this list of 16, 10 drugs are classed as treatment drugs by DrugBank, whilst all but one, tiazofurine, of the remaining 6 drugs are known to affect blood pressure. As tiazofurine is an experimental drug which DrugBank states may have a clinical use in cancer treatment<citet data-cite="noauthor_tiazofurine_nodate"><sup>tiazofurine</sup></citet>, not much information is available on it. However, DrugBank also states that it does interact with other drugs to increase the risk of thrombosis, which may indicate a link to hypertension<citet data-cite="huang_association_2016"><sup>huang_assoc</sup></citet>. 
#
# Conversely, of the 48 drugs that are unique to only one annotation result (not shown in the table), the majority are non-treatment drugs with only four treatment drugs identified - verapamil, captopril and two combination drugs, delapril with manidipine and rosuvastatin with amlodipine and perindopril. 

# \newpage

top_drugs_df = pd.DataFrame(drug_sets_df.transpose().sum(), columns=['Analysis count'])
top_drugs_df = top_drugs_df[top_drugs_df['Analysis count'] > 6]
top_drugs_df = pd.merge(top_drugs_df, ATC_LEVELS, left_index=True, right_index=True)
top_drugs_df.sort_values(by=['Analysis count'], ascending=False, inplace=True)
top_drugs_df = rdisp.add_treatment_state(ht_gwas, top_drugs_df)
top_drugs_df.index.set_names(["ATC code"], inplace=True)

if formatting == "HTML":
    display(HTML(top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            top_drugs_df.to_latex(
                column_format=r"p{3.2cm}p{1.5cm}p{5cm}p{1.5cm}p{1.5cm}",
                caption=r"Table of drugs most frequently enriched across annotations. ATC drug codes are listed "
                        r"alongside the count of analyses in which the drug is significantly enriched, if the drug "
                        r"is used to treat hypertension and if the drug is known to affect blood pressure.",
                label="tab:ht_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ## Hypercholesterolemia results
#
# Similarly to hypertension, the hypercholesterolemia results involving DUGGIE as the source of drug-gene interaction data found the highest number of significant drugs, 129, was achieved with functional annotations composed of QTLs from all tissues, this time without the addition of mQTL data from the frontal cortex brain region, but with positional data included using the 'combined' hybrid approach. This is far greater than the 26 significant drugs identified from the proximity analysis. However just over 10% of the 129 significantly identified drugs were treatment drugs, so the number of significantly enriched non-treatment drugs is high. The confusion matrix summary for all sixteen DUGGIE analyses can be found in table \ref{tab:hc_cmatrix_drugs}. 
#
# The statistics from these analyses are shown in table \ref{tab:hc_deriv_drugs}. The functional-only brain tissue annotation shows the best ability to identify treatment drugs as significantly enriched compared to non-treatment drugs, as evidenced by the very low Fisher p-value of 2.37e-20. There is general evidence for the enrichment of treatment drugs among the significant drugs across the annotations, with only one Fisher p-value and no Mann-Whitney p-values scoring higher than 0.05, the chosen threshold for significance.

# +
data_path = "../../../../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"

hc_gwas = "GA_E78"

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

# All GO genesets
go_genesets_df = pd.read_csv(
    data_path + "magma_analysis/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt",
    header=None,
    sep="\t",
    index_col=0,
)
# All GO + ATC genesets
genesets_df = pd.concat([atc_genesets_df, go_genesets_df])

# Emsembl to HGNC translation table
gene_hgnc_df = pd.read_csv(
    data_path + "wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    header=None,
    names=["Ensembl", "1", "chrom", "chromStart", "chromEnd", "HGNC"],
)[["Ensembl", "HGNC"]]

hc_dti_df_sets = [[], []]

for dti in dtis:
    for annot in annots:
        newlist = rdisp.read_results_files(
            hc_gwas, annot, magma_results_dir, run_id, dti, MAGMA_VER
        )
        newlist = newlist[newlist["NGENES"] >= N_GENE_THRESH]
        hc_dti_df_sets[dtis.index(dti)].append(newlist)

# +
hc_stats_df = pd.DataFrame()
dti = "duggie"
for annot in annots:
    (
        hc_dti_df_sets[dtis.index(dti)][annots.index(annot)],
        group_labels,
        colours,
    ) = rdisp.group_drugs(
        hc_dti_df_sets[dtis.index(dti)][annots.index(annot)], hc_gwas
    )
    hc_stats_df = hc_stats_df.append(
        rdisp.show_treatment_drug_stats(
            list([hc_dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, SIGNIF_THRESH
        )
    )
    
# Some formatting
hc_stats_df["Sens."] = hc_stats_df["Sens."].map("{:,.3f}".format)
hc_stats_df["Spec."] = hc_stats_df["Spec."].map("{:,.3f}".format)
hc_stats_df["FDR"] = hc_stats_df["FDR"].map("{:,.3f}".format)
hc_stats_df.drop(["DTI DB"], axis=1, inplace=True)
hc_stats_df['Annotation'] = display_annots

hc_cmatrix_stats_df = hc_stats_df.iloc[:, 0:7]

if formatting == "HTML":
    display(HTML(hc_cmatrix_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_cmatrix_stats_df.to_latex(
                column_format="p{3.2cm}p{1.1cm}p{2cm}R{1.5cm}p{1.6cm}R{2cm}R{2cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=r"Confusion matrix results from multiple runs of a gene-set analysis pipeline differing "
                        r"in the gene annotation used, analysing hypercholesterolemia using the DUGGIE drug-gene "
                        r"interaction dataset from q-value results at a 5\% significance level.",
                label="tab:hc_cmatrix_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# \newpage

# +
hc_deriv_stats_df = hc_stats_df.iloc[:, [0, 7, 8, 9, 10, 11]]

if formatting == "HTML":
    display(HTML(hc_deriv_stats_df.to_html(index=False)))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_deriv_stats_df.to_latex(
                column_format="p{3cm}p{2cm}p{2cm}p{2cm}p{2cm}p{2cm}",
                multirow=True,
                multicolumn=True,
                index=False,
                caption=r"Derived confusion matrix results from multiple runs of a gene-set analysis "
                        r"pipeline differing in the gene annotation used, analysing hypercholesterolemia using "
                        r"the DUGGIE drug-gene interaction dataset from q-value results at a 5\% "
                        r"significance level.",
                label="tab:hc_deriv_drugs",
                position="htbp",
                escape=False,
            )
        )
    )
    display(Latex(r"\normalsize"))
# -

# ### Hypercholesterolemia - drug-gene interaction dataset comparison
#
# Table \ref{tab:hc_dti_comparison} arranges the results of the 32 gene set analysis pipeline executions for hypercholesteremia to briefly return to the discussion of chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}, which concluded that DUGGIE has an advantage over STITCH in identifying significant drugs.
#
# Again, DUGGIE has an advantage over STITCH in finding significantly associated drugs with almost all annotations. It also achieves the top results for the largest number of significant drugs and significant treatment drugs as well as in obtaining the lowest Fisher and Mann-Whitney p-values overall. DUGGIE fares comparably with STITCH in differentiating between treatment and non-treatment drugs as given by the Mann-Whitney p-values, as each has the top result in eight annotation categories. But as DUGGIE has a high number of significant non-treatment drugs at the stringent 0.05 significance threshold, which is not necessarily undesirable as this may highlight novel treatment drugs. However, the majority of annotations give a better Fisher p-value for STITCH although this may change if a different level of significance is chosen.

# Generate DUGGIE/STITCH comparison table
hc_stats_df = pd.DataFrame()
for dti in dtis:
    for annot in annots:
        (
            hc_dti_df_sets[dtis.index(dti)][annots.index(annot)],
            group_labels,
            colours,
        ) = rdisp.group_drugs(
            hc_dti_df_sets[dtis.index(dti)][annots.index(annot)], hc_gwas
        )
    
        hc_stats_df = hc_stats_df.append(
            rdisp.show_treatment_drug_stats(
                list([hc_dti_df_sets[dtis.index(dti)][annots.index(annot)]]), 3, SIGNIF_THRESH
            )
        )

# +
hc_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'DUGGIE Sig.', 'STITCH Sig.', 
                                        'DUGGIE Sig. Treatment' , 'STITCH Sig. Treatment',
                                        'DUGGIE Mann-Whitney p-value' , 'STITCH Mann-Whitney p-value',
                                        'DUGGIE Fisher p-value' , 'STITCH Fisher p-value',
                                       ])

hc_dgi_stats_df['DUGGIE Sig.'] = hc_stats_df[hc_stats_df['DTI DB'] == 'duggie']['Sig.']
hc_dgi_stats_df['DUGGIE Sig. Treatment'] = hc_stats_df[hc_stats_df['DTI DB'] == 'duggie']['Sig. Treatment']
hc_dgi_stats_df['DUGGIE Fisher p-value'] = hc_stats_df[hc_stats_df['DTI DB'] == 'duggie']['Fisher p-value'].astype(float)
hc_dgi_stats_df['DUGGIE Mann-Whitney p-value'] = hc_stats_df[hc_stats_df['DTI DB'] == 'duggie']['Mann-Whitney p-value'].astype(float)
hc_dgi_stats_df['STITCH Sig.'] = hc_stats_df[hc_stats_df['DTI DB'] == 'stitch']['Sig.']
hc_dgi_stats_df['STITCH Sig. Treatment'] = hc_stats_df[hc_stats_df['DTI DB'] == 'stitch']['Sig. Treatment']
hc_dgi_stats_df['STITCH Fisher p-value'] = hc_stats_df[hc_stats_df['DTI DB'] == 'stitch']['Fisher p-value'].astype(float)
hc_dgi_stats_df['STITCH Mann-Whitney p-value'] = hc_stats_df[hc_stats_df['DTI DB'] == 'stitch']['Mann-Whitney p-value'].astype(float)

hc_dgi_stats_df['Annotation'] = display_annots
hc_dgi_stats_df.set_index('Annotation', inplace=True)

# +
# Highlight max /min values from each pair of values
hc_dgi_stats_style = hc_dgi_stats_df.style.highlight_max(axis=1, subset=['DUGGIE Sig.', 'STITCH Sig.'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_max(axis=1, subset=['DUGGIE Sig. Treatment', 'STITCH Sig. Treatment'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_min(axis=1, subset=['DUGGIE Fisher p-value', 'STITCH Fisher p-value'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_min(axis=1, subset=['DUGGIE Mann-Whitney p-value', 'STITCH Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
hc_dgi_stats_style.highlight_max(axis=None, subset=['DUGGIE Sig.', 'STITCH Sig.'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_max(axis=None, subset=['DUGGIE Sig. Treatment', 'STITCH Sig. Treatment'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_min(axis=None, subset=['DUGGIE Fisher p-value', 'STITCH Fisher p-value'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_min(axis=None, subset=['DUGGIE Mann-Whitney p-value', 'STITCH Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = hc_dgi_stats_style.format("{:,.2e}", subset=['DUGGIE Fisher p-value', 
                                                     'STITCH Fisher p-value', 
                                                     'DUGGIE Mann-Whitney p-value', 
                                                     'STITCH Mann-Whitney p-value'])
# -

if formatting == "HTML":
    display(HTML(hc_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format="p{3.5cm}|p{1cm}p{1cm}|p{1cm}p{1cm}|p{1.2cm}p{1.2cm}|p{1.2cm}p{1.2cm}",
                caption=r"Table of gene set analysis results analysing hypercholesterolemia from multiple runs of a "
                        r"gene-set analysis pipeline differing in the gene annotation used, comparing the "
                        r"drug-gene interaction data set DUGGIE against STITCH."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:hc_dti_comparison",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Hypercholesterolemia - functional, proximal and hybrid annotation comparison
#
# The hypercholesterolemia results were partitioned by annotation scheme for comparison in tables \ref{tab:hc_annotation_comparison1} and \ref{tab:hc_annotation_comparison2}. As for the identical table \ref{tab:ht_annotation_comparison1} for hypertension, the proximity results per drug-gene interaction dataset are duplicated to aid the comparison but these results are not specific to the listed functional annotation, only to the drug-gene interaction dataset used.
#
# From the significantly identified drugs in table \ref{tab:hc_annotation_comparison1}, the basic functional QTL non-hybrid annotations in contrast to the comparable hypertension results, identify greater numbers of significant drugs in 7 out of 10 annotation categories - for hypertension, it did not achieve this once. This appears to also have the effect of also boosting results from the hybrid functional annotations by increasing the number of significant drugs and significant treatment drugs by a large factor over the numbers achieved by the positional annotations, for example in the best case, eQTL all-tissues with DUGGIE annotation, identifying 129 significant drugs compared to 26 significant drugs identified by the proximity annotation using DUGGIE. This effect is to be expected, as the combined approach takes a greater proportion of its annotations from positionally sourced data than the missing approach. Hence the both hybrid schemes are sensitive to the relative contribution of both positional and functional annotations, with the combined scheme having a greater exposure to the proximity annotation than the missing scheme.
#
# In table \ref{tab:hc_annotation_comparison2}, the Fisher p-values highlight the ability of the QTL functional-based annotations to select treatment drugs when a significance threshold of 0.05 is applied, achieving both the smallest Fisher p-value and the most numerous top results per annotation type. However, this is not repeated for the Mann-Whitney p-values, with all DUGGIE-based functional annotations obtaining higher p-values than the corresponding DUGGIE proximity annotation. The STITCH-based proximity annotation analysis also achieves a lower Mann-Whitney p-value in two out of five comparisons, so functional annotations are not as effective as the proximity annotation at differentiating between treatment and non-treatment drugs when no threshold is applied for hypercholesterolemia.

# +
hc_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'Prox Sig.', 'QTL Sig.', 'Missing Sig.', 'Combined Sig.', 
                                        'Prox Sig. Treatment' , 'QTL Sig. Treatment', 'Missing Sig. Treatment',
                                        'Combined Sig. Treatment',
                                       ])

display_annots3 = [
    "eQTL brain DUGGIE",
    "e/mQTL brain DUGGIE",
    "eQTL blood DUGGIE",
    "eQTL all DUGGIE",
    "e/mQTL all DUGGIE",
    "eQTL brain STITCH",
    "e/mQTL brain STITCH",
    "eQTL blood STITCH",
    "eQTL all STITCH",
    "e/mQTL all STITCH",
]

temp = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('prox')]['Sig.'])
hc_dgi_stats_df['Prox Sig.'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
hc_dgi_stats_df['QTL Sig.'] = list(hc_stats_df[~hc_stats_df['Annotation'].str.contains('hybrid|prox')]['Sig.'])
hc_dgi_stats_df['Missing Sig.'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('hybrid$')]['Sig.'])
hc_dgi_stats_df['Combined Sig.'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('both')]['Sig.'])

temp = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('prox')]['Sig. Treatment'])
hc_dgi_stats_df['Prox Sig. Treatment'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
hc_dgi_stats_df['QTL Sig. Treatment'] = list(hc_stats_df[~hc_stats_df['Annotation'].str.contains('hybrid|prox')]['Sig. Treatment'])
hc_dgi_stats_df['Missing Sig. Treatment'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('hybrid$')]['Sig. Treatment'])
hc_dgi_stats_df['Combined Sig. Treatment'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('both')]['Sig. Treatment'])

hc_dgi_stats_df['Annotation'] = display_annots3
hc_dgi_stats_df.set_index('Annotation', inplace=True)

# +
# Highlight max /min values from each pair of values
hc_dgi_stats_style = hc_dgi_stats_df.style.highlight_max(axis=1, subset=['Prox Sig.', 'QTL Sig.',
                                                                         'Missing Sig.', 'Combined Sig.'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_max(axis=1, subset=['Prox Sig. Treatment',  'QTL Sig. Treatment',
                                                 'Missing Sig. Treatment', 'Combined Sig. Treatment'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
hc_dgi_stats_style.highlight_max(axis=None, subset=['Prox Sig.', 'QTL Sig.',
                                                    'Missing Sig.', 'Combined Sig.'],
                                                        props='color:{Red};')
dummy=hc_dgi_stats_style.highlight_max(axis=None, subset=['Prox Sig. Treatment',  'QTL Sig. Treatment',
                                                 'Missing Sig. Treatment', 'Combined Sig. Treatment'],
                                                        props='color:{Red};')
# -

# \newpage

if formatting == "HTML":
    display(HTML(hc_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format=r"p{3.2cm}|p{0.6cm}p{0.6cm}p{0.8cm}p{1.5cm}|p{0.6cm}p{0.6cm}"
                              r"p{0.8cm}p{1.5cm}",
                caption=r"Table of basic gene set analysis results analysing hypercholesterolemia from multiple runs "
                        r" of a gene-set analysis pipeline differing in the gene annotation used, comparing "
                        r"the functional, proximal and hybrid annotations."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:hc_annotation_comparison1",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
hc_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'Prox Mann-Whitney p-value' , 'QTL Mann-Whitney p-value',
                                        'Missing Mann-Whitney p-value', 'Combined Mann-Whitney p-value',
                                        'Prox Fisher p-value', 'QTL Fisher p-value',
                                        'Missing Fisher p-value', 'Combined Fisher p-value', 
                                       ])

display_annots3 = [
    "eQTL brain DUGGIE",
    "e/mQTL brain DUGGIE",
    "eQTL blood DUGGIE",
    "eQTL all DUGGIE",
    "e/mQTL all DUGGIE",
    "eQTL brain STITCH",
    "e/mQTL brain STITCH",
    "eQTL blood STITCH",
    "eQTL all STITCH",
    "e/mQTL all STITCH",
]

temp = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('prox')]['Fisher p-value'].astype(float))
hc_dgi_stats_df['Prox Fisher p-value'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
hc_dgi_stats_df['QTL Fisher p-value'] = list(hc_stats_df[~hc_stats_df[
                                                'Annotation'].str.contains('hybrid|prox')]['Fisher p-value'].astype(float))
hc_dgi_stats_df['Missing Fisher p-value'] = list(hc_stats_df[hc_stats_df[
                                                'Annotation'].str.contains('hybrid$')]['Fisher p-value'].astype(float))
hc_dgi_stats_df['Combined Fisher p-value'] = list(hc_stats_df[hc_stats_df[
                                                'Annotation'].str.contains('both')]['Fisher p-value'].astype(float))

temp = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('prox')]['Mann-Whitney p-value'].astype(float))
hc_dgi_stats_df['Prox Mann-Whitney p-value'] = [temp[0], temp[0], temp[0], temp[0], temp[0],
                                temp[1], temp[1],  temp[1],  temp[1], temp[1],]
hc_dgi_stats_df['QTL Mann-Whitney p-value'] = list(hc_stats_df[~hc_stats_df[
                                                'Annotation'].str.contains('hybrid|prox')]['Mann-Whitney p-value'].astype(float))
hc_dgi_stats_df['Missing Mann-Whitney p-value'] = list(hc_stats_df[hc_stats_df[
                                                'Annotation'].str.contains('hybrid$')]['Mann-Whitney p-value'].astype(float))
hc_dgi_stats_df['Combined Mann-Whitney p-value'] = list(hc_stats_df[hc_stats_df[
                                                'Annotation'].str.contains('both')]['Mann-Whitney p-value'].astype(float))

hc_dgi_stats_df['Annotation'] = display_annots3
hc_dgi_stats_df.set_index('Annotation', inplace=True)

# +
hc_dgi_stats_style = hc_dgi_stats_df.style

# Highlight max /min values from each pair of values
hc_dgi_stats_style.highlight_min(axis=1, subset=['Prox Fisher p-value', 'QTL Fisher p-value', 
                                                 'Missing Fisher p-value',  'Combined Fisher p-value'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_min(axis=1, subset=['Prox Mann-Whitney p-value', 'QTL Mann-Whitney p-value', 
                                                 'Missing Mann-Whitney p-value', 'Combined Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
hc_dgi_stats_style.highlight_min(axis=None, subset=['Prox Fisher p-value', 'QTL Fisher p-value', 
                                                 'Missing Fisher p-value',  'Combined Fisher p-value'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_min(axis=None, subset=['Prox Mann-Whitney p-value', 'QTL Mann-Whitney p-value', 
                                                 'Missing Mann-Whitney p-value', 'Combined Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = hc_dgi_stats_style.format("{:,.2e}", subset=['Prox Fisher p-value', 
                                                     'QTL Fisher p-value', 
                                                     'Missing Fisher p-value', 
                                                     'Combined Fisher p-value', 
                                                     'Prox Mann-Whitney p-value', 
                                                     'QTL Mann-Whitney p-value', 
                                                     'Missing Mann-Whitney p-value',
                                                     'Combined Mann-Whitney p-value'])
# -

# TODO - just derived values here (as table is too big otherwise)
if formatting == "HTML":
    display(HTML(hc_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format=r"p{3.2cm}|p{1.15cm}p{1.15cm}p{1.15cm}p{1.5cm}|p{1.15cm}p{1.15cm}"
                              r"p{1.15cm}p{1.5cm}",
                caption=r"Table of derived gene set analysis results analysing hypercholesterolemia from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"functional, proximal and hybrid annotations."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:hc_annotation_comparison2",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Hypercholesterolemia - functional annotation tissue sample size comparison
#
# As with hypertension, for hypercholesterolemia the set of analysis results were rearranged to compare the annotations created from different tissue sample sizes, as listed in table \ref{tab:hc_tissue_comparison} which separates the result columns into blood, brain and all tissues, with no blood results for mQTL based annotations as before.
#
# Again blood tissue results, obtained using the highest number of samples, do not in any annotation category or any metric achieve a more favourable result in comparison to brain and all tissue based annotations. For the principal measurement of significant drugs identified, the all tissue results are the best performing in all but one of the annotation categories, achieving the largest result of 129 drugs found to be significantly associated with hypercholesterolemia. 
#
# This trend also translates to the identification by tissue type of treatment drugs, but with the brain tissue based analyses also achieving the highest number of identified drugs many times, identifying just under half of the possible 40 treatment drugs with the eQTL/DUGGIE annotation combination. Brain tissue analyses also broadly made fewer false discoveries of non-treatment drugs than that of all tissue analyses, meaning the Fisher and Mann-Whitney p-values mostly favour brain tissue results. 
#
# This evidence suggests that using QTL data obtained from larger numbers of tissue samples does not increase the number of significant drugs and treatment drugs.

# +
hc_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'Blood Sig.', 'Brain Sig.', 'All Sig.', 
                                        'Blood Sig. Treatment' , 'Brain Sig. Treatment', 'All Sig. Treatment',
                                        'Blood Mann-Whitney p-value' , 'Brain Mann-Whitney p-value',
                                        'All Mann-Whitney p-value',
                                        'Blood Fisher p-value' , 'Brain Fisher p-value', 'All Fisher p-value',
                                       ])

three_na = pd.Series([np.NaN, np.NaN, np.NaN])

display_annots2 = [
    "eQTL DUGGIE",
    "eQTL-miss. DUGGIE",
    "eQTL-comb. DUGGIE",
    "e/mQTL DUGGIE",
    "e/mQTL-miss. DUGGIE",
    "e/mQTL-comb. DUGGIE",
    "eQTL STITCH",
    "eQTL-miss. STITCH",
    "eQTL-comb. STITCH",
    "e/mQTL STITCH",
    "e/mQTL-miss. STITCH",
    "e/mQTL-comb. STITCH",
]

temp = hc_stats_df[hc_stats_df['Annotation'].str.contains('blood')]['Sig.']
hc_dgi_stats_df['Blood Sig.'] = list(temp[0:3].append(three_na).append(temp[3:6]).append(three_na))
hc_dgi_stats_df['Brain Sig.'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('brain')]['Sig.'])
hc_dgi_stats_df['All Sig.'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains('tissues')]['Sig.'])

temp = hc_stats_df[hc_stats_df['Annotation'].str.contains('blood')]['Sig. Treatment']
hc_dgi_stats_df['Blood Sig. Treatment'] = list(temp[0:3].astype(int).append(three_na).append(temp[3:6]).append(three_na))
hc_dgi_stats_df['Brain Sig. Treatment'] = list(hc_stats_df[hc_stats_df[
                                            'Annotation'].str.contains('brain')]['Sig. Treatment'])
hc_dgi_stats_df['All Sig. Treatment'] = list(hc_stats_df[hc_stats_df[
                                            'Annotation'].str.contains('tissues')]['Sig. Treatment'])

temp = hc_stats_df[hc_stats_df['Annotation'].str.contains('blood')]['Fisher p-value'].astype(float)
hc_dgi_stats_df['Blood Fisher p-value'] = list(temp[0:3].append(three_na).append(temp[3:6]).append(three_na))
hc_dgi_stats_df['Brain Fisher p-value'] = list(hc_stats_df[hc_stats_df[
                                            'Annotation'].str.contains('brain')]['Fisher p-value'].astype(float))
hc_dgi_stats_df['All Fisher p-value'] = list(hc_stats_df[hc_stats_df[
                                            'Annotation'].str.contains('tissues')]['Fisher p-value'].astype(float))

temp = hc_stats_df[hc_stats_df['Annotation'].str.contains('blood')]['Mann-Whitney p-value'].astype(float)
hc_dgi_stats_df['Blood Mann-Whitney p-value'] = list(temp[0:3].append(three_na).append(temp[3:6]).append(three_na))
hc_dgi_stats_df['Brain Mann-Whitney p-value'] = list(hc_stats_df[hc_stats_df[
                                            'Annotation'].str.contains('brain')]['Mann-Whitney p-value'].astype(float))
hc_dgi_stats_df['All Mann-Whitney p-value'] = list(hc_stats_df[hc_stats_df[
                                            'Annotation'].str.contains('tissues')]['Mann-Whitney p-value'].astype(float))

hc_dgi_stats_df['Annotation'] = display_annots2
hc_dgi_stats_df.set_index('Annotation', inplace=True)

# +
hc_dgi_stats_style = hc_dgi_stats_df.style

# Highlight max /min values from each pair of values
hc_dgi_stats_style = hc_dgi_stats_df.style.highlight_max(axis=1, subset=['Blood Sig.', 'Brain Sig.', 'All Sig.'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_max(axis=1, subset=['Blood Sig. Treatment', 
                                                 'Brain Sig. Treatment', 'All Sig. Treatment'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_min(axis=1, subset=['Blood Fisher p-value',
                                                 'Brain Fisher p-value', 'All Fisher p-value'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_min(axis=1, subset=['Blood Mann-Whitney p-value', 
                                                 'Brain Mann-Whitney p-value', 'All Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
hc_dgi_stats_style.highlight_max(axis=None, subset=['Blood Sig.', 'Brain Sig.', 'All Sig.'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_max(axis=None, subset=['Blood Sig. Treatment', 
                                                    'Brain Sig. Treatment', 'All Sig. Treatment'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_min(axis=None, subset=['Blood Fisher p-value',
                                                    'Brain Fisher p-value', 'All Fisher p-value'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_min(axis=None, subset=['Blood Mann-Whitney p-value', 
                                                    'Brain Mann-Whitney p-value', 'All Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = hc_dgi_stats_style.format("{:,.2e}", subset=['Blood Fisher p-value', 
                                                     'Brain Fisher p-value', 
                                                     'All Fisher p-value', 
                                                     'Blood Mann-Whitney p-value', 
                                                     'Brain Mann-Whitney p-value', 
                                                     'All Mann-Whitney p-value'], na_rep = 'N/A')
# format ints
dummy = hc_dgi_stats_style.format("{:.0f}", subset=['Blood Sig.', 
                                                    'Blood Sig. Treatment',
                                                   ], na_rep = 'N/A')
# -

# \newpage

if formatting == "HTML":
    display(HTML(hc_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\tiny"))
    display(
        Latex(
            hc_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format=r"l|p{0.35cm}p{0.25cm}p{0.25cm}|p{0.4cm}p{0.4cm}p{0.55cm}|"
                              r"p{1.01cm}p{1.01cm}p{1.01cm}|p{0.95cm}p{0.95cm}p{0.95cm}",
                caption=r"Table of basic gene set analysis results analysing hypercholesterolemia from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"tissue specific functional annotations."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:hc_tissue_comparison",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ### Hypercholesterolemia - contribution of mQTL data to the enrichment of drugs
#
# As a comparison between blood tissue results with and without mQTL data included was not possible due to there being no blood-specific mQTL data sources involved in any annotation, results for the remaining annotations were compiled to make the comparison between annotations with and without mQTL data, found in table \ref{tab:hc_mqtl_comparison}. Both with-mQTL and without-mQTL annotation comparisons show an advantage over the other in many annotations and in each of the four measured metrics, but the number of instances in which the no-mQTL annotations achieve the best result (i.e. green values from table \ref{tab:hc_mqtl_comparison}) is greater for those without-mQTL data, including achieving the top results in all four tabled categories - significant drugs, significant treatment drugs, Mann-Whitney p-values and Fisher p-values. However, including the mQTL data does have some notable significant drug and treatment drug results, mainly for those annotations using brain tissue QTL data.  

# +
hc_dgi_stats_df = pd.DataFrame(columns=[
                                        'Annotation',
                                        'mQTL Sig.', 'no mQTL Sig.', 
                                        'mQTL Sig. Treatment' , 'no mQTL Sig. Treatment',
                                        'mQTL Mann-Whitney p-value' , 'no mQTL Mann-Whitney p-value',
                                        'mQTL Fisher p-value' , 'no mQTL Fisher p-value',
                                       ])

display_annots4 = [
    "Brain DUGGIE",
    "Brain-miss. DUGGIE",
    "Brain-comb. DUGGIE",
    "All DUGGIE",
    "All-miss. DUGGIE",
    "All-comb. DUGGIE",
    "Brain STITCH",
    "Brain-miss. STITCH",
    "Brain-comb. STITCH",
    "All STITCH",
    "All-miss. STITCH",
    "All-comb. STITCH",
]

hc_dgi_stats_df['mQTL Sig.'] = hc_stats_df[~hc_stats_df['Annotation'].str.contains("eqtl|blood|prox")]['Sig.']
hc_dgi_stats_df['mQTL Sig. Treatment'] = hc_stats_df[~hc_stats_df['Annotation'].str.contains(
                                                    "eqtl|blood|prox")]['Sig. Treatment']
hc_dgi_stats_df['mQTL Fisher p-value'] = hc_stats_df[~hc_stats_df['Annotation'].str.contains(
                                                    "eqtl|blood|prox")]['Fisher p-value'].astype(float)
hc_dgi_stats_df['mQTL Mann-Whitney p-value'] = hc_stats_df[~hc_stats_df['Annotation'].str.contains(
                                                    "eqtl|blood|prox")]['Mann-Whitney p-value'].astype(float)

hc_dgi_stats_df['no mQTL Sig.'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~hc_stats_df['Annotation'].str.contains("blood")]['Sig.'])
hc_dgi_stats_df['no mQTL Sig. Treatment'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~hc_stats_df['Annotation'].str.contains("blood")]['Sig. Treatment'])
hc_dgi_stats_df['no mQTL Fisher p-value'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~hc_stats_df['Annotation'].str.contains("blood")]['Fisher p-value'].astype(float))
hc_dgi_stats_df['no mQTL Mann-Whitney p-value'] = list(hc_stats_df[hc_stats_df['Annotation'].str.contains("eqtl") &
                                                   ~hc_stats_df['Annotation'].str.contains("blood")]['Mann-Whitney p-value'].astype(float))

hc_dgi_stats_df['Annotation'] = display_annots4
hc_dgi_stats_df.set_index('Annotation', inplace=True)

# +
# Highlight max /min values from each pair of values
hc_dgi_stats_style = hc_dgi_stats_df.style.highlight_max(axis=1, subset=['mQTL Sig.', 'no mQTL Sig.'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_max(axis=1, subset=['mQTL Sig. Treatment', 'no mQTL Sig. Treatment'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_min(axis=1, subset=['mQTL Fisher p-value', 'no mQTL Fisher p-value'],
                                                        props='color:{Green};')
hc_dgi_stats_style.highlight_min(axis=1, subset=['mQTL Mann-Whitney p-value', 'no mQTL Mann-Whitney p-value'],
                                                        props='color:{Green};')

# Highlight max/min value from all pairs of columns
hc_dgi_stats_style.highlight_max(axis=None, subset=['mQTL Sig.', 'no mQTL Sig.'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_max(axis=None, subset=['mQTL Sig. Treatment', 'no mQTL Sig. Treatment'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_min(axis=None, subset=['mQTL Fisher p-value', 'no mQTL Fisher p-value'],
                                                        props='color:{Red};')
hc_dgi_stats_style.highlight_min(axis=None, subset=['mQTL Mann-Whitney p-value', 'no mQTL Mann-Whitney p-value'],
                                                        props='color:{Red};')

# format floats
dummy = hc_dgi_stats_style.format("{:,.2e}", subset=['mQTL Fisher p-value', 
                                                     'no mQTL Fisher p-value', 
                                                     'mQTL Mann-Whitney p-value', 
                                                     'no mQTL Mann-Whitney p-value'])
# -

if formatting == "HTML":
    display(HTML(hc_dgi_stats_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            hc_dgi_stats_style.to_latex(
                environment="longtable",
                hrules=True,
                column_format="p{3.5cm}|p{1cm}p{1cm}|p{1cm}p{1cm}|p{1.2cm}p{1.2cm}|p{1.2cm}p{1.2cm}",
                caption=r"Table of basic gene set analysis results analysing hypercholesterolemia from multiple runs of a gene-set "
                        r"analysis pipeline differing in the gene annotation used, comparing the "
                        r"annotations including, and not including mQTL data."
                        r"\normalfont{ The best result for each metric is coloured red, and the best result for "
                        r"each metric in each annotation category is coloured green.}",
                label="tab:hc_mqtl_comparison",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# +
drug_sets_df = pd.DataFrame()

for dti in dtis:
    file_postfix = dti + "_qvals.tsv"
    for annot in annots:
        label = dti.upper() + " " + display_annots[annots.index(annot)]
        pval_file = [
            file
            for file in os.listdir(summary_results_dir)
            if file.endswith(file_postfix)
            and "found-" + hc_gwas + "-" + run_id + "-" + annot + "-" in file
        ]
        table = pd.read_csv(os.path.join(summary_results_dir, pval_file[0]), sep="\t")
        table2 = table[table["Q " + annot + " " + dti] < SIGNIF_THRESH]
        drug_set = table2["ATC_CODE"]
        # Add series as column to get a matrix of boolean columns
        drug_sets_df[label] = False
        for drug in drug_set:
            if drug not in drug_sets_df.index:
                drug_sets_df = drug_sets_df.append(pd.Series([False], name=drug))
            drug_sets_df.at[drug, label] = True
        
drug_sets_df = drug_sets_df.fillna(False)
drug_sets_df = drug_sets_df.drop(0, axis=1)
# -

# ### Hypercholesterolemia - overlap of significant drugs found per annotation
#
# UpSet plots were again created for hypercholesterolemia showing the overlap of significantly identified drugs, this time all 32 annotations found at least one significant drug. With significant drug tallies for hypercholesterolemia being higher than those for hypertension, the number of overlaps is also higher so only the head and tail of the resulting UpSet plot are shown, ordered by degree of overlap as before. Figure \ref{fig:hc_upset_0_2} displays the overlaps of 2 degrees or less, whilst figure \ref{fig:hc_upset_9_plus} displays overlaps of 9 degrees or more. The unplotted area between these two shows no surprising data, such as large numbers of overlaps occurring only between the results sets of the least populous drug sets, which would imply that a subgroup of annotations are picking up different drugs unique only to the subgroup and there is no broad consensus over all annotations. 

bold_caption = (
    r"UpSet plot showing the overlap of drug sets identified to be significantly associated with hypercholesterolemia "
    r" from multiple runs of a gene-set analysis pipeline differing in the gene annotation used, for up to 2 "
    r" degrees of overlap."
)
normal_caption = (
    r"Vertical bars show the gene set details per named annotation, with a histogram showing the significant"
    r" gene tallies per annotation. The horizontal histogram shows the size of overlap between the sets indicated"
    " horizontally in the dot matrix. Single dots represent the tally of genes unique to the named annotation gene"
    r" set."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

usp_figure = plt.figure(dpi=200)
dummy = usp.plot(usp.from_indicators(drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=0,
                 max_degree=2)

rdisp.end_caption(formatting, r"fig:hc_upset_0_2")

# \newpage

# Half, 16, of these annotation sets have unique drugs not found by any other annotation, and the size of these unique drug sets are approximately proportional to the total number of significant drugs identified by that annotation, with the most successful annotation, all-tissue eQTLs with combined proximity scheme and DUGGIE drug-gene interaction dataset, finding 129 drugs of which 34 are unique to the annotation. The least successful 6 annotations identify no unique significant drugs.

bold_caption = (
    r"UpSet plot showing the overlap of drug sets identified to be significantly associated with "
    r"hypercholesterolemia from multiple runs of a gene-set analysis pipeline differing in the gene annotation "
    r"used, for 9 or more degrees of overlap."
)
normal_caption = (
    r"Vertical bars show the gene set details per named annotation, with a histogram showing the significant"
    r" gene tallies per annotation. The horizontal histogram shows the size of overlap between the sets indicated"
    " horizontally in the dot matrix. Single dots represent the tally of genes unique to the named annotation gene"
    r" set."
)
rdisp.start_caption(formatting, bold_caption, normal_caption)

usp_figure = plt.figure(dpi=200)
dummy = usp.plot(usp.from_indicators(drug_sets_df), 
                 fig=usp_figure,
                 subset_size='count', 
                 show_counts=True, 
                 orientation='vertical',
                 min_degree=9,
                 max_degree=100)

rdisp.end_caption(formatting, r"fig:hc_upset_9_plus")

# Figure \ref{fig:hc_upset_9_plus} of the highest degree overlaps for hypercholesterolemia shows a similar trend to hypertension, where the bottom left of the matrix has a high density of dots indicate that the drugs occurring most commonly appear in the sets with the highest tallies. This indicates that the most successful annotations are building upon the significant drug discoveries of the less successful annotations and that the analysis method itself is fairly robust.

# ### Hypercholesterolemia - drug enrichment across annotation results
#
# Table \ref{tab:hc_most_frequent_drugs} lists the drugs that are most frequently enriched across results from the pipeline executions using the set of 32 different annotations. Those drugs appearing in 15 or more analysis results are listed, giving 20 drugs in total with the most frequent, Probucol, appearing in 30 of the 32 analyses. From this list of 20, 7 drugs are classed as treatment drugs by DrugBank whilst evidence that the drug affects cholesterol was not found for only two significant drugs - procainamide and azelastine.
#
# At the other end of this enrichment frequency scale, there are 82 drugs unique to only one of the set of 32 annotation results, three of which are classed as treatment - Atorvastatin,  Lovastatin with nicotinic acid and Atorvastatin with perindopril.

top_drugs_df = pd.DataFrame(drug_sets_df.transpose().sum(), columns=['Analysis count'])
top_drugs_df = top_drugs_df[top_drugs_df['Analysis count'] > 15]
top_drugs_df = pd.merge(top_drugs_df, ATC_LEVELS, left_index=True, right_index=True)
top_drugs_df.sort_values(by=['Analysis count'], ascending=False, inplace=True)
top_drugs_df = rdisp.add_treatment_state(hc_gwas, top_drugs_df)
top_drugs_df.index.set_names(["ATC code"], inplace=True)

if formatting == "HTML":
    display(HTML(top_drugs_df.to_html()))
if formatting == "LaTeX":
    display(Latex(r"\scriptsize"))
    display(
        Latex(
            top_drugs_df.to_latex(
                column_format=r"p{3.2cm}p{1.5cm}p{5cm}p{1.5cm}p{1.5cm}",
                caption=r"Table of drugs most frequently enriched across annotations. ATC drug codes are listed "
                        r"alongside the count of analyses in which the drug is significantly enriched, if the drug "
                        r"is used to treat hypercholesterolemia and if the drug is known to affect cholesterol levels.",
                label="tab:hc_most_frequent_drugs",
                position="htbp",
            )
        )
    )
    display(Latex(r"\normalsize"))

# ## Discussion
#
# This chapter set out to identify annotation methods providing the most appropriate set of drugs for drug repurposing when used as input into a gene-set analysis pipeline. It did this by testing several hypotheses, which broadly state that significant drugs and treatment drugs are enriched in gene set analysis pipeline results by annotating genes functionally with QTL data obtained from the largest samples, and augmenting these functional annotations with positional data.
#
# Each of these hypotheses were considered independently by partitioning the same multi-dimensional matrix of results according to features relevant to that hypothesis. The matrix was also used to support previous conclusions that treatment drugs are generally enriched in the analysis pipeline results and that the larger DUGGIE drug-gene interaction dataset has the greatest potential to uncover novel drug repurposing opportunities.
#
# Fisher p-values were used to indicate if treatment drugs are enriched among the significant drugs. Specifically, the Fisher p-value measures the ability of the classification of significant/non-significant drugs at a 0.05 Storey FDR threshold to differentiate between treatment and non-treatment drugs. Only a minority of the annotations, around 15%, did not achieve a significant Fisher p-value, 7 of 32 results for hypertension (table \ref{tab:ht_dti_comparison}) and 3 of 32 results for hypercholesterolemia (table \ref{tab:hc_dti_comparison}). This is good evidence that the gene set analysis pipeline can identify known treatment drugs, and therefore increases confidence that any significant non-treatment drugs discovered may be novel therapeutics and are good candidates for further investigation.
#
# Once more DUGGIE was seen to identify a greater number of significantly associated drugs than STITCH for both hypertension and hypercholesterolemia and obtaining the largest number of such drugs for any annotation in both diseases. As in chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}, DUGGIE and STITCH found similar numbers of treatment drugs. DUGGIE identified a greater number of significant non-treatment drugs than STITCH, which translated into poorer Fisher p-values overall. DUGGIE fared better with regard to discerning between treatment and non-treatment drugs with no threshold applied as given by the Mann-Whitney p-values, achieving lower p-values in the majority of annotation comparisons and obtaining the lowest Mann-Whitney p-value overall for both hypertension and hypercholesterolemia. Despite DUGGIE having an advantage over STITCH in significantly associating more drugs, the ability of STITCH alone to identify treatment drugs means that it may be useful to run analyses using both drug-gene interaction datasets, as concluded in chapter \ref{appraisal-of-the-pipeline-performance-with-hypertension-and-hypercholesterolemia}.
#
# On considering the enrichment of drugs and treatment drugs between the functional and proximal annotations for each disease as shown in tables \ref{tab:ht_annotation_comparison1} and \ref{tab:hc_annotation_comparison1} it can be seen that hypertension and hypercholesterolemia results have some notable differences. For hypertension, the proximity annotation finds a greater number of significant drugs than the functional annotations in most cases, whilst for hypercholesterolemia the functional annotations are more successful. This difference is repeated for the treatment drug statistics as seen in tables \ref{tab:ht_annotation_comparison2} and \ref{tab:hc_annotation_comparison2} in the Fisher p-values and to a lesser extent the Mann-Whitney p-values. With some specific functional annotations outperforming their proximally annotated counterparts it can be said that there is evidence supporting the hypothesis that functionally annotating genes increases the enrichment of treatment drugs, but with the caveat that careful consideration of the functional annotation is important. However, no further specific recommendations can be made as there is no obvious pattern from the diseases and annotations studied so far.
#
# The same can also be said towards the increase in enrichment of treatment drugs by functional annotations augmented with positional data. At least one of the missing or combined hybrid annotations consistently shows a marked improvement over their respective QTL-only functional annotation in almost all cases. So again, there is evidence supporting the hypothesis that augmenting functional annotations with positional annotations increases the enrichment of treatment drugs but careful consideration of the functional annotation is important to achieve this.
#
# From investigating the effect of larger functional annotation tissue sample sizes on hypertension (table \ref{tab:ht_tissue_comparison}) and hypercholesterolemia (table \ref{tab:hc_tissue_comparison}), it is plain that the results do not show an increase in significant drugs or significant treatment drugs for the blood tissue based annotations which have the highest sample sizes. Even when the number of annotated genes involved in the analysis is corrected to be of similar number, as is the case for the missing and combined hybrid approaches (see table \ref{tab:db_annotations}), identification of significant drug and treatment drug numbers by the higher sample size blood tissue annotations is still lacking compared to the brain and all tissue annotations, which are obtained from lower tissue sample sizes. There is no evidence to support the hypothesis that using QTL data obtained from larger numbers of tissue samples increases the number of significant drugs and treatment drugs in analysis pipeline results, and potentially this is evidence against the hypothesis. However, one confounding effect could be that the tissue types chosen on which to base the comparison - blood, brain and all tissues - have QTLs that as a set, may differ substantially with respect to their association to disease genes and this is likely to bias the results in a currently unknown way. Further work, out of scope for this project, could be carried out to determine which tissue type QTLs show such a bias and correct it if possible.
#
# Finally the effect of the extra brain-tissue specific mQTL data can be considered. These results, shown in tables \ref{tab:ht_mqtl_comparison} and \ref{tab:hc_mqtl_comparison}, are mixed but there are many which are favourable towards annotations involving mQTL data. These are inconclusive in supporting the hypothesis that including specific brain tissue mQTL data alongside eQTL data increases the enrichment of treatment drugs.
#
# In moving on to consider neuropsychiatric diseases, there is no clear choice between many of the annotation variations explored here. It can be seen that the analysis pipeline is useful in selecting drugs known to treat the studied disease, that DUGGIE gives a richer possibility of novel treatment drugs and that functional annotations augmented by positional data can outperform a more basic proximity annotation. However, the choice of positional augmentation scheme or tissue specific functional data is not clear.
#
# As the UpSet plots of drug set overlaps (figures \ref{fig:ht_upset_4_plus} and \ref{fig:hc_upset_9_plus}) indicate that the more successful annotations broadly capture the same drugs as the less successful annotations, it would be advisable to perform executions of the pipeline on neuropsychiatric diseases covering both DUGGIE and STITCH, brain and all tissues (with and without mQTL data) using both the missing and combined annotation schemes. Without requiring a detailed aetiological analysis and relying on the assumption that the most powerful annotations are the most successful in identifying the greatest number of significant drugs, this would allow the results from the best performing pipeline execution results to be selected and analysed further. 
#
# Furthermore, there are some drawbacks to this approach as any significant disease-associated drug identified may be backed by an array of association q-values, one from each pipeline execution, that are obtained from more than one set of drug-gene annotations i.e. DUGGIE and STITCH. Hence identifying the genes driving the result is more complex, but still possible.
#
# However, drugs appearing in more analysis pipeline results using different annotations were seen to more likely be treatment drugs and drugs only found by one analysis pipeline were more likely to be non-treatment drugs. Hence the selection of drugs that are potential repurposing candidates could be aided by ranking significant drugs by the count of pipeline execution results in which they appear.
