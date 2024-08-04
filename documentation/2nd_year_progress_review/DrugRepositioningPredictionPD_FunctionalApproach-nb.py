# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
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

# #!pip install seaborn
import seaborn as sns
import pylab
import scipy.stats as stats
from scipy.stats import chi2, norm

# from latex_envs.latex_envs import figcaption

# Also generate tables in latex format, with formatting opts
pd.set_option("display.latex.repr", True)
pd.set_option("display.latex.longtable", True)
pd.set_option("display.latex.multirow", True)
pd.set_option("display.latex.multicolumn", True)
pd.set_option("display.max_colwidth", 135)

# For plots
fig_unit = 16
fig_width = 3
fig_height = 2
# -

# \newpage
# # Abstract
#
# Gene set analysis is an illuminating method of identifying groups of genes most associated with a disease by testing for association between a gene set mapped from Genome Wide Association Study (GWAS) results and each gene set of interest.
#
# Based on the hypothesis that drugs targeting proteins associated with genes that increase risk to disease will reveal novel therapeutic opportunities, and that these can be used to inform and improve drug-based treatment, the association between groups of drug protein targets and genome regions associated with Parkinson's disease was investigated using gene set analysis.
#
# Utilising the results from a recent Parkinson's GWAS meta study, a gene-wise association analysis using the MAGMA tool<citet data-cite="leeuw_magma:_2015"><sup>magma</sup></citet> was performed with functional Quantitative Trait Locus (QTL) Single Nucleotide Polymorphism (SNP) to gene mappings to obtain a set of gene associations for Parkinson's disease. These association results were used to carry out a gene-set analysis against sets of drug target genes grouped using the Anatomical Therapeutic Chemical (ATC) classification system<citet data-cite="noauthor_anatomical_nodate"><sup>atc</sup></citet> which hierarchically classifies approved drugs according to their pharmacological, therapeutic and chemical characteristics.
#
# Several significant drug group-gene associations were observed, mainly in the immunoglobulins and immunomodulating agent ATC categories suggesting that the immune system has a role in the aetiology of Parkinson's disease, and possible approved-drug repurposing possibilities in this drug group were identified.
#
# The automated code used for this study is provided<citet data-cite="einon_genomic_nodate"><sup>gnt-gitlab</sup></citet> and can be easily modified to analyse other disease GWAS results against the provided, or other functional & classification gene sets of interest. All code is GPLv2 licensed.

# \newpage
# # Introduction
#
# Parkinson's is the second most common neurodegenerative disease<citet data-cite="corti_what_2011"><sup>corti</sup></citet> and is expected to increase in an aging population, escalating an already large social and economic burden on societies. This debilitating condition, understood to be caused by a slow loss of nerve cells in a part of the midbrain known as the substantia nigra, leading to a drop in the production of dopamine and associated degradation of motor control<citet data-cite="dauer_parkinsons_2003"><sup>Dauer</sup></citet>. A lack of detailed knowledge regarding the mechanism of nerve cell loss is currently hampering the progression of effective treatments.
#
# Previous studies have demonstrated that for genes associated with a disease, if the proteins for which they encode are associated with targets for existing medications<citet data-cite="ruderfer_polygenic_2016"><sup>ruderfer</sup></citet> then it follows that these genes are clear targets for drug studies<citet data-cite="stein_effect_2012,navarese_effects_2015"><sup>stein,navarese</sup></citet> and offer opportunities to re-purpose existing approved drugs. It has also recently been demonstrated that drugs supported by genetic associations are more likely to successfully advance through the later stages of clinical trials<citet data-cite="nelson_support_2015"><sup>nelson</sup></citet>.
#
# A gene set analysis requires at least two sets of genes - one obtained from disease GWAS SNP data with p-values mapped to genes, and another from a gene set of interest - in this case a set of drug target genes defined by one of the ATC classes. Each ATC gene set undergoes a competitive gene-set analysis, comparing the mean association of genes in the gene set with Parkinson's disease genes, against the association of genes not in the gene set with Parkinson's disease genes. Several covariates are taken into account which are described further below. A low gene-set p-value implies an increased association of the gene-set with Parkinson's disease against all genes in the analysis and the most significant associations could signify a potential drug repurposing opportunity.
#
# Hence to prepare for a gene set analysis first the SNP-based p-values obtained from a GWAS must be translated to gene based p-values - a two step process consisting of a gene annotation step and a gene analysis step.
#
# ## Gene annotation
#
# Gene annotation in this case is the assignment of SNPs to genes, based on the SNP's effect or function. A large number of SNPs involved in GWAS are non-coding but functionally relevant, and assignment of these compared to assignment of coding SNPs to genes is less straightforward.
#
# The most accessible strategy for annotation is a positional or proximity based method where a SNP is assigned to the nearest gene on the same chromosome. However, this approach doesn't take into account any further known biological evidence and has been often found to assign SNPs to genes that differ from other sources of evidence such as QTL expression measurements that can, in many cases, implicate the greatest effect on a gene close, but not adjacent to the SNP being considered.
#
# To attempt to improve the biological basis for SNP annotation, QTL data can be used. A QTL is a region of DNA whose variation correlates with the variation of a quantitative trait such as an mRNA expression level (eQTL) or CpG methylation of a DNA region (mQTL). These DNA loci are identified by a direct association test between SNPs and a quantitative measurement of the phenotype (e.g. methylation or gene expression) from tens or hundreds of individuals<citet data-cite="nica_expression_2013"><sup>nica</sup></citet>, and provide more direct biological evidence than associating SNPs positionally. As these measured traits can vary by tissue type, QTL data sets are often stratified thus and this provides an opportunity to select a subset of tissue types if enough aetiology is known of the disease being investigated to allow - this reduces the multiple testing required, preserving statistical power. For Parkinson's, brain specific QTLs are an obvious choice and are used here, but more recent research<citet data-cite="killinger_vermiform_2018"><sup>killinger</sup></citet> also indicates that the gut and appendix tissues may be of interest to include.
#
# As QTL data may not include results significant enough to identify SNPs for all known genes, a *hybrid* approach where genes not annotated by any QTL derived SNPs are annotated using positionally derived SNP lists - this increases the coverage of available genes, giving more opportunity to find novel associations. A further refinement to this strategy (not used yet in this study), could be to attempt to provide SNP annotations from QTL data for tissue types not chosen for the disease-pertinent set (i.e. brain tissues for Parkinson's), before attempting to fill in the missing data from positional annotations in a 3-tier hybrid approach.
#
# ## Gene analysis
#
# The annotation output consists of a list of SNPs for each gene that can be used to translate a set of disease-SNP association p-values, as generated from a GWAS, into a set of disease-gene association p-values. There are many methods and tools to perform this step - the SNP-mean and SNP-top multiple model (or *multi=snp-wise*) as implemented by the MAGMA tool was used for this study. For each gene, this model takes the list of SNPs annotated to it, along with the SNP-disease association p-values and calculates both the mean p-value (which MAGMA documentation refers to as the *snp-wise=mean* model), and the top 1 p-value (referred to as *snp-wise=top*). The two resultant p-values, once from each model, are combined to give a joint p-value that represents the gene's association with the disease. MAGMA converts the joint p-value to a Z-score using a probit transformation. The set of such Z-scores for all annotated genes can be used as input to a gene set analysis.
#
# The mean model on it's own measures the mean SNP association but tends to skew towards associations in areas of higher LD within a gene, while the top model on it's own, only considers the most significant SNP is very sensitive to genes with few SNP associations. By aggregating these scores into a single p-value, both these unwanted effects are blunted. Alongside the gene p-values, a MAGMA gene analysis also calculates correlations between all pairs of genes which are used to account for Linkage Disequilibrium (LD) in any subsequent gene set analysis. Linkage disequilibrium is the non-random association of allele loci along a genome. Loci are in LD when the frequency of association of the alleles is not equal to that expected if they were random - this can occur due to genetic linkage, mutation, genetic drift and
# population structure.
#
# ## ATC drug / target gene sets
#
# The ATC classification system (along with it's associated Defined Daily Dose, DDD) is designed to exchange and compare data on drug utilisation worldwide, and managed by the World Health Organization (WHO). It catalogues a non-exhaustive set of approved, licensed drugs by a 7 character code which is split into 5 levels:
#
# 1) Anatomical/pharmacological main group
#
# 2) Pharmacological/therapeutic main group
#
# 3) and 4) Chemical/therapeutic/pharmacological subgroup
#
# 5) Chemical substance
#
# The 2nd, 3rd and 4th levels are often used to identify pharmacological subgroups where appropriate, instead of therapeutic or chemical subgroups. For example, the level 5 ATC code (chemical substance) for *metaformin* is **A10BA02**. This substance also belongs to 4 other groups, one at each level, which can also be identified from it's level 5 code: **A** is it's 1st level, anatomical main group code - *Alimentary tract and metabolism*. **A10** is the 2nd level code, therapeutic subgroup - *Drugs used in diabetes*. **A10B** is the 3rd level, pharmacological subgroup - *Blood glucose lowering drugs, excl. insulin* and finally **A10BA** is the 4th level code, chemical subgroup - *Biguanides*.
#
# In effect, if all level 5 coded substances are taken as a set, this is also the complete set of substances described by all level 4 codes, or level 3 codes etc, just grouped differently. It follows that the set of protein/gene targets of drugs at level 5 is also the same set of genes, grouped differently, at each of the other 4 levels.
#
# Using available drug-target databases it is possible to populate a table of level 5 chemical substances with known protein / gene targets of each drug.
#
# As drugs can have several different therapeutic or pharmacological applications (e.g. ethanol can be used topically or orally), the same substance can appear in different groups at the same level. For subsequent analysis using gene sets based on the ATC, each level is de-duplicated so identical drug sets (i.e. ATC codes having matching gene sets) only appear once to reduce the effects of multiple testing.
#
# ## Gene set analysis
#
# The MAGMA tool can once again be used to perform the final step in the pipeline, a competitive gene set analysis implemented as a linear regression of a gene data matrix, where each gene is marked either as a member of the test set or outside this set and number of SNPs, Z-score, gene size and gene-gene correlations (calculated as part of the previous gene analysis) included as standard covariates. The calculated gene-set p-value is the result of a test that the mean association of genes in the gene set is greater that that of genes not in the gene set, conditioned broadly on the standard covariates listed.
#
# A separate gene set analysis is carried out for each ATC level, as each level consists of the same gene set but grouped according to different criteria.
#
# The effect of the Major Histocompatibility Complex (MHC) on the gene set analysis may also be investigated. The MHC is a region on chromosome 6 having a high density of important immune related genes, and any statistical association attributed to a region in close linkage disequilibrium with the MHC may be falsely attributed to MHC genes which are know to have a high association with brain-related diseases such as Parkinson's and Schizophrenia. To gauge the effect of the MHC on the gene set analysis, a second run of the analysis pipeline against the ATC genesets excluding genes within the MHC region can be made.
#
# ## Gene Ontology pathways
#
# The Gene Ontology<citet data-cite="ashburner_gene_2000,thegeneontologyconsortium_expansion_2017"><sup>go1,go2</sup></citet> (GO) is a comprehensive resource for information regarding the functions of genes and gene products. It describes functional relationships between genes, known as *GO terms*, that describe molecular function (molecular-level activities performed by gene products), cellular component (the locations relative to cellular structures in which a gene product performs a function), and biological process (processes accomplished by multiple molecular activities). As each GO term is associated with a set of interacting genes, also known as a *pathway*, it follows that GO terms can therefore be used in a gene set analysis with each term a set to be analysed against all genes not in the set, being the collection of all other existing GO term genes.

# \newpage
# # Materials and Methods
# Summary statistics from a recent meta-analysis of Parkinson's disease GWAS<citet data-cite="chang_meta-analysis_2017"><sup>chang</sup></citet> were used as input to perform the prerequisite gene analysis along with an annotated hybrid set of genes generated from functional QTLs, then populating any missing genes with positional annotations if available.
#
# All genome datasets used are based on the GRCh37 genome build.
#
# The MAGMA gene analysis tool v1.06 was used to perform both gene and gene-set analyses which were orchestrated using bash shell scripts. Data wrangling was mainly performed using the Python 3 pandas, numpy and scipy libraries implemented in Jupyter Notebooks. Linux HPC clustered hardware was used to perform the analysis, either a departmental cluster ('Rocks') running the SGE scheduler or a university cluster ('Hawk'), managed by ARCCA, running the slurm scheduler - scripts suitable for running on both of these, or similar systems, are also provided.
#
# ## Annotate genes with functional variant SNPs

# <img src="./create_annotation.png">

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\includegraphics{create_annotation}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:annotate} - Overview of software pipeline steps to annotate genes. \normalfont{Grey boxes and green diamonds represent processing steps implemented by the named jupyter notebooks. Blue boxes represent input/output data sets.}}"
    )
)
display(Latex(r"\label{fig:annotate}"))
display(Latex(r"\end{figure}"))

# Several Quantitative Trait Loci (QTL) datasets related to the brain were taken - all 13 GTEx<citet data-cite="carithers_novel_2015"><sup>gtex</sup></citet> brain expression QTLs (eQTLs) from the following regions: putamen basal ganglia, cervical c1 spinal cord, substantia nigra, amygdala, anterior cingulate cortex, caudate basal ganglia, cerebellar hemisphere, cerebellum, cortex, frontal cortex, hippocampus, hypothalamus and the nucleus accumbens basal ganglia; the Common Mind Consortium (CMC, https://www.synapse.org/#!Synapse:syn4622659) set of dorsal lateral pre-frontal cortex brain eQTLs and a set of (CNGG) pre-frontal cortex methylation QTLs (mQTLs) were also used. QTLs with absolute effect sizes less that 0.1 were discarded as having a negligible effect as well as SNPs on the X chromosome to reduce unnecessary statistical multiple testing. SNP chromosome positions were converted to rs numbers using a lookup table provided by GTEx. Each QTL dataset was transformed to a table indexed by gene with a set of associated significant SNP variants listed for that gene. A final list of functional gene annotations was created by combining all such QTL dataset gene tables and removing duplicate SNPs from the list of associations for each gene.
#
# Additionally, where a gene had no gene-variant SNP associations functionally mapped but had a positional mapping present in the MAGMA provided default annotated gene list; the positional mapping was added to the gene annotation table to create a hybrid set of gene to variant mappings.
#
# ## Curate a gene / drug ontology database
#
# The full DrugBank database v5.1.0<citet data-cite="noauthor_drugbank_nodate"><sup>drugbank</sup></citet> was downloaded locally, which explicitly annotates drugs with lists of protein targets and ATC code (if one exists). Drugs that had at least one protein target with an ATC code were filtered and the drug name, protein targets and ATC code extracted. Protein identifiers were converted to genes using UniProt. ATC codes with less than 5 protein targets were discarded - this could be increased to the more accepted minimum of 10 genes per set if more sources of ATC drug gene targets are found. ATC codes at the same ATC level having duplicate target gene lists were removed to leave only one in order to reduce multiple testing. This analysis was repeated using the set of ATC codes at each of the ATC levels 1-5 to obtain 5 separate datasets.
#
# The effect of Major Histocompatibility Complex (MHC) genes on the analysis was also considered, looking at the role of any genes that cover the region on chromosome 6 between bases 26000000 and 33000000.
#
# In an effort to include as many possible ATC coded drugs in the analysis, a more exploratory approach was employed by performing a fuzzy text comparison between the drug names of drugbank entries with protein targets but no ATC code, and ATC level 5 code drug names using the python fuzzywuzzy package. This package was used to perform a pairwise comparison between all possible pairs from the two drug name sets and scores the similarity on a scale between 0 and 100, with 100 being given to a perfect text match and 0 to text with no similarity. Any confirmed matches or alternate spellings found, for example, were added to the gene / drug database.
#
# ## Associate genes with Parkinson's disease

# <img src="./magma_analysis.png">

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\includegraphics{magma_analysis}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:magma} - Overview of MAGMA analysis pipeline steps. \normalfont{Grey boxes and green diamonds represent processing steps implemented by the named jupyter notebooks. Blue boxes represent input/output data sets.}}"
    )
)
display(Latex(r"\label{fig:magma}"))
display(Latex(r"\end{figure}"))

# The gene analysis step takes GWAS summary statistics, in this case a set of Parkinson's disease / SNP association p-values and transforms them into a set of Parkinson's disease / gene Z-scores suitable for use in a gene set analysis.
#
# The Parkinson's disease GWAS meta-data summary results were converted to use SNP rs numbers before being passed to MAGMA and as summary statistics are used as opposed to raw genotype data, a reference for Linkage Disequilibrium (LD) calculations between genes is required which was provided by the 1000 genomes project phase 3 european panel data set (downloaded from the MAGMA website).
#
# Some exploratory investigations were carried out to determine the most suitable MAGMA gene analysis settings (see the relevant results section below). For the reported gene set results, a MAGMA gene-model setting of *multi=snp-wise* and a *gene-settings fixed-permp* option of 20,000 permutations were used. The gene set analysis was run using the MAGMA batch mode in batches of between 10 and 30, depending on the hardware specification of the Linux cluster running the analysis.
#
# The *multi=snp-wise* setting assigned a gene p-value using a model that combines the mean and top (1) p-values of all SNPs annotated to that gene.
#
# The *fixed-permp* gene settings is used to modify the permutation-based model to calculate a gene's top p-value. By default the permutation model is adaptive, where permutations are run in small batches and stopped when a threshold of permutations greater than the observed p-value is counted so the number of permutations can very between successive runs with the same inputs. The fixed permutations setting instead runs a fixed number of permutations to calculate the empirical p-value.
#
# At this stage the suitability of the choice of QTLs used to create gene annotations along with the hybrid annotation approach used was investigated by comparing p-value histograms and correlations between the three annotation approaches - functional, positional and hybrid functional-positional.
#
# ##  Identify associations between pharmacologically related sets of drug targets and Parkinson's disease
#
# The drug ontology database and disease association gene data were used as inputs in a gene enrichment analysis performed again using MAGMA, testing whether the genes included in each pharmacologically related gene-set are more strongly associated with the disease than others in a competitive gene-set analysis. The structure of the ATC hierarchy allows sets of genes characterised by anatomical, chemical, pharmacological and therapeutic relationships to be tested for enrichment. As the gene-sets are redefined from the same pool of gene for each of these properties, one per ATC level, an independent gene-set analysis was conducted for each of the 5 ATC levels to minimise multiple testing. An additional gene set analysis was carried out using a set of annotated Gene Ontology (GO) terms.
#
# Conditional effects such as gene size and correlation between neighbouring genes caused by LD were accounted for by MAGMA using values calculated during the previous gene analysis step.
#
# The output gene set association p-values for each gene set analysis were corrected for multiple testing by calculating q-values (Storey)<citet data-cite="storey_direct_2002"><sup>storey</sup></citet>, which adjust the p-values based on an optimised False Discovery Rate (FDR) approach. A significance threshold of less than 5% was chosen for the resulting q-value to be listed in the summary results.
#
# We then summarised the results over all 5 ATC levels in a table indexed by the ATC code of the gene set and ranked by q-value. The top genes (truncated to 50 at most for brevity) of each significant gene set were listed by ensembl and HGNC code <citet data-cite="noauthor_custom_nodate"><sup>HGNC</sup></citet> in decreasing absolute Z-score order.
#
# Additionally, the effect of using the hybrid functional-positional approach on the final gene set results was measured at this stage by running a more usual pathway gene set analysis alongside the ATC gene sets, and also running an analysis of all 5 ATC gene sets and pathways against the functionally and positionally annotated GWAS gene p-values.
#
# All the results described above were also listed ordered by the uncorrected p-values, filtered by a 5% cut-off in order to track gene set results between approaches that weren't significant in all.

# \newpage
# # Results

# +
data_path = "../../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"
pd.set_option("display.max_colwidth", -1)

gwas = "PD"
qtls = "brain"
run_id = "-no_kaviar"

## Summary table functions
# Returns ATC data (ATC code & description) in a dataframe, given an ATC table file of the same
def get_atc_file(filename):
    result = pd.DataFrame()
    with open(filename, "r") as file:
        lines = file.readlines()
        for line in lines:
            code = line.split(" ")[0]
            desc = " ".join(line.split(" ")[1:]).strip()
            result = result.append(
                {"atc_code": code, "Description": desc}, ignore_index=True
            )

    return result


# Returns ATC data from all 5 ATC levels in one dataframe
def get_atc_levels():
    atc_level_files = [
        data_path + "atc/level-" + str(i) + "-atc.txt" for i in range(1, 6)
    ]

    # convert each atc_level_file to a 2-column data frame, and cat all of them together
    atc_levels = pd.concat(
        [get_atc_file(atc_level_file) for atc_level_file in atc_level_files]
    )
    atc_levels.set_index("atc_code", inplace=True)
    return atc_levels


# Read all results from a magma geneset analysis for this gwas and annotation type
def read_results_files(gwas, annot_type):
    prefix = (
        magma_results_dir
        + gwas
        + "/results/"
        + "magma_geneset_result-"
        + gwas
        + "-"
        + annot_type
        + run_id
    )
    results_fileset = [prefix + "-atc" + str(i) + ".sets.out" for i in range(1, 6)]

    return [
        pd.read_csv(file, comment="#", delim_whitespace=True)
        for file in results_fileset
        if os.path.exists(file)
    ]


# Summarise results by applying significance and number of gene thresholds, and store in a file
# Generates results for both P and Q values.
def summarise_drug_results(
    gwas, annot_type, atc_levels, n_gene_thresh=5, signif_thresh=0.05
):
    results = read_results_files(gwas, annot_type)

    # only consider classes/drugs with Qval < QVAL_THRESH and NGENES >= N_GENE_THRESH
    significants = [result[result["NGENES"] >= n_gene_thresh] for result in results]

    if not significants:
        print("returning - no significants for " + gwas + " " + annot_type)
        return

    # Calculate the q-values for the local results per GWAS (per ATC code)
    for result in significants:
        result.loc[:, "Q"] = qvalue.estimate(np.array(result["P"]))

    # Put all individual ATC results into one big dataframe
    all_significant = pd.concat(significants)

    q_significant = all_significant[all_significant["Q"] < signif_thresh]
    q_final = pd.merge(
        q_significant, atc_levels, right_index=True, left_on="SET"
    ).sort_values("Q")
    q_final.to_csv(
        os.path.join(
            summary_results_dir,
            "drugs_found-" + gwas + run_id + "-" + annot_type + "_qvals.tsv",
        ),
        sep="\t",
        index=False,
    )

    p_significant = all_significant[all_significant["P"] < signif_thresh]
    p_final = pd.merge(
        p_significant, atc_levels, right_index=True, left_on="SET"
    ).sort_values("P")
    p_final.to_csv(
        os.path.join(
            summary_results_dir,
            "drugs_found-" + gwas + run_id + "-" + annot_type + "_pvals.tsv",
        ),
        sep="\t",
        index=False,
    )


def summarise_gopath_results(gwas, annot_type, n_gene_thresh=5, signif_thresh=0.05):
    prefix = magma_results_dir + gwas + "/results/"
    file = (
        prefix
        + "magma_geneset_result-"
        + gwas
        + "-"
        + annot_type
        + run_id
        + "-gopaths.sets.out"
    )

    if not os.path.exists(file):
        display("File not found: " + file)
        return

    result = pd.read_csv(file, comment="#", delim_whitespace=True)

    # only consider classes/drugs with Qval < QVAL_THRESH and NGENES >= N_GENE_THRESH
    significants = result[result["NGENES"] >= n_gene_thresh]

    if significants.empty:
        return

    significants.rename(columns={"SET": "SHORT_NAME", "FULL_NAME": "SET"}, inplace=True)

    # Calculate the q-values for the local results per GWAS (per ATC code)
    significants.loc[:, "Q"] = qvalue.estimate(np.array(result["P"]))

    q_significant = significants[significants["Q"] < signif_thresh]
    q_significant.to_csv(
        os.path.join(
            summary_results_dir,
            "pathways_found-" + gwas + run_id + "-" + annot_type + "_qvals.tsv",
        ),
        sep="\t",
        index=False,
    )

    p_significant = significants[significants["P"] < signif_thresh]
    p_significant.to_csv(
        os.path.join(
            summary_results_dir,
            "pathways_found-" + gwas + run_id + "-" + annot_type + "_pvals.tsv",
        ),
        sep="\t",
        index=False,
    )


def __main__():
    annot_type_list = ["func-" + qtls, "prox", "hybrid-" + qtls]
    N_GENE_THRESH = 2
    SIGNIF_THRESH = 0.05
    ATC_LEVELS = get_atc_levels()
    for annot_type in annot_type_list:
        summarise_drug_results(
            gwas, annot_type, ATC_LEVELS, N_GENE_THRESH, SIGNIF_THRESH
        )
        summarise_gopath_results(gwas, annot_type, N_GENE_THRESH, SIGNIF_THRESH)


if __name__ == "__main__":
    __main__()


def get_annot_results(annot, gwas):
    file = os.path.join(
        magma_results_dir,
        gwas,
        "results",
        "magma_gene_result-" + gwas + "-" + annot + run_id + ".genes.out",
    )
    df = None

    if os.path.exists(file):
        df = pd.read_csv(file, delim_whitespace=True)
    else:
        display("File not found: " + file)

    return df


# All ATC genesets
atc_genesets_df = pd.concat(
    [
        pd.read_csv(
            data_path + "drugbank/dbank_gene_set-atc" + str(x + 1) + ".txt",
            header=None,
            sep="\t",
            index_col=0,
        )
        for x in range(5)
    ]
)
# All GO genesets
go_genesets_df = pd.read_csv(
    data_path + "/magma_analysis/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt",
    header=None,
    sep="\t",
    index_col=0,
)
# All GO + ATC genesets
genesets_df = pd.concat([atc_genesets_df, go_genesets_df])
# Emsembl to HGNC translation table
gene_hgnc_df = pd.read_csv(
    data_path + "/wrangling/NCBI37.3.ensembl.gene.loc",
    sep="\t",
    header=None,
    names=["Ensembl", "1", "chrom", "chromStart", "chromEnd", "HGNC"],
)[["Ensembl", "HGNC"]]
# Hybrid (func + prox) gene analysis results
hybriddf = get_annot_results("hybrid-" + qtls, gwas)
# Functional gene analysis results
funcdf = get_annot_results("func-" + qtls, gwas)
# Proximity gene analysis results
proxdf = get_annot_results("prox", gwas)
# -

# ## Annotated functional / hybrid  genes
#
# In all 20,000 genes were annotated using the hybrid method and functional brain QTLs, with a range of SNPs mapped to each gene numbering between 2 and 14,092. Of the 20,000 genes, 15,952 were annotated using functional data and the remaining 4,048 annotated using positional data. In comparison, the standard MAGMA proximity annotation contains 19,863 genes with SNP number ranges between 2 and 14,849.
#
# ## ATC gene ontology database
#
# The number of entries at each ATC level 1-5 are listed in table \ref{tab:atc_levels} below, with the corresponding number that have gene annotation data following curation, used as input into the MAGMA gene set analysis. The range of gene set sizes are also listed for reference, as there are known biases related to the size of a gene set<citet data-cite="mooney_gene_2015"><sup>mooney</sup></citet>, the common range limit being 10-200.
#
#  |ATC level       |1    |2    |3    |4    |5|
#  | --- | --- | --- | --- | --- | --- |
#  |Defined codes  |14   |93  |267  |885 |4823|
#  |Annotated      |14   |83  |181  |318  |473|
#  |Gene set size range |49-668|5-226|5-180|5-149|5-124|

display(Latex(r"\caption{Table \ref{tab:atc_levels} - ATC level annotation}"))
display(Latex(r"\label{tab:atc_levels}"))
display(Latex(r"\end{longtable}%"))

#  Removing MHC genes resulted in the gene *ENSG00000198704* being removed from the analysis as the only difference, this only affected the drug with ATC level 5 code *V03AB32* (an antidote and antioxidant), which along with it's constituent ATC codes for levels 1-4, was not present in any significant results obtained. This gene and drug were therefore not removed from the analysis.

# ## GWAS gene association
#
# For each of the three annotations - hybrid, functional and positional, a gene analysis was run to determine the effect of the hybrid method over the functional and to compare both with the standard positional approach to annotating genes.
#
# Below are the p-value histograms<citet data-cite="breheny_p-value_2018"><sup>Breheny</sup></citet> generated for these showing a distinct peak near zero of significant non-null features for all three:
#

# +
# Plotting functions
def plot_pval_histogram(gwas_annot, df):
    if df is not None:
        plt.subplot(fig_height, fig_width, plot_index, sharex=ax1, sharey=ax1)
        plt.hist(df["P_JOINT"], bins=55)
        plt.xlabel("p-value")
        plt.ylabel("Frequency")
        plt.title(gwas_annot + " gene p-value histogram")
        return plot_index + 1
    else:
        return plot_index


def plot_scatter_annotations(gwas, annot1, annot1df, annot2, annot2df):
    # merge frames on SET, giving P_left and P_right, only where the SET value is in both
    merged = pd.merge(annot1df, annot2df, how="inner", on="GENE")[
        ["GENE", "P_JOINT_x", "P_JOINT_y"]
    ]
    merged.columns = ["GENE", annot1 + " -log(Pgene)", annot2 + " -log(Pgene)"]

    # take -log10 of Pgene
    merged[annot1 + " -log(Pgene)"] = -np.log10(merged[annot1 + " -log(Pgene)"])
    merged[annot2 + " -log(Pgene)"] = -np.log10(merged[annot2 + " -log(Pgene)"])

    return merged


def plot_scatter_correlations(gwas, annot1, df1, annot2, df2):
    df = plot_scatter_annotations(gwas, annot1, df1, annot2, df2)

    # linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df[annot1 + " -log(Pgene)"], df[annot2 + " -log(Pgene)"]
    )

    # plotting
    sp = plt.subplot(fig_height, fig_width, plot_index)
    sp.scatter(df[annot1 + " -log(Pgene)"], df[annot2 + " -log(Pgene)"])
    sp.plot([0, 140], [0, 140], "r--", label="Equivalent p-value")
    plt.xlabel(annot1 + " -log(Pgene)")
    plt.ylabel(annot2 + " -log(Pgene)")
    plt.title(gwas + " " + annot1 + "/" + annot2 + " correlation")
    plt.legend(["Equivalent p-value"])
    ax2 = plt.gca()
    plt.text(
        0,
        -0.25,
        "r_value = " + "{0:.3f}".format(r_value) + ", err " + "{0:.3f}".format(std_err),
        transform=ax2.transAxes,
    )

    return plot_index + 1


# -

display(Latex(r"\vspace{-2mm}"))

# +
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

plot_index = 1
ax1 = plt.gca()
plot_index = plot_pval_histogram("PD hybrid", hybriddf)
plot_index = plot_pval_histogram("PD func", funcdf)
plot_index = plot_pval_histogram("PD prox", proxdf)
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-18mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd_histograms} - Gene analysis p-value histograms for association with Parkinson's disease, for experiments run using proximity annotated genes, functionally annotated genes and a hybrid of the two.}"
    )
)
display(Latex(r"\label{fig:pd_histograms}"))
display(Latex(r"\end{figure}"))
display(Latex(r"\newpage"))

# But although the positional histogram peak is higher than the functional, indicating more highly associated genes explainable partly by the difference in annotated gene numbers between the two, the hybrid results show a higher peak again, justifying it's use as a more biologically significant and potentially more powerful dataset in further analysis.
#
# The correlations between the three sets of results are of interest, as they also highlight the differences in p-value for the same gene between pairs of results:

# +
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

plot_index = plot_scatter_correlations(gwas, "hybrid", hybriddf, "prox", proxdf)
plot_index = plot_scatter_correlations(gwas, "func", funcdf, "prox", proxdf)
plot_index = plot_scatter_correlations(gwas, "hybrid", hybriddf, "func", funcdf)
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd_correlation} - Gene analysis p-value correlations between experiments run using proximity annotated genes, functionally annotated genes and a hybrid of the two.}"
    )
)
display(Latex(r"\label{fig:pd_correlation}"))
display(Latex(r"\end{figure}"))

# +
# Calculate how many genes are below the red equivalence line, to see if func approach gives higher p-vals
small_hybrid_df = hybriddf[["GENE", "ZSTAT"]]
small_func_df = funcdf[["GENE", "ZSTAT"]]
small_prox_df = proxdf[["GENE", "ZSTAT"]]

mrgd = pd.merge(small_hybrid_df, small_prox_df, on="GENE", how="inner")
total = mrgd.count()
better_func = mrgd[mrgd["ZSTAT_x"] > mrgd["ZSTAT_y"]].count()
equal = mrgd[mrgd["ZSTAT_x"] == mrgd["ZSTAT_y"]].count()
# display("Percentage of genes with higher hybrid p-vaules than prox " + str(better_func['GENE'] / total['GENE']))
# -

# Counting the proportion of genes below the line, for which the functional p-value of that gene is lower than it's proximity gene p-value indicates an improvement, giving a ratio of 53%.
#
# ### Choice of gene analysis parameters
#
# During subsequent runs of the processing pipeline as part of the development process it was quickly noticed that gene set results would vary considerably with no significant changes to the input data or pipeline itself, sometimes the difference being between no significant results and more than 10 significant gene sets on repeat runs. Noting that a default adaptive permutation MAGMA setting is involved when using the *multi=snp-wise* joint model, and that the MAGMA manual v1.06 alludes to a statistical effect of the fixed permutation setting with:
#
#  '*...this simple approach is not very efficient however, since a large number of permutations only has added value if the p-value is very low...*'
#
# without being specific as to the meaning of 'low' or 'added value', an exploratory analysis of the effect of the permutation setting on the result was undertaken.
#
# The *gene-settings* MAGMA flag allows two modifiers that affect the computation of an empirical p-value which is amalgamated subsequently into a gene's Z-score, and thus affecting the gene set analysis results. These modifiers are *adap-permp* and *fixed-permp*. *adap-permp* is the default taking optional maximum permutations, minimum permutations and stopping criteria parameters; whilst *fixed-permp* takes an optional number of permutations to run.
#
# Initially, 30 runs of an identical pipeline with invariant inputs were run using the default adaptive permutation setting on the (2014) Parkinson's disease GWAS. This scenario was then repeated three more times, only varying the permutation settings to 5,000, 10,000 and 20,000 fixed permutations subsequently for each set of 30 runs. The results, consisting of the set of p-values obtained for each run, expressed as a negative log of the gene set p-value are shown below as box plots for three significant gene sets (J06, V10ZA and L04AA23) and three non-significant gene sets (C03, V03AB and A01AB09).
#
# #### Significant gene sets
#

# + {"caption": "somecaption", "label": "fig:somelabel", "widefigure": true}
perms = ["adapt", "5k", "10k", "20k"]
idx = range(0, 30)


def plot_adapt_fixed_boxplot(setname, gwas, plotnum, min_ylim, max_ylim):
    set_df = pd.DataFrame(index=idx, columns=perms)
    for perm in perms:
        set_df[perm] = pd.read_csv(
            "~/source/magma-perms-analysis/"
            + gwas
            + "/"
            + "magma-perms-"
            + gwas
            + "-"
            + perm
            + "-"
            + setname
            + ".txt",
            header=None,
            delim_whitespace=True,
        )
    set_df = set_df.apply(np.log10).abs()

    fig.add_subplot(fig_height, fig_width, plotnum)
    graph = sns.boxplot(data=set_df).set(
        xlabel="Permutation setting",
        ylabel="-log(p-value)",
        title=gwas + " p-value variation for gene set " + setname,
        ylim=(min_ylim, max_ylim),
    )
    return plotnum + 1


fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("J06", gwas, plotnum, 1, 12)
plotnum = plot_adapt_fixed_boxplot("V10XA", gwas, plotnum, 1, 12)
plotnum = plot_adapt_fixed_boxplot("L04AA23", gwas, plotnum, 1, 12)
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd_sigperms} - Comparison of the Parkinson's disease association -log(p-value) variance between different gene analysis permutation settings for three ATC codes significantly associated with Parkinson's.}"
    )
)
display(Latex(r"\label{fig:pd_sigperms}"))
display(Latex(r"\end{figure}"))

# The box plots in figure \ref{fig:pd_sigperms} show a marked discrepancy between the adaptive permutation and other fixed permutation results - in the case of ATC gene set *J06*, the range of values are almost exclusive with very little overlap.
#
# Looking at the QQ plots for ATC gene set *J06* in figure \ref{fig:pd_qq_sigperms}, it appears that the discrepancy between the adaptive and fixed permutation results is not due to skewed distributions, as they are approximately normal.

# +
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

sets_df = pd.DataFrame(index=idx, columns=perms)
for perm in perms:
    sets_df[perm] = pd.read_csv(
        "~/source/magma-perms-analysis/PD/magma-perms-PD-" + perm + "-J06.txt",
        header=None,
        delim_whitespace=True,
    )
    sets_df[perm] = sets_df[perm].apply(np.log10).abs()

ax = fig.add_subplot(fig_height, fig_width, 1)
stats.probplot(sets_df["adapt"], dist="norm", plot=pylab)
ax.set_title("PD/J06 Adaptive p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 2)
stats.probplot(sets_df["10k"], dist="norm", plot=pylab)
ax.set_title("PD/J06 10k p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 3)
stats.probplot(sets_df["20k"], dist="norm", plot=pylab)
ax.set_title("PD/J06 20k p-value QQ plot against Normal")
pylab.show()
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd_qq_sigperms} - Comparison of QQ plots of various permutation settings for Parkinson's disease association with ATC gene set J06. \normalfont{Each blue dot represents a p-value association of the gene set with the disease resulting from one gene analysis run. The red line indicates the normal distribution.}}"
    )
)
display(Latex(r"\label{fig:pd_qq_sigperms}"))
display(Latex(r"\end{figure}"))

# \newpage
# #### Non-significant gene sets

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("C03", gwas, plotnum, 0, 0.56)
plotnum = plot_adapt_fixed_boxplot("V03AB", gwas, plotnum, 0, 0.56)
plotnum = plot_adapt_fixed_boxplot("A01AB09", gwas, plotnum, 0, 0.56)

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd_notsigperms} - Comparison of the Parkinson's disease association -log(p-value) variance between different gene analysis permutation settings for three ATC codes not significantly associated with Parkinson's.}"
    )
)
display(Latex(r"\label{fig:pd_notsigperms}"))
display(Latex(r"\end{figure}"))

# Comparing the significant (figure \ref{fig:pd_sigperms}) and non-significant box plots (figure \ref{fig:pd_notsigperms}), the results show that for gene set p-values, an effect was observed where between the default adaptive permutation setting and fixed permutation settings, for multiple runs of the same dataset there is an increased variance of the adaptive results and higher mean (and therefore lower -log(mean)) compared to fixed permutation values > 5,000. The separation of the adaptive and fixed resultant -log(p-value) ranges is greater for more significant gene sets.
#
# To query if the discrepancy id due to the properties of the gene sets used, the experiment was run again but this time instead of using the ATC defined gene sets a set of 500 randomly chosen gene sets of size 10-1000 were used, where the constituent genes are selected randomly from the ATC set used previously. Three of the gene set results are shown in figure \ref{fig:pd_randperms} for comparison.
#
# #### 2017 meta-gwas PD analysis

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("J06", "PD2", plotnum, 2.5, 10)
plotnum = plot_adapt_fixed_boxplot("V10XA", "PD2", plotnum, 2.5, 10)
plotnum = plot_adapt_fixed_boxplot("L04AA23", "PD2", plotnum, 2.5, 10)

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd2_notsigperms} - Comparison of the 2017 Parkinson's disease GWAS meta study association -log(p-value) variance between different gene analysis permutation settings for three ATC codes not significantly associated with Parkinson's.}"
    )
)
display(Latex(r"\label{fig:pd2_notsigperms}"))
display(Latex(r"\end{figure}"))

# +
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

sets_df = pd.DataFrame(index=idx, columns=perms)
for perm in perms:
    sets_df[perm] = pd.read_csv(
        "~/source/magma-perms-analysis/PD2/magma-perms-PD2-" + perm + "-J06.txt",
        header=None,
        delim_whitespace=True,
    )
    sets_df[perm] = sets_df[perm].apply(np.log10).abs()

ax = fig.add_subplot(fig_height, fig_width, 1)
stats.probplot(sets_df["adapt"], dist="norm", plot=pylab)
ax.set_title("PD2/J06 Adaptive p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 2)
stats.probplot(sets_df["10k"], dist="norm", plot=pylab)
ax.set_title("PD2/J06 10k p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 3)
stats.probplot(sets_df["20k"], dist="norm", plot=pylab)
ax.set_title("PD2/J06 20k p-value QQ plot against Normal")
pylab.show()
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd2_qq_sigperms} - Comparison of QQ plots of various permutation settings for 2017 Parkinson's disease GWAS meta study association with ATC gene set J06. \normalfont{Each blue dot represents a p-value association of the gene set with the disease resulting from one gene analysis run. The red line indicates the normal distribution.}}"
    )
)
display(Latex(r"\label{fig:pd2_qq_sigperms}"))
display(Latex(r"\end{figure}"))

# #### Randomly allocated gene sets

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("set_3_177", gwas, plotnum, 0.5, 5.5)
plotnum = plot_adapt_fixed_boxplot("set_92_218", gwas, plotnum, 0.5, 5.5)
plotnum = plot_adapt_fixed_boxplot("set_445_499", gwas, plotnum, 0.5, 5.5)

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:pd_randperms} - Comparison of the Parkinson's disease association -log(p-value) variance between different gene analysis permutation settings for three random genesets.}"
    )
)
display(Latex(r"\label{fig:pd_randperms}"))
display(Latex(r"\end{figure}"))
display(Latex(r"\vspace{-10mm}"))

# +
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

sets_df = pd.DataFrame(index=idx, columns=perms)
for perm in perms:
    sets_df[perm] = pd.read_csv(
        "~/source/magma-perms-analysis/PD/magma-perms-PD-" + perm + "-set_445_499.txt",
        header=None,
        delim_whitespace=True,
    )
    sets_df[perm] = sets_df[perm].apply(np.log10).abs()

ax = fig.add_subplot(fig_height, fig_width, 1)
stats.probplot(sets_df["adapt"], dist="norm", plot=pylab)
ax.set_title("PD/set_445_499 Adaptive p-value QQ plot vs Normal")
ax = fig.add_subplot(fig_height, fig_width, 2)
stats.probplot(sets_df["10k"], dist="norm", plot=pylab)
ax.set_title("PD/set_445_499 10k p-value QQ plot vs Normal")
ax = fig.add_subplot(fig_height, fig_width, 3)
stats.probplot(sets_df["20k"], dist="norm", plot=pylab)
ax.set_title("PD/set_445_499 20k p-value QQ plot vs Normal")
pylab.show()
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-15mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:rand_qq_sigperms} - Comparison of QQ plots of various permutation settings for Parkinson's disease association with random gene set set\_445\_499. \normalfont{Each blue dot represents a p-value association of the gene set with the disease resulting from one gene analysis run. The red line indicates the normal distribution.}}"
    )
)
display(Latex(r"\label{fig:rand_qq_sigperms}"))
display(Latex(r"\end{figure}"))
display(Latex(r"\newpage"))

# These show the larger variance for the adaptive setting, and tending towards a difference in means for more significant results as was seen with the ATC data set. This variation range in some cases also straddles the -log(p-value) of 1.3, which corresponds to a p-value of 5% - the significance threshold chosen. This gives evidence to the speculation that the effect is due to the MAGMA implementation, as opposed to an artifact of specific real-world input data sets.

# #### Comaprison with other GWAS results - Huntington's
#
# The ATC genes sets were also analysed against different gene-annotated GWAS results sets, in this first case for Huntington's disease (HD). This is a lower powered study, with a small, questionable p-value peak near zero for the three annotation approaches used previously for Parkinson's disease:

# +
run_id = ""
hd_hybriddf = get_annot_results("hybrid-" + qtls, "HD")
hd_funcdf = get_annot_results("func-" + qtls, "HD")
hd_proxdf = get_annot_results("prox", "HD")
run_id = "-no_kaviar"

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

plot_index = 1
ax1 = plt.gca()
plot_index = plot_pval_histogram("HD hybrid", hd_hybriddf)
plot_index = plot_pval_histogram("HD func", hd_funcdf)
plot_index = plot_pval_histogram("HD prox", hd_proxdf)
# -

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:hd_histograms} - Gene analysis p-value histograms for association with Huntington's disease, for experiments run using proximity annotated genes, functionally annotated genes and a hybrid of the two.}"
    )
)
display(Latex(r"\label{fig:hd_histograms}"))
display(Latex(r"\end{figure}"))

# And the boxplots comparing p-value variation for the adaptive, 5k, 10k and 20k fixed settings for three of the more significantly associated ATC gene sets (M05, V03AF and A11AA02) are:

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("M05", "HD", plotnum, 1.5, 9)
plotnum = plot_adapt_fixed_boxplot("V03AF", "HD", plotnum, 1.5, 9)
plotnum = plot_adapt_fixed_boxplot("A11AA02", "HD", plotnum, 1.5, 9)

display(Latex(r"\vspace{-15mm}"))
display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(
    Latex(
        r"\caption{Figure \ref{fig:hd_perms} - Comparison of the Huntington's disease association -log(p-value) variance between different gene analysis permutation settings for three random genesets.}"
    )
)
display(Latex(r"\label{fig:hd_perms}"))
display(Latex(r"\end{figure}"))

# These show no real difference in means between the fixed and adaptive approaches, but the difference in variance is still present.
#
# #### Comaprison with other GWAS results - Schizophrenia
#
# The Schitzophrenia GWAS is a much more high powered analysis, as can be seen by the p-value histograms in figure \ref{fig:sz_histograms}. In this case, the adaptive results look far more acceptible, with only a slight increase in variation when comparing with the fixed permutation results.

# +
# << TODO >> Look at SZ
run_id = "-1_20k-prox"
sz_hybriddf = get_annot_results("hybrid-" + qtls, "SZ")
sz_funcdf = get_annot_results("func-" + qtls, "SZ")
sz_proxdf = get_annot_results("prox", "SZ")
run_id = "-no_kaviar"

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

plot_index = 1
ax1 = plt.gca()
plot_index = plot_pval_histogram("SZ hybrid", sz_hybriddf)
plot_index = plot_pval_histogram("SZ func", sz_funcdf)
plot_index = plot_pval_histogram("SZ prox", sz_proxdf)
# -

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:sz_histograms} - Gene analysis p-value histograms for association with Schizophrenia, for experiments run using proximity annotated genes, functionally annotated genes and a hybrid of the two.}"
    )
)
display(Latex(r"\label{fig:sz_histograms}"))
display(Latex(r"\end{figure}"))

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("C08", "SZ", plotnum, 8, 23)
plotnum = plot_adapt_fixed_boxplot("C08DA", "SZ", plotnum, 8, 23)
plotnum = plot_adapt_fixed_boxplot("C08DA51", "SZ", plotnum, 8, 23)

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:sz_sigperms} - Comparison of the Schizophrenia association -log(p-value) variance between different gene analysis permutation settings for three ATC codes significantly associated with Schizophrenia.}"
    )
)
display(Latex(r"\label{fig:sz_sigperms}"))
display(Latex(r"\end{figure}"))
display(Latex(r"\newpage"))

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("N05", "SZ", plotnum, 0.4, 12)
plotnum = plot_adapt_fixed_boxplot("N05CD", "SZ", plotnum, 0.4, 12)
plotnum = plot_adapt_fixed_boxplot("N05CD06", "SZ", plotnum, 0.4, 12)

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:sz_notsigperms} - Comparison of the Schizophrenia association -log(p-value) variance between different gene analysis permutation settings for three ATC codes not significantly associated with Schizophrenia.}"
    )
)
display(Latex(r"\label{fig:sz_notsigperms}"))
display(Latex(r"\end{figure}"))

# #### SZ PGC1
#
# Comparing against the SZ drug geneset analysis for the PGC1 GWAS with proximity mapped genes, again 30 runs are used.

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("C08", "SZ1", plotnum, 0, 7)
plotnum = plot_adapt_fixed_boxplot("C08DA", "SZ1", plotnum, 0, 7)
plotnum = plot_adapt_fixed_boxplot("C08DA51", "SZ1", plotnum, 0, 7)

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:sz1_sigperms} - Comparison of the Schizophrenia (PGC1) association -log(p-value) variance between different gene analysis permutation settings for three ATC codes most significantly associated with Schizophrenia.}"
    )
)
display(Latex(r"\label{fig:sz1_sigperms}"))
display(Latex(r"\end{figure}"))
display(Latex(r"\newpage"))

fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(wspace=0.2)
plotnum = 1
plotnum = plot_adapt_fixed_boxplot("N05", "SZ1", plotnum, 0, 4)
plotnum = plot_adapt_fixed_boxplot("N05CD", "SZ1", plotnum, 0, 4)
plotnum = plot_adapt_fixed_boxplot("N05CD06", "SZ1", plotnum, 0, 4)

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-10mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:sz1_notsigperms} - Comparison of the (PGC1) Schizophrenia association -log(p-value) variance between different gene analysis permutation settings for three ATC codes not significantly associated with Schizophrenia.}"
    )
)
display(Latex(r"\label{fig:sz1_notsigperms}"))
display(Latex(r"\end{figure}"))

# +
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.5, wspace=0.2)

sets_df = pd.DataFrame(index=idx, columns=perms)
for perm in perms:
    sets_df[perm] = pd.read_csv(
        "~/source/magma-perms-analysis/SZ1/magma-perms-SZ1-" + perm + "-C08DA.txt",
        header=None,
        delim_whitespace=True,
    )
    sets_df[perm] = sets_df[perm].apply(np.log10).abs()

ax = fig.add_subplot(fig_height, fig_width, 1)
stats.probplot(sets_df["adapt"], dist="norm", plot=pylab)
ax.set_title("SZ1/C08DA Adaptive p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 2)
stats.probplot(sets_df["10k"], dist="norm", plot=pylab)
ax.set_title("SZ1/C08DA 10k p-value QQ plot against Normal")
ax = fig.add_subplot(fig_height, fig_width, 3)
stats.probplot(sets_df["20k"], dist="norm", plot=pylab)
ax.set_title("SZ1/C08DA 20k p-value QQ plot against Normal")
pylab.show()
# -

display(Latex(r"\begin{figure}[h]"))
display(Latex(r"\centering"))
display(Latex(r"\vspace{-15mm}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:sz1_qq_sigperms} - Comparison of QQ plots of various permutation settings for Schizophrenia association with ATC gene set C08DA. \normalfont{Each blue dot represents a p-value association of the gene set with the disease resulting from one gene analysis run. The red line indicates the normal distribution.}}"
    )
)
display(Latex(r"\label{fig:sz1_qq_sigperms}"))
display(Latex(r"\end{figure}"))
display(Latex(r"\newpage"))

# ## Gene set analysis results
#
# ### Significant ATC associations identified
#
# Table \ref{tab:pd_genesetresults} below shows the 26 significantly associated ATC gene sets that survive correction for multiple testing, having a q-value greater than 0.05. Each is listed with the number of genes contributing to the set in Ensembl notation along with the equivalent HGNC gene name and it's z-score of association with Parkinson's. The text '-?-' is displayed either if the HGNC gene name is not available or the gene was not annotated (in which case it did not take part in the gene set analysis). Genes are listed in descending Z-score order for each ATC gene set, the Z-score giving the strength of the gene's association with Parkinson's disease.
#
# \scriptsize

# +
def add_genesets(table):
    table2 = pd.DataFrame(columns=table.columns)
    for index, row in table.iterrows():
        aset = genesets_df.loc[index].str.split(" ").tolist()
        setgenes_df = pd.DataFrame(aset).melt()[["value"]]
        setgenes_df.columns = ["Ensembl"]
        if not setgenes_df.empty:
            setgenes2_df = pd.merge(
                setgenes_df, hybriddf, how="left", left_on="Ensembl", right_on="GENE"
            )[["Ensembl", "ZSTAT"]]
            hugo_df = pd.merge(setgenes2_df, gene_hgnc_df, how="left", on="Ensembl")

            # Limit to top 50 for now
            hugo_df = hugo_df.head(50)

            hugo_df["absZSTAT"] = hugo_df["ZSTAT"].abs()
            hugo_df.sort_values(by=["absZSTAT"], ascending=False, inplace=True)
            hugo_df.drop("absZSTAT", axis=1, inplace=True)
            hugo_df["ZSTAT"] = hugo_df["ZSTAT"].round(5)
            hugo_df.fillna("-?-", inplace=True)
            hugo_df.reset_index(inplace=True)

            table2 = table2.append(row)
            table2.at[index, "Ensembl"] = hugo_df["Ensembl"][0]
            table2.at[index, "HGNC"] = hugo_df["HGNC"][0]
            table2.at[index, "ZSTAT"] = hugo_df["ZSTAT"][0]
            table2.at[index, "SET"] = index
            table2 = table2.append(hugo_df[1:])

    # Hack for gosets
    if "SHORT_NAME" in table2.columns:
        table2.rename({"SET": "Description", "SHORT_NAME": "SET"}, axis=1, inplace=True)

    table2.fillna(" ", inplace=True)
    table2 = table2.reindex_axis(
        ["SET", "NGENES", "P", "Q", "Description", "Ensembl", "HGNC", "ZSTAT"], axis=1
    )

    return table2


# Display a set of results tables for the given GWAS list, file type and column (P or Q)
def display_drug_tables(gwas, file_postfix, column, table_id):
    summary_dir = os.path.join(summary_results_dir)
    pval_files = [
        file for file in os.listdir(summary_dir) if file.endswith(gwas + file_postfix)
    ]
    pval_files.sort()

    table = pd.read_csv(os.path.join(summary_dir, pval_files[table_id]), sep="\t")
    if not table.empty:
        table.sort_values(column, inplace=True)
        table.set_index("SET", inplace=True)
        table["Q"] = table["Q"].round(6)
        table["P"] = table["P"].round(6)
        if len(table) != 0:
            table = add_genesets(table)
    return table


# + {"caption": "somecaption", "label": "fig:somelabel", "widefigure": true}
table = display_drug_tables(
    gwas, run_id + "-" + "hybrid-" + qtls + "_qvals.tsv", "Q", 0
)
display(
    Latex(
        table.to_latex(
            column_format="p{2cm}p{1cm}p{1cm}p{0.8cm}p{5cm}cp{1cm}c",
            multirow=True,
            multicolumn=True,
            index=False,
        )
    )
)
display(
    Latex(
        r"\caption{Table \ref{tab:pd_genesetresults} - Significant Parkinson's ATC gene sets}"
    )
)
display(Latex(r"\label{tab:pd_genesetresults}"))
display(Latex(r"\end{longtable}%"))
table.to_csv("~/table.csv")
# -

# \normalsize
#
#
# The union of the 26 significant gene sets comprising of 39 genes are listed in table \ref{tab:pd_genesetgenes} below, ordered by Z-score of the gene's association with Parkinson's disease.
#
# \scriptsize

# +
geneset = set(table["Ensembl"])
tmp = pd.DataFrame(list(geneset), columns=["Ensembl"])
tmp2 = table.drop_duplicates(subset=["Ensembl"])

result = pd.merge(tmp, tmp2, how="left", validate="1:1")
result = result.reindex_axis(["Ensembl", "HGNC", "ZSTAT"], axis=1)

result["absZSTAT"] = pd.to_numeric(result["ZSTAT"], errors="coerce").abs()
result.sort_values(by=["absZSTAT"], ascending=False, inplace=True)
result.drop("absZSTAT", axis=1, inplace=True)
result.fillna("-?-", inplace=True)
# -

display(
    Latex(
        result.to_latex(
            column_format="lll", multirow=True, multicolumn=True, index=False
        )
    )
)
display(
    Latex(
        r"\caption{Table \ref{tab:pd_genesetgenes} - Union of genes of ATC code genesets showing significant association with Parkinson's}"
    )
)
display(Latex(r"\label{tab:pd_genesetgenes}"))
display(Latex(r"\end{longtable}%"))

# \normalsize
# ### Significant GO pathway associations identified
#
# Table \ref{tab:pd_pathwayresults} below shows the 12 significantly associated GO term gene sets that survive correction for multiple testing, having a q-value greater than 0.05. Each is listed with the top 50 genes contributing to the set in Ensembl notation along with the equivalent HGNC gene name and it's z-score of association with Parkinson's. The text '-?-' is displayed either if the HGNC gene name is not available or the gene was not annotated (in which case it did not take part in the gene set analysis). Genes are listed in descending Z-score order for each gene set, the Z-score giving the strength of the gene's association with Parkinson's disease.
#
# \scriptsize

# +
# pathway table
pwtable = display_drug_tables(
    gwas, run_id + "-" + "hybrid-" + qtls + "_qvals.tsv", "Q", 1
)
# display(table)
pwtable["Description"] = pwtable["Description"].str.replace("_", " ")
pwtable["Description"] = pwtable["Description"].str.replace("^GO.........", "")
pwtable["SET"] = pwtable["SET"].str[:10]

display(
    Latex(
        pwtable.to_latex(
            column_format="p{2cm}p{1cm}p{1cm}p{0.8cm}p{5cm}cp{1cm}c",
            multirow=True,
            multicolumn=True,
            index=False,
        )
    )
)
display(
    Latex(
        r"\caption{Table \ref{tab:pd_pathwayresults} - Significant Parkinson's pathway gene sets}"
    )
)
display(Latex(r"\label{tab:pd_pathwayresults}"))
display(Latex(r"\end{longtable}%"))
# -

# \newpage
# \normalsize
#
# # Discussion
#
# The association of a set of existing approved drugs with Parkinson's disease was investigated using functional gene set analysis. First, three gene annotation sets were curated using currently available up-to-date data sources - one positionally annotated, sourced from the MAGMA tool auxiliary data, a second functionally annotated gene set generated using several brain-specific QTL data sets and a third that is a hybrid of the two to allow maximum coverage of known genes. These annotation sets were then used in three separate parallel gene analyses, one for each set, to provide association p-values with Parkinson's disease. Next the ATC, a drug classification system, was used to generate a group of gene sets stratified into 5 levels containing gene drug targets, according to the 5 levels of the ATC. These were in turn input to a final gene set analysis step against the previous three annotation datasets - positional, functional and hybrid. An additional gene set analysis was performed for all three annotations against a group of gene sets obtained from Gene Ontology (GO) terms.
#
# During the gene analysis stage some extra exploratory work was carried out in order to overcome perceived issues in the default MAGMA settings for the chosen gene analysis model. This resulted in the default top-1 gene adaptive permutation setting being discarded from the main analysis pipeline, replaced by a fixed permutation of 50,000 in order to obtain a more consistent simulated p-value for the top-1 gene.
#
# With 20,000 genes having an annotation for the final hybrid annotation set, this is within the current best estimate limit of 19,000 - 20,000 protein coding genes in the human genome<citet data-cite="ezkurdia_multiple_2014"><sup>mezkurdia</sup></citet> and a good outcome. But with just over 4,000 of them being annotated using positional data, more effort to increase the functionally described portion of genes could still provide an improvement.
#
# Only around 10% of the defined level 5 ATC codes were annotated sufficiently to pass quality control (473 out of 4823), using the Drugbank data set. Although 4823 drugs is the theoretical maximum, not all of these will have protein targets (the current proportion is not known without more detailed analysis) but it is thought this proportion could be improved upon by incorporating data from other sources.
#
# The gene analysis histograms given in figure \ref{fig:pd_histograms} indicate that the hybrid approach gives the better result, this approach also could have added extra noise and increased the false discovery rate - measuring the differences in q-values for these results may be a more reliable measurement for comparison. Would extending this comparison by performing multiple gene analyses of all combinations of available tissue type QTLs to find which would give the best association be an interesting exploratory experiment?
#
# Some complimentary drug classes were found to be highly associated with Parkinson's, as revealed by the 26 significant ATC code associations discovered; these are (taking their ATC level 2 descriptions as an overview), **J06** - *immune sera and immunoglobulins* and **L04** - *Immunosuppressants* both relating to the immune system. Several ATC codes were also indicated in the **V10** level 2 code, the signal for which appears to be mainly due to the drug Tositumomab (**V10AX53**), an iodine based radioisotope used to treat non-Hodgkin's lymphoma - but the drugbank also lists two related iodine isotopes that are used in the treatment and diagnosis of Parkinson's - Altropane and Ioflupane I-123.
#
# For these drug-disease relationships, it is not determined in which direction the drug affects the disease, whether it exacerbates or alleviates the condition. One broadly useful measure for this could be the change of expression of the genes involved - generally an increase in expression could be an indicator of a drug useful for fighting disease.
#
#
#
# Discuss the connection between these genes (via pathways?), suggest network analysis as the next step that does this in a more automated way - a tool is needed to do this?

# +
# HTML('''<script>
# code_show=true;
# function code_toggle() {
# if (code_show){
# $('div.input').hide();
# } else {
# $('div.input').show();
# }
# code_show = !code_show
# }
# $( document ).ready(code_toggle);
# </script>
# The raw code for this IPython notebook is by default hidden for easier reading.
# To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>.''')
# -

gitspec = subprocess.check_output(
    ["git", "describe", "--all", "--dirty", "--long"]
).strip()
# print('This notebook was run on git revision ' + str(gitspec, 'utf-8'))
