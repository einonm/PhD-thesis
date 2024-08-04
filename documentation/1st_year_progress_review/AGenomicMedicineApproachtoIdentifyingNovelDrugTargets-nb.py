# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
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

# !pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt

# !pip install scipy
module_path = os.path.abspath(os.path.join("../..", "summaries"))
if module_path not in sys.path:
    sys.path.append(module_path)
import qvalue
from IPython.display import display, Latex, display_latex

from scipy import stats
import warnings

warnings.filterwarnings("ignore")

# Also generate tables in latex format, with formatting opts
pd.set_option("display.latex.repr", True)
pd.set_option("display.latex.longtable", True)
pd.set_option("display.max_colwidth", 135)
# -

# \newpage
# # Introduction
# Studies have demonstrated that for genes which are linked to a disease by genome-wide association studies, if the proteins for which they encode are associated with targets for existing medications<citet data-cite="ruderfer_polygenic_2016"><sup>ruderfer</sup></citet> then it follows that these genes are clear targets for drug development studies<citet data-cite="stein_effect_2012,navarese_effects_2015"><sup>stein,navarese</sup></citet> and offer opportunities to re-purpose existing approved drugs. It has also recently been demonstrated that drugs supported by genetic associations are more likely to successfully advance through the later stages of clinical trials<citet data-cite="nelson_support_2015"><sup>nelson</sup></citet>.
# Genome-wide association studies (GWAS) have now collectively identified common susceptibility variants at a large number of loci that increase risk for neuropsychiatric and neurological disease. Large-scale exome and whole genome sequencing studies are also starting to identify low frequency genetic mutations that modify the functional properties of the encoded proteins<citet data-cite="olgiati_dnajc6mutations_2016,farlow_whole-exome_2016,lesage_loss_2016,jaberi_mutation_2016"><sup>olgiati,farlow,lesage,jaberi</sup></citet>. To evaluate the clinical impact of these susceptibility alleles, it is now essential to understand their biological significance. Such insight might improve patient diagnosis and prognosis, and offer the potential to develop novel or improved therapies.

# ## Aims of the research
# This project is based on the hypothesis that identifying drugs targeting proteins, which are associated with genes that increase risk to neuropsychiatric and neurological diseases will reveal novel therapeutic opportunities that can be used to inform and improve treatment of the diseases.
#
# To explore this hypothesis, an analytical pipeline has been designed that performs gene set analyses between a set of functionally annotated genes generated from GWAS summary statistics and another gene set obtained from a curated database of known drug gene targets. Any significant results from this analysis could be used to highlight novel opportunities for therapeutic interventions, and would be a candidate for further investigation using methods such as network-based proximity analysis, as used by Guney et al.<citet data-cite="guney_network-based_2016"><sup>guney</sup></citet> within the interactome of each drug set - as planned for exploration in the later stages of the project.
#
#
# This pipeline can be described further in three parts:
#
# ### Generate a drug targets/ontology database
# Using several independent and publicly available databases of drug data, a set of drugs and their drug protein/gene targets will be gathered, indexed according to the Anatomical Therapeutic Chemical (ATC) classification system<citet data-cite="noauthor_anatomical_nodate"><sup>ATC</sup></citet>. The genes encoding the protein targets will be established with reference to the human genome reference sequence (hg37). The ATC pharmaceutical coding system divides drugs according to the mechanism of action as well as their therapeutic and chemical characteristics. Drugs are classified according to 5 levels: 1) anatomical main group, 2) therapeutic main group, 3) therapeutic/pharmacological subgroup, 4) chemical/therapeutic/pharmacological subgroup, and 5) chemical substance.
#
# As an example, the level 5 ATC code *A11GA01* is described as *Ascorbic acid (vitamin C)*. This also belongs to four other pharmacological groups according to the ATC system, one for each of the other levels, for the levels 1-4 these are namely *A*, for *Alimentary tract and metabolism*, *A11* for *vitamins*, *A11G* for *Ascorbic acid (Vitamin C), incl. combinations* and finally *A11GA* for *Ascorbic acid (vitamin C), plain*.
#
# Throughout this document, the terms *pharmacological group* and *ATC code* are used interchangably for a group of drugs at any ATC level.
#
# ### Curate a genomic map of functional genetic variants
# A collection of available functional Quantative Trait Loci (QTL) mappings will be compiled, including from methylation QTLs and expression QTLs for a range of tissue types that functionally link SNPs to genes. This functional mapping can be merged with GWAS results for a particular disease by the MAGMA software, an established method of assigning SNPs to genes<citet data-cite="leeuw_magma:_2015"><sup>leeuw</sup></citet>, taking the aggregate p-value from both the mean SNP and top SNP gene analysis models, and fully accounts for LD between SNPs.
#
# ### Identify pharmacologically related sets of drug targets associated with disease
# The drug target and genetic association data will then be further analysed using the gene enrichment analysis capabilities of MAGMA, which will perform a competitive gene-set analysis to test whether the genes included in each pharmacologically related gene-set are more strongly associated with the disease than other genes. This approach includes conditional effects such as gene size in the analysis and will also account for the correlations between neighbouring genes caused by LD. As the gene-sets are being objectively defined by the ATC classification system an independent analysis will be conducted for each of the 5 ATC levels.
# The output gene-set association p-values for each gene-set analysis will be corrected for multiple testing by calculating q-values (Storey), which adjust the p-values based on an optimised False Discovery Rate (FDR) approach. A significance threshold of 5% is chosen, so any association with a q-value below a value of 0.05 is deemed to be significant.

# ## Analytical pipeline implementation and methods
#
# For results to be generated in a repeatable and reproducible manner, the implemented software pipeline and modifications must be recorded and the version referenced along with the results. To facilitate this, the git distributed version control tool<citet data-cite="noauthor_git_nodate"><sup>git</sup></citet> and Jupyter notebooks<citet data-cite="thomas_jupyter_2016"><sup>jupyter</sup></citet> are used. The widely used Python language is chosen to implement the bulk of the pipeline. All notebooks and scripts were run on local HPC cluster resources, and bash scripting used to orchestrate running code in batches and on the HPC cluster where applicable.
#
# A Jupyter notebook was created for each of the main pipeline processing steps (figure \ref{fig:toplevel}), performing the tasks of exploratory data analysis, implementing data transformations and displaying visualisations. Where applicable, notebooks were saved as python scripts to run on data in batches. Separate Jupyter notebooks were created for important data wrangling activities, such as gene ID conversions and assigning genomic positions to SNPs and genes.
#
# In order to provide additional validation checks of the functional gene annotation implementation, MAGMA was again used to generate a positionally annotated gene-set, where GWAS SNPs are mapped to genes according to their genomic location, with the closest gene to a SNP being mapped to it. Throughout this document, this position-based annotation will be referred to as a *proximity gene annotation* and the analytical pipeline, when run using proximity gene annotation data, as a *proximity gene analysis*. Any proximity gene analysis is run parallel to the functional gene analysis, differing only in the gene annotation file used to generate the GWAS gene p-values. This is indicated in grey on Figure \ref{fig:toplevel}.
#
# The pipeline was implemented and optimised to run on a local HPC cluster, giving a run time for producing a set of summary results for both functional and proximity based analysis on several GWAS summary datasets of around an hour.
#
# Ensembl gene IDs<citet data-cite="aken_ensembl_2016"><sup>ensembl</sup></citet> are used throughout the pipeline to reference genes and rs numbers to reference SNPs. All gene and SNP references are based on the hg37 genome reference sequence (also known as Build 37 or hg19).
#
# The major histocompatibility complex (MHC), or other large areas of genetic association were not removed from the GWAS results.

# Insert the toplevel pipeline diagram
display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(Latex(r"\includegraphics{toplevel-pipeline}"))
display(
    Latex(
        r"\caption{Figure \ref{fig:toplevel} - Overview of pipeline steps. \normalfont{Green boxes represent data processing activities. Cylinders are data inputs and outputs - usually blue, grey represents validation data and yellow the final output.}}"
    )
)
display(Latex(r"\label{fig:toplevel}"))
display(Latex(r"\end{figure}"))

# # Results
#
# ## Drug-target ATC indexed database
# As an initial development data set, only one drug-target source was used, the DrugBank<citet data-cite="law_drugbank_2014"><sup>drugbank</sup></citet>. The full DrugBank database was downloaded locally and for each drug that has protein targets and a valid ATC code the drug name, the protein targets and ATC code were extracted. Protein identifiers were converted to genes using UniProt<citet data-cite="noauthor_uniprot:_2017"><sup>uniprot</sup></citet>. ATC codes with less than 5 protein targets were discarded. ATC codes where another ATC code at the same ATC level has a matching target gene list were also removed. This analysis was repeated using the set of ATC codes at each of the ATC levels 1-5.
#
# The results are summarised in table \ref{tab:atc_db}, which also gives the total number of ATC codes filtered out at each stage.

# +
stages = [
    "Defined ATC codes",
    "ATC codes from DrugBank",
    "After removing ATC codes < 5 genes",
    "After removing duplicate gene lists (final tally)",
]

df = pd.DataFrame(
    {
        "ATC level": range(1, 6),
        stages[0]: [14, 93, 267, 885, 4823],
        stages[1]: [14, 87, 203, 506, 1560],
        stages[2]: [14, 75, 133, 212, 348],
        stages[3]: [14, 75, 133, 210, 315],
    }
)

df.set_index("ATC level", inplace=True)
df = df.transpose()
df = df.reindex([stages[0], stages[1], stages[2], stages[3]])
display(Latex(df.to_latex(column_format="p{8cm}rrrrr")))
# -

display(
    Latex(
        r"\caption{Table \ref{tab:atc_db} - Table of ATC code tallies during drug-target / ATC database creation.}"
    )
)
display(Latex(r"\label{tab:atc_db}"))
display(Latex(r"\end{longtable}%"))

# ## Functional gene annotation dataset
# To generate a combined functional gene annotation dataset over all QTLs, each individual QTL dataset is analysed separately, and the resultant set of annotations combined into one. The GTEx expression QTL data sets<citet data-cite="carithers_novel_2015"><sup>gtex</sup></citet> of 44 tissue types, the set of Common Mind Consortium brain eQTLS<citet data-cite="fromer_gene_2016"><sup>fromer</sup></citet> and a set of brain-specific methylation QTL data sets currently available within the MRC Centre for Neuropsychiatric Genetics and Genomics were used. For each QTL a slope/$\beta$ limit of less than +/- 0.10 was applied, discarding any QTLs below this as having a negligible effect. SNP positions were converted to rs numbers using the Kaviar database<citet data-cite="glusman_kaviar:_2011"><sup>kaviar</sup></citet> and gene positions were annotated using the BioMart data mining tool<citet data-cite="smedley_biomart_2015"><sup>biomart</sup></citet>.
#
# The individual functional gene annotations were then combined by amalgamating the SNP list for the same gene appearing in more than one annotation and removing any duplicate SNPs from each gene's SNP list.
#
# In all 17,775 genes were functionally annotated. In comparison, the proximity annotation contains 19,071 genes.
#
# For each GWAS analysed, the combined functional gene annotation data set, GWAS summary statistics and the 1000 Genomes European genotype data set (as provided by the MAGMA project) are used as input to a MAGMA gene analysis, resulting in a list of genes with computed p-values assigned using the combined mean and top SNP-wise gene analysis model.

# ## Drug target / disease associations
#
# Pharmacological groups at each ATC level that are significantly associated with disease were generated using functional gene annotation datasets from several disease GWAS summary statistics. This provides up to 5 lists of pharmacologically related sets of drug targets for each disease GWAS - i.e. 5 results per ATC level.
#
# For easier reference, any pharmacological groups significantly associated with disease were collated across all ATC levels and placed in ascending q-value order. Over the course of the project so far, 16 different disease GWAS were analysed over a wide range of diseases, some with more well understood aetiology than others. Of the 16 GWAS analysed, those for diabetes, high blood pressure and tuberculosis (TB) are selected here, taken from the UK Biobank GWAS results<citet data-cite="noauthor_uk_nodate"><sup>ukbiobank</sup></citet>. These results are show in tables \ref{tab:Diabetesfunc} to \ref{tab:TBfunc}.

# +
# Get all ATC classifications into one table
atc_level_files = [
    "../../data/atc/level-1-atc.txt",
    "../../data/atc/level-2-atc.txt",
    "../../data/atc/level-3-atc.txt",
    "../../data/atc/level-4-atc.txt",
    "../../data/atc/level-5-atc.txt",
]

magma_results_dir = "../../data/magma_analysis/magma_results/"
summary_results_dir = "../../data/summaries/"
gwas_list = ["PD"]  # , 'HighBloodPressure', 'Diabetes']

# Set up some plotting variables
fig_unit = 17
fig_width = 2
fig_height = 4

# +
def get_atc_file(filename):
    file = open(filename, "r")
    lines = file.readlines()
    result = pd.DataFrame()
    for line in lines:
        code = line.split(" ")[0]
        desc = " ".join(line.split(" ")[1:]).strip()
        result = result.append(
            {"atc_code": code, "Description": desc}, ignore_index=True
        )

    file.close()
    return result


def get_atc_levels():
    # convert each atc_level_file to a 2-column data frame, and cat all of them together
    atc_levels = pd.concat(
        [get_atc_file(atc_level_file) for atc_level_file in atc_level_files]
    )
    atc_levels.set_index("atc_code", inplace=True)
    return atc_levels


def read_results_files(gwas, annot_type):
    # read in results files from a MAGMA gene set analysis
    prefix = magma_results_dir + gwas + "/results/"
    results_fileset = [
        prefix + "magma_geneset_result-" + gwas + "-" + annot_type + "-atc1.sets.out",
        prefix + "magma_geneset_result-" + gwas + "-" + annot_type + "-atc2.sets.out",
        prefix + "magma_geneset_result-" + gwas + "-" + annot_type + "-atc3.sets.out",
        prefix + "magma_geneset_result-" + gwas + "-" + annot_type + "-atc4.sets.out",
        prefix + "magma_geneset_result-" + gwas + "-" + annot_type + "-atc5.sets.out",
    ]

    results = [
        pd.read_csv(
            file, comment="#", delim_whitespace=True, usecols=["SET", "NGENES", "P"]
        )
        for file in results_fileset
    ]
    return results


# -


def summarise_drug_results(
    gwas, annot_type, atc_levels, n_gene_thresh=5, signif_thresh=0.05
):
    results = read_results_files(gwas, annot_type)

    # only consider classes/drugs with Qval < QVAL_THRESH and NGENES >= N_GENE_THRESH
    significants = [result[result["NGENES"] >= n_gene_thresh] for result in results]

    # Calculate the q-values for the local results per GWAS (per ATC code)
    for result in significants:
        result.loc[:, "q-value"] = qvalue.estimate(np.array(result["P"]))

    # Put all individual ATC results into one big dataframe
    all_significant = pd.concat(significants)

    q_local_significant = all_significant[all_significant["q-value"] < signif_thresh]
    q_local_final = pd.merge(
        q_local_significant, atc_levels, right_index=True, left_on="SET"
    ).sort_values("q-value")
    q_local_final.to_csv(
        summary_results_dir
        + "drugs_found-"
        + gwas
        + "-"
        + annot_type
        + "_local_qvals.tsv",
        sep="\t",
        index=False,
    )

    p_significant = all_significant[all_significant["P"] < signif_thresh]
    p_final = pd.merge(
        p_significant, atc_levels, right_index=True, left_on="SET"
    ).sort_values("P")
    p_final.to_csv(
        summary_results_dir + "drugs_found-" + gwas + "-" + annot_type + "_pvals.tsv",
        sep="\t",
        index=False,
    )


# +
def __main__():
    annot_type_list = ["func-brain-no_kaviar"]
    N_GENE_THRESH = 2
    SIGNIF_THRESH = 0.05
    ATC_LEVELS = get_atc_levels()

    for annot_type in annot_type_list:
        for gwas in gwas_list:
            summarise_drug_results(
                gwas, annot_type, ATC_LEVELS, N_GENE_THRESH, SIGNIF_THRESH
            )


if __name__ == "__main__":
    __main__()

# +
pval_files = [
    file
    for file in os.listdir(summary_results_dir)
    if file.endswith("_local_qvals.tsv")
]
pval_files.sort()

# boolean to put a newline after first table, for formatting
first_table = True

for file in pval_files:
    table = pd.read_csv(summary_results_dir + file, sep="\t", index_col=0)
    table.sort_values("q-value", inplace=True)
    table.columns = ["Gene count", "p-value", "q-value", "Description"]
    table.index.rename("ATC code", inplace=True)
    if len(table) != 0:
        gwas = file[12:-21]
        atype = file[-20:-16]
        display(Latex(table.to_latex(column_format="rrrrp{8cm}")))
        display(
            Latex(
                r"\caption{Table \ref{tab:"
                + gwas
                + atype
                + r"} - "
                + gwas
                + r" functional significant pharmacological group q-values.}"
            )
        )
        display(Latex(r"\label{tab:" + gwas + atype + r"}"))
        # display(Latex(r'\end{longtable}%'))
        if first_table:
            display(Latex(r"\newpage"))
            first_table = False
# -

# For each disease GWAS, correlations between the gene p-values resulting from the proximity and functional gene analyses were calculated, and correlation coefficients used to monitor the correctness of the functional annotation pipeline. Also for each disease GWAS, a histogram of p-values was plotted (Figure \ref{fig:functional_pvals}). Histograms without a well defined tail indicate a low powered GWAS and hence high FDR.

# ### Gene p-value histograms

# +
fig_height = 2
fig = plt.figure(figsize=(fig_unit, fig_unit * fig_height / fig_width))
fig.subplots_adjust(hspace=0.35, wspace=0.3)

for i, gwas in enumerate(gwas_list):
    prefix = magma_results_dir + gwas + "/results/"
    results_func = pd.read_csv(
        prefix + "magma_gene_result-" + gwas + "-func-brain-no_kaviar.genes.out",
        delim_whitespace=True,
        comment="#",
    )

    plt.subplot(fig_height, fig_width, i + 1)
    plt.hist(results_func["P_JOINT"], bins=500)
    plt.xlabel("p-value")
    plt.ylabel("Frequency")
    plt.title(gwas + " Functional genes")
# -

display(Latex(r"\begin{figure}[htpb]"))
display(Latex(r"\centering"))
display(
    Latex(
        r"\caption{Figure \ref{fig:functional_pvals} - Functional analysis gene p-value histograms.}"
    )
)
display(Latex(r"\label{fig:functional_pvals}"))
display(Latex(r"\end{figure}"))

# # Discussion
#
#
# Since defining this project, other research has been published with similar directions. For example, Gaspar and Breen<citet data-cite="gaspar_drug_2017"><sup>gaspar</sup></citet> find drug repurposing candidates in schizophrenia, mapping SNPs to genes by choosing the closest gene to a SNP. The FUMA (Functional mapping and annotation of genetic associations) tool<citet data-cite="watanabe_functional_2017"><sup>watanabe</sup></citet> uses functional data to map SNPs to genes, but does not involve any drug target database with which to perform gene-set analysis. This project differs from either approach by using both functional SNP-gene annotations and a drug target database in order to investigate novel drug repositioning opportunities.
#
# An encouraging set of significant pharmacological groups were produced for the tested disease GWAS.
# The gene p-value histograms shown in Figure \ref{fig:functional_pvals} shows one of the GWAS studies is underpowered (tuberculosis), which gives a large FDR and makes the pharmacological group results questionable. Parkinson's, schizophrenia and diabetes have much better gene p-value histograms with a well defined, low tail that result in more robust q-values.
#
# Good correlations (r > 0.5) between the functional and proximity analysis results for each disease GWAS, as well as several identical significant hits for some disease GWAS also indicates that the functional analysis pipeline created has no major errors or omissions.
#
# A cursory investigation of the ATC level 5 drugs associated with the disease GWAS provides corroborating evidence that they are not accidental, for example:
#
# * The diabetes results highlight the drug acitretin, which, according to the British Association of Dermatologists<citet data-cite="noauthor_acitretin_nodate"><sup>BAD</sup></citet>, "Acitretin can worsen diabetes...", and a study was also found investigating the effect of acitretin on glucose tolerance<citet data-cite="hartmann_effect_1992"><sup>hartmann</sup></citet>.
# * The high blood pressure results indicate a link to phenylpropanolamine, clenbuterol, labetalol and thiazides. Labetalol and thiazides are used to treat hypertension, and studies were found linking phenylpropanolamine and clenbuterol to changes in blood pressure<citet data-cite="salerno_impact_2005"><sup>salerno</sup></citet>.
# * Both quaternary ammonium compounds<citet data-cite="byrne_efficacy_1999"><sup>byrne</sup></citet> and Helicobacter pylori drugs<citet data-cite="perry_infection_2010"><sup>perry</sup></citet> were found to have evidence associating them with tuberculosis.
#
# These results also highlight the fact that this type of study does not differentiate between drugs that worsen, and drugs that aid a condition.

# # Appendix
#
# ## Thesis Outline
#
# * **Chapter 1: Introduction**
# * **Chapter 2: Curation of a functional gene annotation data set**
#     * 2.1 Introduction
#     * 2.2 Materials and methods
#     * 2.3 Results
#     * 2.4 Discussion
# * **Chapter 3: Generation of a drug-target / ontology database**
#     * 3.1 Introduction
#     * 3.2 Materials and methods
#     * 3.3 Results
#     * 3.4 Discussion
# * **Chapter 4: An analysis pipeline to identify pharmacologically related sets of drug targets associated with disease**
#     * 4.1 Introduction
#     * 4.2 Materials and methods
#     * 4.3 Implementation
#     * 4.4 Results
#     * 4.5 Discussion
# * **Chapter 5: A protein interactome network analysis of genes linked to drugs that are significantly associated with disease**
#     * 5.1 Introduction
#     * 5.2 Materials and methods
#     * 5.3 Implementation
#     * 5.4 Results
#     * 5.5 Discussion
# * **Chapter 6: Discussion & Conclusions**
#
#
# ## Research Timeline
# As is typical with modern software development, an agile, iterative approach will be used to develop the analysis pipeline further, as opposed to a traditional 'waterfall' approach. It is estimated that each iteration will take 6 months (indicating ~5 more iterations will take place over the lifetime of the project). Each iteration will seek to incorporate one or more of the following features; whilst culminating in a working deliverable that has been debugged, optimised and tested to a high degree of confidence:
#
# * Increasing the number of drugs analysed - identify more sources of drug - protein target data to incorporate
# * Increasing the number of genes considered - not all genes have a QTL, so alternative approaches will be considered in addition to incorporating QTLs of more tissue types. One alternative could be to use a hybrid approach of positional and functional SNP-gene associations
# * Incorporating exonic sequencing results, which are anticipated to become available in the near future
# * The effect of removing MHC and other large association areas on results
# * Using data within the pipeline to describe the magnitude and direction of expression changes associated with a drug
# * A protein interactome network analysis of genes linked to drugs that are sigificantly associated with disease
# * Increasing the number and types of visualisation
# * A graphical web front-end for wider use within the research community.
#

from IPython.display import display, HTML

HTML(
    """<script>
code_show=true; 
function code_toggle() {
   if (code_show){$('div.input').hide();} else {$('div.input').show();} code_show = !code_show} 
   $( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>."""
)
