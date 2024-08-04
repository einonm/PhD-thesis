#! /usr/bin/env bash

# Use the scripts directory to begin searching the directory tree
ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"/..

# Tell LaTeX to use this directory to look for style files
export BSTINPUTS=${ROOTDIR}/tplx/
export TEX=${ROOTDIR}/tex/

# Use the python environment set up for jupyter use
# From https://gitlab.com/einonm/jupyternb-setup
source /opt/anaconda3/etc/profile.d/conda.sh
JDIR=~/jupyter-notebooks
conda activate $JDIR

# The list of individual chapter and appendix notebooks
SECTIONS=("${TEX}/preliminary_pages/Declarations"
          "${TEX}/preliminary_pages/Acknowledgements"
          "${TEX}/preliminary_pages/Summary"
          "${TEX}/preliminary_pages/Abbreviations"
          "${TEX}/chapter01/Chapter01_Introduction"
          "${TEX}/chapter02/Chapter02_AnalysisPipelineToIdentifyRelatedSetsOfDrugTargetsAssociatedWithDisease"
          "${TEX}/chapter03/Chapter03_AppraisalOfThePipelineWithHypertensionAndHyperchloresterolemia"
          "${TEX}/chapter04/Chapter04_UsingFunctionalQTLGeneAnnotationToBetterIdentifyGenesAssociatedWithDisease"
          "${TEX}/chapter05/Chapter05_ApplicationOfTheAnalysisPipelineToNeuropsychiatricDiseases"
          "${TEX}/chapter06/Chapter06_DiscussionAndConclusions"
          "${TEX}/chapter07/Chapter07_SupplementaryInformation"
          "${TEX}/appendix_dgis/Appendix_DGIs"
          "${TEX}/appendix_magma_perm_issues/Appendix_magma_perm_issues"
          "${TEX}/appendix_disease_results/Appendix_disease_results"
          "${TEX}/appendix_pipeline_scripts/Appendix_pipeline_scripts")

for FILE in "${SECTIONS[@]}"; do
    cd $(dirname "${FILE}")
    CFILE=$(basename "${FILE}")

    # Convert chapter notebook to bare-bones latex to use as 'input' in Thesis.tex
    jupyter nbconvert --to latex "${CFILE}.ipynb" --template ${BSTINPUTS}/thesis_chapter.tplx
done

    cd ${TEX}/thesis/

    # Nasty hack to copy images to get thesis to generate
    cp -a ${TEX}/chapter02/Chapter02_AnalysisPipelineToIdentifyRelatedSetsOfDrugTargetsAssociatedWithDisease_files .
    cp -a ${TEX}/chapter03/Chapter03_AppraisalOfThePipelineWithHypertensionAndHyperchloresterolemia_files .
    cp -a ${TEX}/chapter04/Chapter04_UsingFunctionalQTLGeneAnnotationToBetterIdentifyGenesAssociatedWithDisease_files .
    cp -a ${TEX}/chapter05/Chapter05_ApplicationOfTheAnalysisPipelineToNeuropsychiatricDiseases_files .
    cp -a ${TEX}/appendix_magma_perm_issues/Appendix_magma_perm_issues_files .
    cp -a ${TEX}/appendix_disease_results/Appendix_disease_results_files .

    pdflatex thesis
    # Run twice to get TOC etc.
    bibtex  thesis
    pdflatex thesis
    pdflatex thesis

    # Remove any old temporary files
    for EX in aux bbl blg log out toc lof lot;
    do
        rm thesis.$EX
    done
