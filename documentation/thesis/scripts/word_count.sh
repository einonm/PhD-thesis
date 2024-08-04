#! /bin/sh
ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"/..
export TEX=${ROOTDIR}/tex/
echo "texdir=${TEX}"

# Accurate word count (excluding summary, acknowledgements, contents pages, appendices, tables, diagrams and figures, references, bibliography, footnotes and endnotes).
cd ${TEX}; echo chapter01/Chapter01_Introduction.tex \
    chapter02/Chapter02_AnalysisPipelineToIdentifyRelatedSetsOfDrugTargetsAssociatedWithDisease.tex \
    chapter03/Chapter03_AppraisalOfThePipelineWithHypertensionAndHyperchloresterolemia.tex \
    chapter04/Chapter04_UsingFunctionalQTLGeneAnnotationToBetterIdentifyGenesAssociatedWithDisease.tex \
    chapter05/Chapter05_ApplicationOfTheAnalysisPipelineToNeuropsychiatricDiseases.tex \
    chapter06/Chapter06_DiscussionAndConclusions.tex \
    chapter07/Chapter07_SupplementaryInformation.tex \
    | xargs texcount
    #appendix_dgis/Appendix_DGIs.tex \
#find ${TEX} -name "*.tex" | xargs texcount

# Word count
WRDCNT=`pdftotext ${TEX}/thesis/thesis.pdf - | tr -d '.' | wc -w`
echo "Thesis PDF has $WRDCNT words (Gross, includes tables, appendicies and preliminary pages)"
