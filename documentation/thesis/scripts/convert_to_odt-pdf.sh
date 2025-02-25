#! /usr/bin/env bash

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export BSTINPUTS=${BINDIR}

#strip file extension
FILE=`echo "$1" | cut -d"." -f1`

source /opt/anaconda3/etc/profile.d/conda.sh
JDIR=~/jupyter-notebooks
conda activate $JDIR

for EX in aux bbl blg log out tex toc;
do
    rm ${FILE}.$EX
done

jupyter nbconvert --to latex "${PWD}/${FILE}.ipynb" --template ${BINDIR}/latex_thesis.tplx #--debug

htlatex ${FILE}

# Run twice to get TOC etc.
bibtex  ${FILE}
htlatex ${FILE}
htlatex ${FILE}

soffice --headless \
        --convert-to odt \
        ${FILE}.html
