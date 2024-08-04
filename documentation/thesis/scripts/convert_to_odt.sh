#! /usr/bin/env bash

FILE=`echo "$1" | cut -d"." -f1`

source /opt/anaconda3/etc/profile.d/conda.sh
JDIR=~/jupyter-notebooks
conda activate $JDIR

jupyter nbconvert --to html "${FILE}.ipynb" --no-input

soffice --headless \
        --convert-to odt \
        ${FILE}.html
