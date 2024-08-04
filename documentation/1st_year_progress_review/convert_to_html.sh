#! /usr/bin/env bash

#strip file extension
FILE=`echo "$1" | cut -d"." -f1`

echo "$FILE"

jupyter nbconvert --to latex "${FILE}.ipynb" --template latex.tplx

sed -i 's/\\end{longtable}$//' "${FILE}.tex"

latex "${FILE}.tex"
bibtex  "${FILE}.aux"

# Run twice to get TOC etc.
htlatex "${FILE}.tex"
htlatex "${FILE}.tex"

htlatex "${FILE}.tex"

# View
google-chrome "${FILE}.html" &
