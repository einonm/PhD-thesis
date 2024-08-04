#! /usr/bin/env bash

source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate ~/jupyter-notebooks

# use all eQTL GTEX files
FILES=`find /neurocluster/databank/QTLs/GTEx/v7/ -name "*.signifpairs.txt"`
#FILES=$1
OUTDIR='../data/qtl/qtls-all'

if [ ! -e ${OUTDIR}/logs/ ]; then
    mkdir -p ${OUTDIR}/logs/
fi

for FILE in $FILES:
do
    # use the file part of the filename as the gene.annot file
    ANNOT=`echo $FILE | sed 's/.signifpairs.txt//' | sed 's/.*\///'`
    /opt/gridengine/bin/linux-x64/qsub \
      -q all.q \
      -cwd \
      -j y \
      -b y \
      -V \
      -l h_vmem=40G \
      -o ${OUTDIR}/logs/ \
      -e ${OUTDIR}/logs/ \
      -M einonm@cf.ac.uk \
      -m a \
      -N Parse_eQTL_data-nb \
      python ./Parse_eQTL_data-nb.py $FILE ${OUTDIR}/$ANNOT.gene.annot
done

