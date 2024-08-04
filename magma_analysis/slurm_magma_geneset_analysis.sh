#! /bin/sh

MAGMA_VER=$1
GWAS=$2
ANNOT_TYPE=$3
RUNID=$4

# Removed option: 'nostitch'
for GENESET in 'duggie' 'stitch' 'gopaths'; do
	${PWD}/magma_geneset_analysis.sh ${MAGMA_VER} ${GWAS} ${ANNOT_TYPE} ${GENESET} ${RUNID}
done
