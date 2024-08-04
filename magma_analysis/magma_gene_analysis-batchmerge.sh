#! /usr/bin/env bash
#
# Script to merge a batch run of magma.
# params:
# $1 - The MAGMA version string
# $2 - GWAS name (corresponding to data directory name)
# $3 - annotation type, proximity (prox), functional (qtl) or hybrid
# $4 - RUNID - ID unique to this run

MAGMA_VER=$1
GWAS=$2
TYPES=$3
RUNID=$4

if [ ${MAGMA_VER} == '1.06' ]; then
    module load raven
    module load magma
else
    module load magma/${MAGMA_VER}
fi

MAGMA="magma"
BATCH_PREFIX="../data/magma_analysis/magma_results/${GWAS}/${RUNID}/magma_gene_result-${GWAS}-${TYPES}-${RUNID}"

${MAGMA} \
  --merge "${BATCH_PREFIX}" \
  --out "${BATCH_PREFIX}"
