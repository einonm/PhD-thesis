#! /usr/bin/env bash

if [[ $# -ne 5 ]]; then
    printf "Script not run. Please supply five parameters:\n  The MAGMA version\n  The disease GWAS\n  The annotation type\n  The geneset annotation\n  The run ID\n"
  exit 1
fi

MAGMA_VER=$1
GWAS=$2
ANNOT_TYPE=$3
GENESET=$4
RUNID=$5

if [ ${MAGMA_VER} == '1.06' ]; then
    module load raven
    module load magma
else
    module load magma/${MAGMA_VER}
fi

MAGMA="magma"

RESULTS_DIR="../data/magma_analysis/magma_results/${GWAS}/${RUNID}"

if [ ${GENESET} == 'gopaths' ]; then
  SETANNOT="../data/magma_analysis/GO_ALL_PC_20-2000_MAGMA.LH.ensembl.druggable.txt"
elif [ ${GENESET} == 'stitch' ]; then
  SETANNOT="../data/target-dbs/stitch_dgi_targets_atc_ensembl.csv"
elif [ ${GENESET} == 'duggie' ]; then
  SETANNOT="../data/target-dbs/all_dgi_targets_atc_ensembl.csv"
elif [ ${GENESET} == 'nostitch' ]; then
  SETANNOT="../data/target-dbs/nostitch_dgi_targets_atc_ensembl.csv"
fi

${MAGMA} \
  --gene-results ${RESULTS_DIR}/magma_gene_result-${GWAS}-${ANNOT_TYPE}-${RUNID}.genes.raw \
  --set-annot ${SETANNOT} \
  --out ${RESULTS_DIR}/magma_geneset_result-${GWAS}-${ANNOT_TYPE}-${RUNID}-${GENESET}

