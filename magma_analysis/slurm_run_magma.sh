#! /usr/bin/env bash

# The MAGMA version, 1.06/1.08/1.08b/1.09b
MAGMA_VER='1.09b'

# Perform Magma gene analyses
ANALYSE_GENES=true
# perform Magma geneset analyses
ANALYSE_GENESETS=true

# SLURM options
QSUB="sbatch --parsable --account=scw1349 --time=2-00:00:00" #--qos=maxjobs500"
Q="c_compute_neuro1_long"
BJOBS=22 # Number of batch jobs to split gene analysis

RUNID="$1"
GWAS_LIST=${@:2}

#ANNOT_TYPE_LIST=('prox' 'qtl-brain' 'qtl-all' 'hybrid-brain' 'hybrid-all' 'hybridboth-brain' 'hybridboth-all')
#ANNOT_TYPE_LIST=('prox' 'qtl-all' 'hybrid-all' 'hybridboth-all' 'newqtl' 'newhybrid' 'newhybridboth' 'everyqtl' 'everyhybrid' 'everyhybridboth')
ANNOT_TYPE_LIST=('prox' 'allbrain' 'alleqtlbrain' 'alleqtlblood' 'alleqtltissues' 'alltissues' 'allbrainhybrid' 'alleqtlbrainhybrid' 'alleqtlbloodhybrid' 'alleqtltissueshybrid' 'alltissueshybrid' 'allbrainhybridboth' 'alleqtlbrainhybridboth' 'alleqtlbloodhybridboth' 'alleqtltissueshybridboth' 'alltissueshybridboth')

# alternate hybrid containing all SNPs from positional + QTL
#ANNOT_TYPE_LIST=('hybridboth-brain' 'hybridboth-all')

for GWAS in ${GWAS_LIST[@]}; do
  echo Submitting $RUNID for GWAS $GWAS
  RESULTS_DIR="../../../data/gnt-data/magma_analysis/magma_results/${GWAS}/${RUNID}"
  if [ ! -d "${RESULTS_DIR}" ]; then
    mkdir -p "${RESULTS_DIR}"
  fi

  GWAS_LOGDIR="${RESULTS_DIR}/logs"
  if [ ! -d "${RESULTS_DIR}/logs" ]; then
    mkdir -p "${RESULTS_DIR}/logs"
  fi

 JID=2
 JID2=2
  for ANNOT_TYPE in ${ANNOT_TYPE_LIST[@]}; do
    if ${ANALYSE_GENES}; then
        CMMD="$QSUB \
         -o ${GWAS_LOGDIR}/gene-%A_%a.out \
  	     -p ${Q} \
         --mem=12G \
         -a 1-${BJOBS} \
             $PWD/magma_gene_analysis.sh ${MAGMA_VER} ${GWAS} ${ANNOT_TYPE} ${RUNID} ${BJOBS}"
        echo $CMMD
        JID=$($CMMD)
      CMMD="$QSUB \
         -o ${GWAS_LOGDIR}/genemerge-%A_%a.out \
	     -p ${Q} \
         --mem=6G \
         --dependency=afterok:${JID} \
             $PWD/magma_gene_analysis-batchmerge.sh ${MAGMA_VER} ${GWAS} ${ANNOT_TYPE} ${RUNID}"
        echo $CMMD
        JID2=$($CMMD)
    fi
    if ${ANALYSE_GENESETS}; then
        CMMD="$QSUB \
            -o ${GWAS_LOGDIR}/geneset-%A_%a.out \
	        -p ${Q} \
            --mem=6G \
            --mail-user=einonm@cf.ac.uk \
            --mail-type=END \
            --dependency=afterok:${JID2} \
               ${PWD}/slurm_magma_geneset_analysis.sh ${MAGMA_VER} ${GWAS} ${ANNOT_TYPE} ${RUNID}"
        echo $CMMD
        $CMMD
    fi
  done
done
