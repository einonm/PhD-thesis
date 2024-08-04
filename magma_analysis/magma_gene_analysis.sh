#! /usr/bin/env bash
#
# Script to aid the easy running of a magma gene analysis on SGE and Slurm.
# params:
# $1 - The MAGMA version string to use
# $2 - GWAS name (corresponding to data directory name)
# $3 - annotation type, proximity (prox), hybrid or functional (qtl), with optionally the QTL_TYPE
#      following: -brain or -all
# $4 - run number, for multiple runs
# $5 - The total number of batch jobs being used

MAGMA_VER=$1
GWAS=$2
TYPES=$3
RUNID=$4
BJOBS=$5

if [ ${MAGMA_VER} == '1.06' ]; then
    module load raven
    module load magma
else
    module load magma/${MAGMA_VER}
fi

MAGMA="magma"

BATCH_OPT="--batch $SLURM_ARRAY_TASK_ID $BJOBS"
UKBB_ROOT="/neurocluster/databank/ukbb_gwas_results"
BFILE="../data/magma_auxiliary_files/g1000_eur/g1000_eur"

USE="use=SNP,P"
RESULTS_DIR="../data/magma_analysis/magma_results/${GWAS}/${RUNID}"

if [ "${GWAS}" == "PD" ]; then
  PVAL="../data/gwas/PD_2014_META-asleon.txt"
  N="N=9000"
  USE="use=MarkerName,P.value"
elif [ "${GWAS}" == "SZ1" ]; then
  PVAL="../data/gwas/pgc.scz.full.2012-04.txt"
  N="N=9000"
  USE="use=snpid,pval"
elif [ "${GWAS}" == "SZ" ]; then
  PVAL="../data/gwas/daner_PGC_SCZ52_0513a.txt"
  N="N=9000"
elif [ "${GWAS}" == "PDpain" ]; then
  PVAL="../data/gwas/Log.McGill.assoc"
  N="N=9000"
elif [ "${GWAS}" == "PDcluster1" ]; then
  PVAL="../data/gwas/Cluster1.assoc"
  N="N=9000"
elif [ "${GWAS}" == "PDcluster2" ]; then
  PVAL="../data/gwas/Cluster2.assoc"
  N="N=9000"
elif [ "${GWAS}" == "PDcluster3" ]; then
  PVAL="../data/gwas/Cluster3.assoc"
  N="N=9000"
elif [ "${GWAS}" == "PDcluster4" ]; then
  PVAL="../data/gwas/Cluster4.assoc"
  N="N=9000"
elif [ "${GWAS}" == "Insomnia" ]; then
  PVAL="../data/gwas/gwas/GWAS_1.txt"
  N="N=9000"
elif [ "${GWAS}" == "Glaucoma" ]; then
  PVAL="../data/gwas/gwas/GWAS_2.txt"
  N="N=9000"
elif [ "${GWAS}" == "TB" ]; then
  PVAL="../data/gwas/gwas/GWAS_3.tsv"
  N="N=9000"
  USE="use=rsid,pval"
elif [ "${GWAS}" == "IPF" ]; then
  PVAL="../data/gwas/gwas/GWAS_4.tsv"
  N="N=9000"
  USE="use=rsid,pval"
elif [ "${GWAS}" == "HighBloodPressure" ]; then
  PVAL="../data/gwas/gwas/GWAS_5.tsv"
  N="N=9000"
  USE="use=rsid,pval"
elif [ "${GWAS}" == "Cancer" ]; then
  PVAL="../data/gwas/GWAS_6.tsv"
  N="N=9000"
  USE="use=rsid,pval"
elif [ "${GWAS}" == "Diabetes" ]; then
  PVAL="../data/gwas/GWAS_7.tsv"
  N="N=9000"
  USE="use=rsid,pval"
elif [ "${GWAS}" == "22q" ]; then
  PVAL="../data/gwas/GWAS_22Qonly.assoc.logistic"
  N="N=9000"
elif [ "${GWAS}" == "22qSigPC" ]; then
  PVAL="../data/gwas/GWAS_22Qonly_SigPC.assoc.logistic"
  N="N=9000"
elif [ "${GWAS}" == "HD" ]; then
  PVAL="../data/gwas/GWA12345_summary_rs_snps.txt"
  USE="use=snp_id,pval"
  N="N=9000"
  BFILE="../data/gwas/GWA4_12345-rs"
elif [ "${GWAS}" == "AD" ]; then
  PVAL="../data/gwas/IGAP_stage_1.txt"
  USE="use=MarkerName,Pvalue"
  N="N=9000"
elif [ "${GWAS}" == "PD2" ]; then
  PVAL="../data/gwas/PD-23am/METAANALYSIS_PdGeneAnd23AndMe.USE.assoc"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_E10" ]; then # Diabetes type I
  PVAL="../data/gwas/imputed.allWhites.clinical_c_E10-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_E11" ]; then
  PVAL="../data/gwas/imputed.allWhites.clinical_c_E11-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_I10" ]; then
  PVAL="../data/gwas/imputed.allWhites.clinical_c_I10-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_E78" ]; then
  PVAL="../data/gwas/imputed.allWhites.clinical_c_E78-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_N17" ]; then
  PVAL="../data/gwas/imputed.allWhites.clinical_c_Block_N17-N19-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_D65" ]; then
  PVAL="../data/gwas/imputed.allWhites.clinical_c_Block_D65-D69-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "PGC3_SZ" ]; then
  PVAL="../data/gwas/PGC3_SZ/daner_PGC_SCZ_w3_90_0418b-rs_pos_pvals"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_M05" ]; then
  PVAL="../data/gwas/imputed.allWhites.clinical_c_M05-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_K50" ]; then # Crohn's disease
  PVAL="../data/gwas/imputed.allWhites.clinical_c_K50-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_I25" ]; then # Chronic heart disease
  PVAL="../data/gwas/imputed.allWhites.clinical_c_I25-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_C50" ]; then # Breast cancer
  PVAL="../data/gwas/imputed.allWhites.clinical_c_C50-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "GA_G35" ]; then # Multiple sclerocsis
  PVAL="../data/gwas/imputed.allWhites.clinical_c_G35-filtered.csv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "CROHNS" ]; then # Crohn's disease, https://www.ibdgenetics.org/downloads.html#
  PVAL="../data/gwas/IBD/EUR.CD.gwas_info03_filtered.assoc"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "COVID" ]; then
  PVAL="../data/gwas/COVID19_HGI_ANA5_20200513_rs.txt"
  USE="use=RS,all_inv_var_het_p"
  N="N=9000"
elif [ "${GWAS}" == "AD2022" ]; then
  PVAL="../data/gwas/AD2022.tsv"
  USE="use=snp_id,PVALUE"
  N="N=9000"
elif [ "${GWAS}" == "MDD" ]; then
  PVAL="../data/gwas/MDD/MDD2018_ex23andMe_snp_p.tsv"
  USE="use=SNP,P"
  N="N=9000"
elif [ "${GWAS}" == "BP" ]; then
  PVAL="../data/gwas/BP/pgc-bip2021-all.vcf.tsv"
  USE="use=ID,PVAL"
  N="N=9000"
else
  GZASSOC="${UKBB_ROOT}/GWAS/`echo ${GWAS} | head -c 1`/${GWAS}.assoc.tsv.gz"
  PVAL="../data/gwas/${GWAS}.assoc.tsv"
  if [ ! -f ${PVAL} ]; then
    /bin/gunzip -c ${GZASSOC} > ${PVAL}
  fi
  echo ${PVAL}
  N="N=9000"
  USE="use=rsid,pval"
fi

ANNOT_TYPE=`echo ${TYPES} | sed 's/\(.*\)-.*/\1/'`
QTL_TYPE=`echo ${TYPES} | sed 's/.*-\(.*\)/\1/'`
if [ "${QTL_TYPE}" == "${TYPES}" ]; then
    # no specific QTL type
    QTL_TYPE=
fi

if [ "${ANNOT_TYPE}" == "prox" ]; then
  #GENE_ANNOT="../data/wrangling/magma_gene_atlas.prox.genes.annot.ensemble"
  GENE_ANNOT="../data/wrangling/magma_g1000eur_NCBI37.3.prox.genes.annot.ensemble"
elif [ "${ANNOT_TYPE}" == "qtl" ]; then
  GENE_ANNOT="../data/magma_analysis/combined.gene.annot-${QTL_TYPE}QTLs"
elif [ "${ANNOT_TYPE}" == "hybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/hybrid.gene.annot-${QTL_TYPE}QTLs"
elif [ "${ANNOT_TYPE}" == "hybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/hybridboth.gene.annot-${QTL_TYPE}QTLs"
elif [ "${ANNOT_TYPE}" == "newqtl" ]; then
  GENE_ANNOT="../data/magma_analysis/new_combined.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "newhybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/new_hybrid.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "newhybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/new_hybridboth.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "everyqtl" ]; then
  GENE_ANNOT="../data/magma_analysis/everything_combined.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "everyhybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/everything_hybrid.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "everyhybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/everything_hybridboth.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "eqtlgen" ]; then
  GENE_ANNOT="../data/magma_analysis/eqtlgen_final.tsv"
elif [ "${ANNOT_TYPE}" == "allbrain" ]; then
  GENE_ANNOT="../data/magma_analysis/all-brain.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtlbrain" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-brain.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtlblood" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-blood.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alltissues" ]; then
  GENE_ANNOT="../data/magma_analysis/all-tissues.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtltissues" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-tissues.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "allbrainhybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/all-brain_hybrid.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "allbrainhybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/all-brain_hybridboth.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtlbrainhybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-brain_hybrid.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtlbrainhybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-brain_hybridboth.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtlbloodhybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-blood_hybrid.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtlbloodhybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-blood_hybridboth.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alltissueshybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/all-tissues_hybrid.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alltissueshybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/all-tissues_hybridboth.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtltissueshybrid" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-tissues_hybrid.gene.annot-QTLs"
elif [ "${ANNOT_TYPE}" == "alleqtltissueshybridboth" ]; then
  GENE_ANNOT="../data/magma_analysis/all-eqtl-tissues_hybridboth.gene.annot-QTLs"
fi

#  --seed 1337 \
#  --gene-settings fixed-permp=5000 \
${MAGMA} \
  --bfile "${BFILE}" \
  --gene-annot "${GENE_ANNOT}" \
  --pval "${PVAL}" "${N}" "${USE}" \
  --gene-model multi=snp-wise \
  --gene-settings fixed-permp=20000 \
  ${BATCH_OPT} \
  --out ${RESULTS_DIR}/magma_gene_result-${GWAS}-${TYPES}-${RUNID}

