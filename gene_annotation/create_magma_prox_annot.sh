#! /bin/bash

module load raven
module load magma
MAGMA="magma"

# snp-loc file $1 is filtered to remove those with MAF<=0.01
$MAGMA \
    --snp-loc $1 \
    --gene-loc ../data/magma_auxiliary_files/g1000_eur/NCBI37.3.gene.loc \
    --annotate window=10,10 \
    --out $2
