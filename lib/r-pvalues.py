# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# +
#!/usr/bin/env Rscript

library(qvalue)
# -

# We want to use the Storey R library form calculating p-values, and particularly to obtain pi0 for a set of p values, the estimated FDR.

# +
args = commandArgs(trailingOnly=TRUE)

filename = args[1]

data <- read.table('../data/magma_analysis/magma_results/GA_E78/results/magma_gene_result-GA_E78-hybrid-all-final_drugs.genes.out',
                   header = TRUE)
# -

data

result = qvalue(data$P_JOINT)

result


