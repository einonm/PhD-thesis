# ---
# jupyter:
#   jupytext:
#     formats: py:nomarker
#     text_representation:
#       extension: .py
#       format_name: nomarker
#       format_version: '1.0'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import os
import sys
import pandas as pd
import numpy as np

import math

sys.path.append("../lib/")

# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm
import scipy.stats as stats

from IPython.display import display, HTML

import subprocess

import results_display as rdisp

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score

data_path = "../data/"
magma_results_dir = data_path + "magma_analysis/magma_results/"
summary_results_dir = data_path + "summaries/"
file_postfix = "duggie_qvals.tsv"
gwas = "PGC3_SZ"
run_id = "3qtl"

annots = [
    "prox",
    "alleqtlbrain",
    "alleqtlbrainhybrid",
    "alleqtlbrainhybridboth",
    "allbrain",
    "allbrainhybrid",
    "allbrainhybridboth",
    "alleqtlblood",
    "alleqtlbloodhybrid",
    "alleqtlbloodhybridboth",
    "alleqtltissues",
    "alleqtltissueshybrid",
    "alleqtltissueshybridboth",
    "alltissues",
    "alltissueshybrid",
    "alltissueshybridboth",
]

drug_sets = []

for annot in annots:
    pval_file = [
        file
        for file in os.listdir(summary_results_dir)
        if file.endswith(file_postfix)
        and "found-" + gwas + "-" + run_id + "-" + annot + "-" in file
    ]
    print(pval_file)
    table = pd.read_csv(os.path.join(summary_results_dir, pval_file[0]), sep="\t")
    table2 = table[table["Q " + annot + " duggie"] < 0.05]
    drug_set = set(table2["ATC_CODE"])
    drug_sets.append(drug_set)

len(drug_sets)

results_overlap = pd.DataFrame(columns=annots, index=annots, dtype=float)

for annot in annots:
    for annot2 in annots:
        results_overlap.at[annot, annot2] = float(
            rdisp.jaccard(drug_sets[annots.index(annot)], drug_sets[annots.index(annot2)])
        )

results_overlap = results_overlap.drop("prox").drop("alltissueshybridboth", axis=1)

results_overlap

mask_ut = np.triu(np.ones(results_overlap.shape)).astype(np.bool)
np.fill_diagonal(mask_ut, 0)

plt.figure(figsize=(16, 14))
sns.heatmap(results_overlap, mask=mask_ut, cmap="viridis", annot=True)
plt.show()


