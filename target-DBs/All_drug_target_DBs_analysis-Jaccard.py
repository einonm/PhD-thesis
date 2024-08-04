# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import os
import numpy as np
import pandas as pd
from IPython.display import Image, display, HTML, Latex, display_latex, display_html

# #!pip install venn
# %matplotlib inline
import matplotlib.pyplot as plt
import venn
import csv

# #!pip install seaborn
import seaborn as sns

# #!pip install pandas_profiling
#import pandas_profiling

# +
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects import numpy2ri
from rpy2.robjects.lib import grdevices
from rpy2.robjects.vectors import IntVector

from rpy2.robjects.conversion import localconverter

rggp = importr("ggplot2")
graphics = importr("graphics")


# +
# Temp -recipe for showing R graphics inline
# ith grdevices.render_to_bytesio(grdevices.jpeg, width=1024, height=896, res=150) as img:
#   graphics.barplot(IntVector((1,3,2,5,4)), ylab="Value")

# isplay(Image(data=img.getvalue(), format='jpeg', embed=True))
# -

# create an all_targets_entrez/len_all_entrez columns of de-duped targets for an ATC drug table
def collate_all_targets(in_df):
    df = in_df.copy()
    # Convert each string of entrez IDs and ATC codes to a set list
    for name in ["targets_entrez_" + dgi for dgi in DGIs]:
        df[name] = [x.split() for x in df[name]]
        df[name] = [set(x) for x in df[name]]
        [x.remove("nan") if "nan" in x else x for x in df[name]]
        df[name] = [list(x) for x in df[name]]

    df["all_targets_entrez"] = (
        df["targets_entrez_TTD"]
        + df["targets_entrez_DrugBank"]
        + df["targets_entrez_DGIdb"]
        + df["targets_entrez_STITCH"]
        + df["targets_entrez_GtoPdb"]
        + df["targets_entrez_T3DB"]
        + df["targets_entrez_DrugCentral"]
        + df["targets_entrez_DsigDB"]
    )

    df["all_targets_entrez"] = [set(x) for x in df["all_targets_entrez"]]
    df["len_all_entrez"] = [len(x) for x in df["all_targets_entrez"]]

    return df


# Format an ATC table ready to convert to CSV useful to show graphs in thesis
def archive_atc_table(in_df):
    df = in_df[
        [
            "atc_code",
            "targets_entrez_TTD",
            "targets_entrez_DrugBank",
            "targets_entrez_DGIdb",
            "targets_entrez_STITCH",
            "targets_entrez_GtoPdb",
            "targets_entrez_T3DB",
            "targets_entrez_DsigDB",
            "targets_entrez_DrugCentral",
            "all_targets_entrez",
            "len_all_entrez",
        ]
    ]

    for dgi in DGIs:
        for index, row in df.iterrows():
            df.at[index, "targets_entrez_" + dgi] = " ".join(
                [str(x) for x in df.at[index, "targets_entrez_" + dgi]]
            )
            if not df.at[index, "targets_entrez_" + dgi]:
                df.at[index, "targets_entrez_" + dgi] = "nan"

    for index, row in df.iterrows():
        df.at[index, "all_targets_entrez"] = " ".join(
            [str(x) for x in df.at[index, "all_targets_entrez"]]
        )

    return df


# define Jaccard similarity function
def jaccard(set1, set2):
    intersection = len(list((set1).intersection(set2)))
    union = (len(list(set1)) + len(list(set2))) - intersection

    return float(intersection) / union


# ## Drug target DB analysis
#
# Here, we want to explore the database that we've created - looking at the contributing data sources and measuring the improvement each one provides - e.g. if no new targets or drugs are uniquely contributed by a database, it could be removed from the analysis.
#
# TODO - when code is hidden, output and notes should make sense.

HTML(
    """<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>."""
)

# Look at the collected data from all contributing databases:

jointall_df = pd.read_csv(
    "../data/target-dbs/jointall_targets-8dbs.tsv",
    sep="\t",
    dtype={
        "PubChem CID": str,
        "ATC_TTD": str,
        "targets_entrez_TTD": str,
        "ATC_DrugBank": str,
        "targets_entrez_DrugBank": str,
        "targets_entrez_DGIdb": str,
        "targets_entrez_STITCH": str,
        "targets_entrez_GtoPdb": str,
        "targets_entrez_T3DB": str,
        "ATC_DsigDB": str,
        "targets_entrez_DsigDB": str,
        "targets_entrez_DrugCentral": str,
        "ATC_DrugCentral": str,
    },
)
jointall_df["targets_entrez_DsigDB"] = jointall_df["targets_entrez_DsigDB"].astype(str)
jointall_df["PubChem CID"] = jointall_df["PubChem CID"].astype(str)
jointall_df.fillna(value="nan", inplace=True)
jointall_df

# ## Process all contributing target lists and amalgamate into one list for each drug

DGIs = ["STITCH", "DGIdb", "DrugCentral", "DrugBank", "DsigDB", "GtoPdb", "TTD", "T3DB"]

# Collect all ATC codes together to make one column
jointall_df["all_ATC"] = (
    jointall_df["ATC_TTD"].str.split()
    + jointall_df["ATC_DrugBank"].str.split()
    + jointall_df["ATC_DsigDB"].str.split()
    + jointall_df["ATC_DrugCentral"].str.split()
)
jointall_df.drop(
    ["ATC_TTD", "ATC_DrugBank", "ATC_DsigDB", "ATC_DrugCentral"], axis=1, inplace=True
)

# +
# Concatenate all ATC values per PubChemID - note may have duplicates for different pubchem IDs!
jointall_df["all_ATC"] = [set(x) for x in jointall_df["all_ATC"]]

[x.remove("nan") if "nan" in x else x for x in jointall_df["all_ATC"]]
jointall_df["all_ATC"] = [sorted(list(x)) for x in jointall_df["all_ATC"]]
jointall_df = jointall_df[jointall_df["all_ATC"].astype(bool)]
# -

# ### Examine ATC codes in detail to increase target lists
#
# The relationships between ATC codes were explored and text processing used to generate three lists covering all assigned ATC codes:
#
#  1. Individual, single drug ATC codes. e.g. prednisolone - A07EA01
#  2. Mixture drug ATC codes, where two or more individual drugs are used together. The ATC codes of the individual drugs used were noted against each mixture ATC code, e.g. prednisolone and promethazine - V03AB05; A07EA01 and D04AA10
#  3. Combination drugs, which are either derivatives of individual drugs or mixture drugs with a broad category and cannot be assigned to individual ATC codes (e.g. psycholeptics). It is assumed this category would either have a target list duplicated with an individual drug or not have an assignable target list unless one is assigned explicitly by a source interaction database. e.g. prednisolone, combinations - A01AC54
#
# We want a unique list of ATC codes with the targets for each - we can collate the member ATC target lists of mixture drugs to increase the size of the mixture drugs' target list, and also remove many combination ATC codes if their target lists match that of the individual drug it derives from.

# First get a table listed by unique ATC code.
atc_exploded_df = jointall_df.explode("all_ATC")

# +
# Gather all records for one ATC code into one row
# This ensures rows with matching ATC codes also have matching/joined target lists
atc_exploded_df.set_index("all_ATC", inplace=True)
atc_exploded_df = atc_exploded_df.groupby("all_ATC").transform(lambda x: " ".join(x))
atc_exploded_df = atc_exploded_df.reset_index()

# duplicate rows with matching ATC/target lists occur after the explode/groupby, drop them
atc_exploded_df = atc_exploded_df.drop_duplicates()

atc_exploded_df
# -

# ### Examine all ATC codes that are individual drugs

single_atc_df = pd.read_csv("../data/atc/single_drugs_atc5.tsv", sep="\t")
single_atc_targets_df = pd.merge(
    single_atc_df, atc_exploded_df, left_on="atc_code", right_on="all_ATC", how="inner"
)
single_atc_targets_df

# ### Examine all ATC codes that are mixture drugs

# TODO - combine the member atc code targets for these
mixture_atc_df = pd.read_csv("../data/atc/mixture_drugs_atc5.tsv", sep="\t")
mixture_atc_targets_df = pd.merge(
    mixture_atc_df, atc_exploded_df, left_on="atc_code", right_on="all_ATC", how="inner"
)
mixture_atc_targets_df

mixture_atc_targets_df["member_atc_codes"] = mixture_atc_targets_df[
    "member_atc_codes"
].str.split()

# Add targets for mixture members, if available
for index, row in mixture_atc_targets_df.iterrows():
    for atc in row["member_atc_codes"]:
        if atc in list(single_atc_targets_df["atc_code"]):
            # single_atc is the row of targets per DB for this single drug ATC code
            single_atc = single_atc_targets_df[single_atc_targets_df["atc_code"] == atc]
            # Add all single_atc code targets to this mixture entry
            for db_name in ["targets_entrez_" + dgi for dgi in DGIs]:
                mixture_atc_targets_df.at[index, db_name] = (
                    mixture_atc_targets_df.at[index, db_name]
                    + " "
                    + single_atc[db_name].item()
                )

# ### Examine all ATC codes that are 'combinations of' drugs

combination_atc_df = pd.read_csv("../data/atc/combination_drugs_atc5.tsv", sep="\t")
combination_atc_targets_df = pd.merge(
    combination_atc_df,
    atc_exploded_df,
    left_on="atc_code",
    right_on="all_ATC",
    how="inner",
)
combination_atc_targets_df

# so 389+176+2499 = 3064, but there are 3086 drugs initially. What are the other 22?
combined_set = set(combination_atc_targets_df["atc_code"]).union(
    set(mixture_atc_targets_df["atc_code"])
)
combined_set = combined_set.union(set(single_atc_targets_df["atc_code"]))
display(len(combined_set))

# ## Measure contribution of  ATC drugs from each DB
#
# ### Contributions of 1 or more targets to drugs from each DB

# +
intermediate_single_df = collate_all_targets(single_atc_targets_df)
intermediate_mix_df = collate_all_targets(mixture_atc_targets_df)

intermediate_combi_df = collate_all_targets(combination_atc_targets_df)
intermediate_atc_df = pd.concat(
    [intermediate_combi_df, intermediate_single_df, intermediate_mix_df],
    ignore_index=True,
)

for index, row in intermediate_atc_df.iterrows():
    intermediate_atc_df.at[index, "entrez"] = " ".join(
        [str(x) for x in intermediate_atc_df.at[index, "all_targets_entrez"]]
    )

# Remove drugs with no targets, also do so for split tables (combi/mix/single) to check numbers
intermediate_combi_df = intermediate_combi_df[
    intermediate_combi_df["all_targets_entrez"] != set()
]
intermediate_single_df = intermediate_single_df[
    intermediate_single_df["all_targets_entrez"] != set()
]
intermediate_mix_df = intermediate_mix_df[
    intermediate_mix_df["all_targets_entrez"] != set()
]
intermediate_atc_df = intermediate_atc_df[
    intermediate_atc_df["all_targets_entrez"] != set()
]

display(
    "Combination ATC drugs with >1 target number: " + str(len(intermediate_combi_df))
)
display(
    "Individual ATC drugs with >1 target number: " + str(len(intermediate_single_df))
)
display("mixture ATC drugs with >1 target number: " + str(len(intermediate_mix_df)))

display(
    "There are " + str(len(intermediate_atc_df)) + " ATC codes with at least 1 target"
)
# -

archive_atc_table(intermediate_atc_df).to_csv(
    "../data/target-dbs/intermediate_atc_targets.tsv", sep="\t", index=None
)
# df.to_csv('../data/target-dbs/intermediate_atc_targets.tsv', sep='\t', index=None)

# +
atcs = [None] * (len(DGIs) + 1)
set_names = []

for i in range(len(DGIs)):
    col_name = "targets_entrez_" + DGIs[i]
    atcs[i] = set(
        [
            x
            for x in intermediate_atc_df[col_name].index
            if len(intermediate_atc_df.at[x, col_name]) >= 1
        ]
    )
    if "nan" in atcs[i]:
        atcs[i].remove("nan")
    set_names.append(DGIs[i] + ": " + str(len(atcs[i])))
    display(
        "From "
        + str(len(intermediate_atc_df))
        + " QUAGMIRE ATC drugs >= 1 target, "
        + str(len(atcs[i]))
        + " have targets from "
        + DGIs[i]
    )

# -

# How much of this gene set is each DB uniquely capturing?
for name in ["targets_entrez_" + dgi for dgi in DGIs]:
    dgi_genes = [intermediate_atc_df.at[x, name] for x in intermediate_atc_df.index]
    all_list = []
    for item in dgi_genes:
        for elem in item:
            all_list.append(elem)
    display(name + ": " + str(len(set(all_list))))

dgi_genes = [
    intermediate_atc_df.at[x, "all_targets_entrez"] for x in intermediate_atc_df.index
]
all_list = []
for item in dgi_genes:
    for elem in item:
        all_list.append(elem)
display("all_targets_entrez: " + str(len(set(all_list))))

# ### Contributions of 5 or more targets to drugs from each DB

# +
intermediate_combi_df = intermediate_combi_df[
    intermediate_combi_df["len_all_entrez"] >= 5
]
intermediate_single_df = intermediate_single_df[
    intermediate_single_df["len_all_entrez"] >= 5
]
intermediate_mix_df = intermediate_mix_df[intermediate_mix_df["len_all_entrez"] >= 5]
intermediate_atc_df = intermediate_atc_df[intermediate_atc_df["len_all_entrez"] >= 5]

display(
    "Combination ATC drugs with >5 target number: " + str(len(intermediate_combi_df))
)
display(
    "Individual ATC drugs with >5 target number: " + str(len(intermediate_single_df))
)
display("mixture ATC drugs with >5 target number: " + str(len(intermediate_mix_df)))

display(
    "There are "
    + str(len(intermediate_atc_df))
    + " ATC coded drugs with 5 or more targets"
)
display(
    "(There are "
    + str(len(intermediate_atc_df.drop_duplicates(subset=["entrez"])))
    + " ATC coded drugs with 5 or more targets excluding duplicates)"
)

# +
atcs = [None] * len(DGIs)
set_names = []

for i in range(len(DGIs)):
    col_name = "targets_entrez_" + DGIs[i]
    atcs[i] = set(
        [
            x
            for x in intermediate_atc_df[col_name].index
            if len(intermediate_atc_df.at[x, col_name]) >= 5
        ]
    )
    if "nan" in atcs[i]:
        atcs[i].remove("nan")
    set_names.append(DGIs[i] + ": " + str(len(atcs[i])))
    display(
        "From "
        + str(len(intermediate_atc_df))
        + " QUAGMIRE ATC drugs >= 5 targets, "
        + str(len(atcs[i]))
        + " have targets from "
        + DGIs[i]
    )

# +
# Create a table of Jaccard co-effs giving overlap targets from all databases and DUGGIE -
# for >5 target drugs

#

# -

# ### Remove combination ATC drug duplicates of individual drugs

# +
# Annotate combi drugs with individual ATC code target lists
final_combi_df = collate_all_targets(combination_atc_targets_df)
final_single_df = collate_all_targets(single_atc_targets_df)
final_mix_df = collate_all_targets(mixture_atc_targets_df)
final_combi_df["single_targets"] = ""

for index, row in final_combi_df.iterrows():
    result = final_single_df[final_single_df["atc_code"] == row["single ATC"]][
        "all_targets_entrez"
    ].values
    if len(result):
        final_combi_df.at[index, "single_targets"] = result[0]
    else:
        final_combi_df.at[index, "single_targets"] = set()

# set operations - difference between combination and single ATC code target lists for each combination drug
final_combi_df["difference"] = (
    final_combi_df["all_targets_entrez"] - final_combi_df["single_targets"]
)

# Number of combination drugs that don't match the corresponding individual ATC's targets list
# But we need to remove drugs with < 5 targets to compare properly here
final_combi_df2 = final_combi_df[final_combi_df["len_all_entrez"] >= 5]
display(
    "We can remove "
    + str(len(final_combi_df2[final_combi_df2["difference"] == set()]))
    + " duplicated combi drugs"
)
display(
    "leaving "
    + str(len(final_combi_df2[final_combi_df2["difference"] != set()]))
    + " combi drugs"
)

# Let's drop matching combination drugs, as they are duplicates (and we don't want them replacing the individual ATC codes when we
# drop all duplicates later on)
final_combi_df = final_combi_df.drop(
    final_combi_df[final_combi_df["difference"] == set()].index
)
final_combi_df = final_combi_df.drop(
    ["single_targets", "difference", "single ATC", "single description"], axis=1
)

all_atc_df = pd.concat(
    [final_combi_df, final_single_df, final_mix_df], ignore_index=True
)

for index, row in all_atc_df.iterrows():
    all_atc_df.at[index, "entrez"] = " ".join(
        [str(x) for x in all_atc_df.at[index, "all_targets_entrez"]]
    )

# Only look at drug with >= 5 targets
all_atc_df = all_atc_df[all_atc_df["len_all_entrez"] >= 5]

all_atc_df = all_atc_df.drop(["all_ATC", "member_atc_codes", "PubChem CID"], axis=1)
all_atc_df = all_atc_df.reset_index()
display(
    "There are "
    + str(len(all_atc_df))
    + " ATC coded drugs with 5 or more targets after removing duplicate combi drugs"
)
# -

all_atc_df

# ### **Convert the entrez IDs to Ensembl, to make comparison with GTEx QTL data later on**

conversion_df = pd.read_csv(
    "/home/mpnme/data/phd/wrangling/Ensembl_prot_entrez.tsv", sep="\t"
)
conversion_df.drop(["Ensembl_prot", "HGNC symbol"], axis=1, inplace=True)
conversion_df.dropna(subset=["Entrez ID"], inplace=True)
conversion_df["Entrez ID"] = conversion_df["Entrez ID"].astype(int)
conversion_df.drop_duplicates(subset="Entrez ID", inplace=True)
conversion_df.set_index("Entrez ID", inplace=True)
conversion_df.sort_index(inplace=True)


all_atc_df["all_ensembl"] = np.empty((len(all_atc_df), 0)).tolist()
for dgi in DGIs:
    all_atc_df[dgi + "_ensembl"] = np.empty((len(all_atc_df), 0)).tolist()
fail_list = []
for index, row in all_atc_df.iterrows():
    ensembl = []
    for entrez in row["all_targets_entrez"]:
        int_entrez = int(entrez)
        if int_entrez in conversion_df.index:
            ensembl.append(conversion_df.at[int_entrez, "Ensembl"])
        else:
            fail_list.append(entrez)
    all_atc_df.at[index, "all_ensembl"] = ensembl

    for dgi in DGIs:
        ensembl = []
        for entrez in row["targets_entrez_" + dgi]:
            int_entrez = int(entrez)
            if int_entrez in conversion_df.index:
                ensembl.append(conversion_df.at[int_entrez, "Ensembl"])
        all_atc_df.at[index, dgi + "_ensembl"] = ensembl


for index, row in all_atc_df.iterrows():
    all_atc_df.at[index, "gene_list"] = " ".join(
        [str(x) for x in sorted(list(set(all_atc_df.at[index, "all_ensembl"])))]
    )

all_atc_df["len_all_ensembl"] = [len(x) for x in all_atc_df["all_ensembl"]]
all_atc_df = all_atc_df[all_atc_df["len_all_ensembl"] >= 5]
final_df = all_atc_df.copy()  # [['atc_code', 'gene_list']]
final_df.drop(
    ["description", "index", "len_all_ensembl", "len_all_entrez", "entrez"],
    axis=1,
    inplace=True,
)
final_df.index.sort_values()
display(final_df.head())

## Test attempt to find the evidence source for a particular drug-gene interaction (HTT, entrz 3064 / ENSG00000197386)
htt_drugs = final_df[final_df['gene_list'].str.contains('ENSG00000130164')]
display(htt_drugs)
for index,row in htt_drugs.iterrows():
    for column in htt_drugs.columns:
        if column.__contains__('ensembl'):
            if row[column].__contains__('ENSG00000130164'):
                print(column)

## Test attempt to find the evidence source for a particular drug-gene interaction
for column in htt_drugs.columns:
    print(column)
    display(htt_drugs.at[422, column])

# ### Calculate jaccard index map of overlap between contributing DB genesets and DUGGIE genesets for each drug

# +
for dgi in DGIs:
    final_df.drop(["targets_entrez_" + dgi], axis=1, inplace=True)

final_df.drop(["all_targets_entrez"], axis=1, inplace=True)
# -

# Gene lists should be sorted for this to work
final_df.drop_duplicates(subset=["gene_list"], inplace=True)
display(final_df)

for dgi in DGIs:
    final_df[dgi] = [
        jaccard(
            set(final_df.at[x, "all_ensembl"]), set(final_df.at[x, dgi + "_ensembl"])
        )
        for x in final_df.index
    ]

final_df.head()

heatmap_df = final_df[
    [
        "STITCH",
        "DGIdb",
        "DrugCentral",
        "DrugBank",
        "DsigDB",
        "GtoPdb",
        "TTD",
        "T3DB",
        "atc_code",
    ]
]
heatmap_df.set_index(["atc_code"], inplace=True)
heatmap_df.sort_values(["STITCH"], inplace=True, ascending=False)
heatmap_df

plt.figure(figsize=(20, 20))
sns.heatmap(heatmap_df, cmap="YlGnBu")
plt.show()

# ### Generate sets of drugs for each database, of DUGGIE drugs to which they contributed

# +
# Non-contributory DBs have a jaccard index of 0 where they don't contribute targets to a drug
db_drugsets = {}
drugset_sizes = {}
db_genesets = {}
geneset_sizes = {}

for dgi in DGIs:
    db_drugsets[dgi] = set(final_df[final_df[dgi] != 0]["atc_code"])
    drugset_sizes[dgi] = len(db_drugsets[dgi])
    db_genesets[dgi] = set(
        final_df[dgi + "_ensembl"].apply(pd.Series).stack().reset_index(drop=True)
    )
    geneset_sizes[dgi] = len(db_genesets[dgi])
# -

# Look at the overlap of each DB with other DBs
db_drug_overlap = pd.DataFrame(columns=DGIs, index=DGIs, dtype=float)

for dgi in DGIs:
    for dgi2 in DGIs:
        db_drug_overlap.at[dgi, dgi2] = float(
            jaccard(db_drugsets[dgi], db_drugsets[dgi2])
        )

db_drug_overlap = db_drug_overlap.drop("STITCH").drop("T3DB", axis=1)
db_drug_overlap.to_csv("../data/target-dbs/all_dgi_drug_jaccard.csv", header=True)
db_drug_overlap

# +
mask_ut = np.triu(np.ones(db_drug_overlap.shape)).astype(np.bool)
np.fill_diagonal(mask_ut, 0)

plt.figure(figsize=(12, 10))
sns.heatmap(db_drug_overlap, mask=mask_ut, cmap="viridis", annot=True)
# plt.savefig('foo.png', bbox_inches='tight')
plt.show()


# +
# from PIL import Image
# im = Image.open('foo.png')
# angle = 45
# out = im.rotate(angle, expand=True)
# out.save('rotate-output.png')

# +
plt.figure(figsize=(10, 10))
plt.bar(list(drugset_sizes.keys()), drugset_sizes.values())

with open("../data/target-dbs/all_dgi_drugset_sizes.csv", "w") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=drugset_sizes.keys())
    writer.writeheader()
    writer.writerow(drugset_sizes)

plt.show()

# +
# Look at the overlap of each DB with other DBs
db_gene_overlap = pd.DataFrame(columns=DGIs, index=DGIs, dtype=float)

# Look at the correlation between contributed gene targets
for dgi in DGIs:
    for dgi2 in DGIs:
        db_gene_overlap.at[dgi, dgi2] = float(
            jaccard(db_genesets[dgi], db_genesets[dgi2])
        )
# -

db_gene_overlap

db_gene_overlap = db_gene_overlap.drop("STITCH").drop("T3DB", axis=1)
db_gene_overlap.to_csv("../data/target-dbs/all_dgi_targets_jaccard.csv", header=True)

# +
mask_ut = np.triu(np.ones(db_gene_overlap.shape)).astype(np.bool)
np.fill_diagonal(mask_ut, 0)

plt.figure(figsize=(12, 10))
sns.heatmap(db_gene_overlap, mask=mask_ut, cmap="viridis", annot=True)
# plt.savefig('foo.png', bbox_inches='tight')
plt.show()
# -

plt.figure(figsize=(10, 10))
plt.bar(list(geneset_sizes.keys()), geneset_sizes.values())
with open("../data/target-dbs/all_dgi_geneset_sizes.csv", "w") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=geneset_sizes.keys())
    writer.writeheader()
    writer.writerow(geneset_sizes)
plt.show()
