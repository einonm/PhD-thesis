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
from IPython.display import display, HTML, Latex, display_latex, display_html

# #!pip install venn
# %matplotlib inline
import matplotlib.pyplot as plt
import venn

# #!pip install seaborn
import seaborn as sns

# #!pip install pandas_profiling
import pandas_profiling


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

# Look at the collected data from each contributing database:

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

# +
# TODO - get a list of PubChem CIDs vs ATC codes (https://pubchem.ncbi.nlm.nih.gov/classification/#hid=79)
#                                                (https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest)
# However, if there are different compounds with the same ATC code, is it right to include all of them in the
# target lists, as only one form will be given?

# +
# G01AE10 (combinations of sulfonamides) has > 30 different drugs associated
# TODO - How to capture all the targets for drugs like this? Is it being done correctly?
# jointall_df[jointall_df['ATC_DrugBank'].str.contains('G01AE10')]
# -

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

atc_exploded_df[atc_exploded_df['all_ATC'] == 'L01DB02']

convert_list = atc_exploded_df[atc_exploded_df['all_ATC'] == 'L01DB02']['targets_entrez_DGIdb'][754].split(' ')
for entrez in convert_list:
    if entrez != 'nan':
        print(conversion_df.loc[int(entrez)])

convert_list = atc_exploded_df[atc_exploded_df['all_ATC'] == 'L01DB02']['targets_entrez_DsigDB'][754].split(' ')
for entrez in convert_list:
    if entrez != 'nan':
        print(conversion_df.loc[int(entrez)])

convert_list = atc_exploded_df[atc_exploded_df['all_ATC'] == 'L01DB02']['targets_entrez_DrugCentral'][754].split(' ')
for entrez in convert_list:
    if entrez != 'nan':
        print(conversion_df.loc[int(entrez)])

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

# ### Drugs with unknown ATC codes
#
# 'Telapre' is an error in the TDD cross-matching database, whilst other ATC codes do not exist in the current ATC database - typos?

different_set = set(atc_exploded_df["all_ATC"]) - combined_set

# None of these exist in the current official ATC set (??)
different_set

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

# used for thesis text
archive_atc_table(intermediate_atc_df).to_csv(
    "../data/target-dbs/intermediate_atc_targets.tsv", sep="\t", index=None
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
            if len(intermediate_atc_df.at[x, col_name]) >= 1
        ]
    )
    if "nan" in atcs[i]:
        atcs[i].remove("nan")
    set_names.append(DGIs[i] + ": " + str(len(atcs[i])))
    display(
        "From "
        + str(len(intermediate_atc_df))
        + " DUGGIE ATC drugs >= 1 target, "
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

# ### Venn diagram of ATC drug targets contributed by the top 5  databases

labels = venn.get_labels([atcs[2], atcs[3], atcs[0], atcs[1], atcs[6]])
fig, ax = venn.venn5(
    labels, names=[set_names[2], set_names[3], set_names[0], set_names[1], set_names[6]]
)

# ### Venn diagram of contributed ATC drugs from the bottom 5 contributing databases

labels = venn.get_labels([atcs[1], atcs[6], atcs[5], atcs[4], atcs[7]])
fig, ax = venn.venn5(
    labels, names=[set_names[1], set_names[6], set_names[5], set_names[4], set_names[7]]
)

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
        + " DUGGIE ATC drugs >= 5 targets, "
        + str(len(atcs[i]))
        + " have targets from "
        + DGIs[i]
    )
# -

# ###  Venn diagram of contributed ATC drugs from the top 5 contributing databases

labels = venn.get_labels([atcs[0], atcs[1], atcs[2], atcs[4], atcs[3]])
fig, ax = venn.venn5(
    labels, names=[set_names[0], set_names[1], set_names[2], set_names[4], set_names[3]]
)

# ### Venn diagram of contributed ATC drugs from the bottom 5 contributing databases

labels = venn.get_labels([atcs[4], atcs[3], atcs[5], atcs[6], atcs[7]])
fig, ax = venn.venn5(
    labels, names=[set_names[4], set_names[3], set_names[5], set_names[6], set_names[7]]
)

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

# ### Look at targets contributing to each drug entry

# How much of the final gene set is each DB capturing (N.B. drugs >=5 targets)
for dgi in DGIs:
    all_atc_df[dgi] = [
        int(
            len(all_atc_df.at[x, "targets_entrez_" + dgi])
            / all_atc_df.at[x, "len_all_entrez"]
            * 100
        )
        for x in all_atc_df.index
    ]

graph = sns.boxplot(data=all_atc_df[DGIs]).set(
    xlabel="Contributing database",
    ylabel="Target gene coverage %",
    title="Target gene coverage by contributing database for each ATC drug",
)

all_atc_df.describe()

# ### How many unique targets are used in total?

# How much of the final gene set is each DB uniquely capturing?
for name in ["targets_entrez_" + dgi for dgi in DGIs]:
    dgi_genes = [all_atc_df.at[x, name] for x in all_atc_df.index]
    all_list = []
    for item in dgi_genes:
        for elem in item:
            all_list.append(elem)
    display(name + ": " + str(len(set(all_list))))

dgi_genes = [all_atc_df.at[x, "all_targets_entrez"] for x in all_atc_df.index]
all_list = []
for item in dgi_genes:
    for elem in item:
        all_list.append(elem)
display("all_targets_entrez: " + str(len(set(all_list))))

# ## What percentage of all targets for a drug is unique to each DB?

# +
unique_df = all_atc_df

unique_df["targets_entrez_no_TTD"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_DrugBank"] = (
    unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_DGIdb"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_STITCH"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_GtoPdb"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_TTD"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_T3DB"]
)

unique_df["targets_entrez_no_T3DB"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_TTD"]
)

unique_df["targets_entrez_no_DsigDB"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_T3DB"]
    + unique_df["targets_entrez_DrugCentral"]
    + unique_df["targets_entrez_TTD"]
)

unique_df["targets_entrez_no_DrugCentral"] = (
    unique_df["targets_entrez_DrugBank"]
    + unique_df["targets_entrez_DGIdb"]
    + unique_df["targets_entrez_STITCH"]
    + unique_df["targets_entrez_GtoPdb"]
    + unique_df["targets_entrez_T3DB"]
    + unique_df["targets_entrez_DsigDB"]
    + unique_df["targets_entrez_TTD"]
)

for dgi in DGIs:
    unique_df["targets_entrez_no_" + dgi] = [
        set(x) for x in unique_df["targets_entrez_no_" + dgi]
    ]


# +
# TODO - this shorter version is buggy

# Now we need a set of columns that have all 'targets not including this DB' in them.
# We'll use it to see what targets are unique to each DB
# unique_df = all_atc_df

# for dgi in DGIs:
#    other_DGIs = DGIs.copy()
#    other_DGIs.remove(dgi)

# for other_dgi in other_DGIs:
#    temp = []
#    unique_df['targets_entrez_no_' + dgi] = [temp + x for x in unique_df['targets_entrez_' + other_dgi]]

# unique_df['targets_entrez_no_' + dgi] = [set(x) for x in unique_df['targets_entrez_no_' + dgi]]
# -

# (#Targets unique to DB) = #DBtargets - (#DBtargets also in other DBs)
# # %age is ((#Targets unique to DB) / total targets) * 100
for index, row in unique_df.iterrows():
    for DGI in DGIs:
        unique_df.at[index, DGI] = (
            len(
                set(unique_df.at[index, "targets_entrez_" + DGI])
                - set(unique_df.at[index, "targets_entrez_" + DGI]).intersection(
                    set(unique_df.at[index, "targets_entrez_no_" + DGI])
                )
            )
            / unique_df.at[index, "len_all_entrez"]
        ) * 100

graph = sns.boxplot(data=unique_df[DGIs]).set(
    xlabel="Contributing database",
    ylabel="Coverage unique to DB %",
    title="Unique target gene contribution by each database for each ATC drug ",
)

# **Convert the entrez IDs to Ensembl, to make comparison with GTEx QTL data later on**

conversion_df = pd.read_csv(
    "/home/mpnme/data/phd/wrangling/Ensembl_prot_entrez.tsv", sep="\t"
)
conversion_df.drop(["Ensembl_prot", "HGNC symbol"], axis=1, inplace=True)
conversion_df.dropna(subset=["Entrez ID"], inplace=True)
conversion_df["Entrez ID"] = conversion_df["Entrez ID"].astype(int)
conversion_df.drop_duplicates(subset="Entrez ID", inplace=True)
conversion_df.set_index("Entrez ID", inplace=True)
conversion_df.sort_index(inplace=True)


all_atc_df["Ensembl"] = np.empty((len(all_atc_df), 0)).tolist()
fail_list = []
for index, row in all_atc_df.iterrows():
    ensembl = []
    for entrez in row["all_targets_entrez"]:
        int_entrez = int(entrez)
        if int_entrez in conversion_df.index:
            ensembl.append(conversion_df.at[int_entrez, "Ensembl"])
        else:
            fail_list.append(entrez)
    all_atc_df.at[index, "Ensembl"] = ensembl

# We failed to convert these entrez into Ensembl
display(len(set(fail_list)))
display(len(fail_list))
set(fail_list)

all_atc_df["unconverted"] = all_atc_df["all_targets_entrez"].map(len) - all_atc_df[
    "Ensembl"
].map(len)

shortlist = all_atc_df[all_atc_df["unconverted"] != 0]
display(shortlist.head())
display(len(shortlist))
dummy = plt.hist(shortlist["unconverted"], bins=7)

# Count the total number of interactions per DB
for dgi in DGIs:
    print(
        dgi
        + " total interactions is : "
        + str(sum([len(x) for x in all_atc_df["targets_entrez_" + dgi]]))
    )
print(
    "DUGGIE total interactions is : "
    + str(sum([len(x) for x in all_atc_df["all_targets_entrez"]]))
)

# Count the total number of targets per DB
for dgi in DGIs:
    print(
        dgi
        + " total interactions is : "
        + str(sum([len(x) for x in all_atc_df["targets_entrez_" + dgi]]))
    )
print(
    "DUGGIE total interactions is : "
    + str(sum([len(x) for x in all_atc_df["all_targets_entrez"]]))
)

for index, row in all_atc_df.iterrows():
    all_atc_df.at[index, "gene_list"] = " ".join(
        [str(x) for x in sorted(list(set(all_atc_df.at[index, "Ensembl"])))]
    )

all_atc_df["len_all_ensembl"] = [len(x) for x in all_atc_df["Ensembl"]]
all_atc_df = all_atc_df[all_atc_df["len_all_ensembl"] >= 5]
final_df = all_atc_df[["atc_code", "gene_list"]]
final_df.index.sort_values()
display(final_df)

# Gene lists should be sorted for this to work
final_df.drop_duplicates(subset=["gene_list"], inplace=True)
display(final_df)

final_df.to_csv(
    "../data/target-dbs/all_dgi_targets_atc_ensembl.csv",
    sep="\t",
    header=None,
    index=None,
)

final_df[final_df["atc_code"] == "N06AA12"]
