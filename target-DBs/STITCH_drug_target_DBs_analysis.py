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

# create an stitch_targets_entrez/len_stitch_entrez columns of de-duped targets for an ATC drug table
def collate_stitch_targets(in_df):
    df = in_df.copy()
    # Convert each string of entrez IDs and ATC codes to a set list
    for name in ["targets_entrez_" + dgi for dgi in DGIs]:
        df[name] = [x.split() for x in df[name]]
        df[name] = [set(x) for x in df[name]]
        [x.remove("nan") if "nan" in x else x for x in df[name]]
        df[name] = [list(x) for x in df[name]]

    df["stitch_targets_entrez"] = df["targets_entrez_STITCH"]

    df["stitch_targets_entrez"] = [set(x) for x in df["stitch_targets_entrez"]]
    df["len_stitch_entrez"] = [len(x) for x in df["stitch_targets_entrez"]]

    return df


# ## Drug target DB analysis, STITCH only
#
# Here, we want to explore the database that we've created - looking at the contributing data sources and measuring the improvement each one provides - e.g. if no new targets or drugs are uniquely contributed by a database, it could be removed from the analysis.
#
# This analysis is for STITCH only, in order to compare it with the all 8 database version.
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
intermediate_single_df = collate_stitch_targets(single_atc_targets_df)
intermediate_mix_df = collate_stitch_targets(mixture_atc_targets_df)

intermediate_combi_df = collate_stitch_targets(combination_atc_targets_df)
intermediate_atc_df = pd.concat(
    [intermediate_combi_df, intermediate_single_df, intermediate_mix_df],
    ignore_index=True,
)

for index, row in intermediate_atc_df.iterrows():
    intermediate_atc_df.at[index, "entrez"] = " ".join(
        [str(x) for x in intermediate_atc_df.at[index, "stitch_targets_entrez"]]
    )

# Remove drugs with no targets, also do so for split tables (combi/mix/single) to check numbers
intermediate_combi_df = intermediate_combi_df[
    intermediate_combi_df["stitch_targets_entrez"] != set()
]
intermediate_single_df = intermediate_single_df[
    intermediate_single_df["stitch_targets_entrez"] != set()
]
intermediate_mix_df = intermediate_mix_df[
    intermediate_mix_df["stitch_targets_entrez"] != set()
]
intermediate_atc_df = intermediate_atc_df[
    intermediate_atc_df["stitch_targets_entrez"] != set()
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

# ### Remove combination ATC drug duplicates of individual drugs

# +
# Annotate combi drugs with individual ATC code target lists
final_combi_df = collate_stitch_targets(combination_atc_targets_df)
final_single_df = collate_stitch_targets(single_atc_targets_df)
final_mix_df = collate_stitch_targets(mixture_atc_targets_df)
final_combi_df["single_targets"] = ""

for index, row in final_combi_df.iterrows():
    result = final_single_df[final_single_df["atc_code"] == row["single ATC"]][
        "stitch_targets_entrez"
    ].values
    if len(result):
        final_combi_df.at[index, "single_targets"] = result[0]
    else:
        final_combi_df.at[index, "single_targets"] = set()

# set operations - difference between combination and single ATC code target lists for each combination drug
final_combi_df["difference"] = (
    final_combi_df["stitch_targets_entrez"] - final_combi_df["single_targets"]
)

# Number of combination drugs that don't match the corresponding individual ATC's targets list
# But we need to remove drugs with < 5 targets to compare properly here
final_combi_df2 = final_combi_df[final_combi_df["len_stitch_entrez"] >= 5]
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

stitch_atc_df = pd.concat(
    [final_combi_df, final_single_df, final_mix_df], ignore_index=True
)

for index, row in stitch_atc_df.iterrows():
    stitch_atc_df.at[index, "entrez"] = " ".join(
        [str(x) for x in stitch_atc_df.at[index, "stitch_targets_entrez"]]
    )

# Only look at drug with >= 5 targets
stitch_atc_df = stitch_atc_df[stitch_atc_df["len_stitch_entrez"] >= 5]

astitch_atc_df = stitch_atc_df.drop(["all_ATC", "member_atc_codes", "PubChem CID"], axis=1)
stitch_atc_df = stitch_atc_df.reset_index()
display(
    "There are "
    + str(len(stitch_atc_df))
    + " ATC coded drugs with 5 or more targets after removing duplicate combi drugs"
)
# -

# ### Look at targets contributing to each drug entry

# How much of the final gene set is each DB capturing (N.B. drugs >=5 targets)
for dgi in DGIs:
    stitch_atc_df[dgi] = [
        int(
            len(stitch_atc_df.at[x, "targets_entrez_" + dgi])
            / stitch_atc_df.at[x, "len_stitch_entrez"]
            * 100
        )
        for x in stitch_atc_df.index
    ]

graph = sns.boxplot(data=stitch_atc_df[DGIs]).set(
    xlabel="Contributing database",
    ylabel="Target gene coverage %",
    title="Target gene coverage by contributing database for each ATC drug",
)

stitch_atc_df.describe()

# ### How many unique targets are used in total?

# How much of the final gene set is each DB uniquely capturing?
for name in ["targets_entrez_" + dgi for dgi in DGIs]:
    dgi_genes = [stitch_atc_df.at[x, name] for x in stitch_atc_df.index]
    all_list = []
    for item in dgi_genes:
        for elem in item:
            all_list.append(elem)
    display(name + ": " + str(len(set(all_list))))

dgi_genes = [stitch_atc_df.at[x, "stitch_targets_entrez"] for x in stitch_atc_df.index]
all_list = []
for item in dgi_genes:
    for elem in item:
        all_list.append(elem)
display("stitch_targets_entrez: " + str(len(set(all_list))))


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


stitch_atc_df["Ensembl"] = np.empty((len(stitch_atc_df), 0)).tolist()
fail_list = []
for index, row in stitch_atc_df.iterrows():
    ensembl = []
    for entrez in row["stitch_targets_entrez"]:
        int_entrez = int(entrez)
        if int_entrez in conversion_df.index:
            ensembl.append(conversion_df.at[int_entrez, "Ensembl"])
        else:
            fail_list.append(entrez)
    stitch_atc_df.at[index, "Ensembl"] = ensembl

# We failed to convert these entrez into Ensembl
display(len(set(fail_list)))
display(len(fail_list))
set(fail_list)

stitch_atc_df["unconverted"] = stitch_atc_df["stitch_targets_entrez"].map(len) - stitch_atc_df[
    "Ensembl"
].map(len)

shortlist = stitch_atc_df[stitch_atc_df["unconverted"] != 0]
display(shortlist.head())
display(len(shortlist))
dummy = plt.hist(shortlist["unconverted"], bins=7)

for index, row in stitch_atc_df.iterrows():
    stitch_atc_df.at[index, "gene_list"] = " ".join(
        [str(x) for x in sorted(list(set(stitch_atc_df.at[index, "Ensembl"])))]
    )

stitch_atc_df["len_stitch_ensembl"] = [len(x) for x in stitch_atc_df["Ensembl"]]
stitch_atc_df = stitch_atc_df[stitch_atc_df["len_stitch_ensembl"] >= 5]
final_df = stitch_atc_df[["atc_code", "gene_list"]]
final_df.index.sort_values()
display(final_df)

# Gene lists should be sorted for this to work
final_df.drop_duplicates(subset=["gene_list"], inplace=True)
display(final_df)

final_df.to_csv(
    "../data/target-dbs/stitch_dgi_targets_atc_ensembl.csv",
    sep="\t",
    header=None,
    index=None,
)
