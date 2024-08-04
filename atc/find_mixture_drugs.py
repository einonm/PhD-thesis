# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ### Find mixture drugs
#
# We want to get separate lists of ATC codes that are:
#
# * singular drugs
# * mixture drugs of one or more substances
#
# It would be good to also construct a hierarchy of drug ATC codes that constitute each mixture ATC code.

import pandas as pd


# +
def get_atc_file(filename):
    file = open(filename, "r")
    lines = file.readlines()
    result = pd.DataFrame()
    for line in lines:
        code = line.split(" ")[0]
        desc = " ".join(line.split(" ")[1:]).strip()
        result = result.append(
            {"atc_code": code, "description": desc}, ignore_index=True
        )

    file.close()
    return result


atc5 = get_atc_file("../data/atc/level-5-atc.txt")
atc5["description"] = atc5["description"].str.lower()
# -

atc5

atc5_2 = atc5
# Strip DG numbers from descriptions
import re

regex = re.compile("(.*?)\[.*?\]")
for index, row in atc5_2.iterrows():
    result = re.findall(regex, atc5_2.at[index, "description"])
    if result != []:
        atc5_2.at[index, "description"] = result[0]

# +
# atc5_2

# +
mixture_acts_df = pd.DataFrame(columns=["description", "atc_code", "member_atc_codes"])
single_acts_df = pd.DataFrame(columns=["description", "atc_code"])
combination_atcs_df = pd.DataFrame(columns=["description", "atc_code"])
atc5_2["description"] = atc5_2["description"].str.strip()

# Find 'and', ',' and 'combinations' in descriptions to identify mixture drugs
for index, row in atc5_2.iterrows():
    if " and " in row["description"]:
        mixture_acts_df = mixture_acts_df.append(
            {"description": row["description"], "atc_code": row["atc_code"]},
            ignore_index=True,
        )
    elif "combinations" not in row["description"]:
        single_acts_df = single_acts_df.append(
            {"description": row["description"], "atc_code": row["atc_code"]},
            ignore_index=True,
        )
    else:
        combination_atcs_df = combination_atcs_df.append(
            {"description": row["description"], "atc_code": row["atc_code"]},
            ignore_index=True,
        )
# -

# Now save list of single drugs
single_acts_df.to_csv("../data/atc/single_drugs_atc5.tsv", sep="\t", index=None)
len(single_acts_df)

len(mixture_acts_df)

len(combination_atcs_df)

# +
# Split this table into up to 3 drugs at ',' and 'and' text points, first by converting ',' into 'and' so
# there is only one term to split on.
mixture_acts_df["split_desc"] = mixture_acts_df["description"].str.replace(",", " and ")
mixture_acts_df["split_desc"] = mixture_acts_df["split_desc"].str.lower()

mixture_acts_df["combi drug 3"] = [
    x[2] if len(x) > 2 else "" for x in mixture_acts_df["split_desc"].str.split("and")
]
mixture_acts_df["combi drug 2"] = [
    x[1] if len(x) > 1 else "" for x in mixture_acts_df["split_desc"].str.split("and")
]
mixture_acts_df["combi drug 1"] = [
    x[0] for x in mixture_acts_df["split_desc"].str.split(" and ")
]

# remove whitespace
mixture_acts_df["combi drug 1"] = mixture_acts_df["combi drug 1"].str.strip()
mixture_acts_df["combi drug 2"] = mixture_acts_df["combi drug 2"].str.strip()
mixture_acts_df["combi drug 3"] = mixture_acts_df["combi drug 3"].str.strip()

# At this point, manually curate the correct ATC codes for combination drugs
# dump the list of individual combi drugs into a list for easier curation
combi_drug_list = (
    mixture_acts_df["combi drug 1"].to_list()
    + mixture_acts_df["combi drug 2"].to_list()
    + mixture_acts_df["combi drug 3"].to_list()
)
# -

combi_drug_list = list(set(combi_drug_list))
combi_drug_list.remove("")
combi_drug_set = set(combi_drug_list)
len(combi_drug_set)

# +
# Find out if ATC codes already exist for these drugs
atc_code_guesses_df = pd.DataFrame(columns=["atc_code", "description", "combi drug"])
for drug in combi_drug_set:
    for index, row in atc5_2.iterrows():
        if drug in row["description"] and " and " not in row["description"]:
            row["combi drug"] = drug
            atc_code_guesses_df = atc_code_guesses_df.append(row, ignore_index=True)

len(atc_code_guesses_df)

# +
trimmed_df = atc_code_guesses_df
trimmed_df["drop"] = False

# Remove 'psycholeptics' and 'combinations', and entires which are not direct matches
for index, row in trimmed_df.iterrows():
    if "psycholeptics" in row["description"] or "combinations" in row["description"]:
        trimmed_df.at[index, "drop"] = True
    elif row["description"] != row["combi drug"]:
        trimmed_df.at[index, "drop"] = True

trimmed_df = trimmed_df[trimmed_df["drop"] == False]
trimmed_df = trimmed_df.drop(["drop"], axis=1)
# remove duplicates, but sort and keep first
trimmed_df.sort_values("description", inplace=True)
trimmed_df["atc_code"] = trimmed_df["atc_code"].astype(str)

trimmed_df.drop_duplicates(subset=["description", "combi drug"], inplace=True)
# -

trimmed_df

# +
mixture_acts_df["drop"] = False

# Now try and fill each combi drug set with these individual contributing ATC codes
for index, row in mixture_acts_df.iterrows():
    if row["combi drug 1"] in list(trimmed_df["description"]):
        mixture_acts_df.at[index, "member_atc_codes"] = trimmed_df[
            trimmed_df["description"] == row["combi drug 1"]
        ]["atc_code"].item()
        if row["combi drug 2"] in list(trimmed_df["description"]):
            mixture_acts_df.at[index, "member_atc_codes"] = (
                mixture_acts_df.at[index, "member_atc_codes"]
                + " "
                + trimmed_df[trimmed_df["description"] == row["combi drug 2"]][
                    "atc_code"
                ].item()
            )
            if row["combi drug 3"]:
                if row["combi drug 3"] in list(trimmed_df["description"]):
                    mixture_acts_df.at[index, "member_atc_codes"] = (
                        mixture_acts_df.at[index, "member_atc_codes"]
                        + " "
                        + trimmed_df[trimmed_df["description"] == row["combi drug 3"]][
                            "atc_code"
                        ].item()
                    )
                else:
                    combination_atcs_df = combination_atcs_df.append(
                        {
                            "description": row["description"],
                            "atc_code": row["atc_code"],
                        },
                        ignore_index=True,
                    )
                    mixture_acts_df.at[index, "drop"] = True
        else:
            combination_atcs_df = combination_atcs_df.append(
                {"description": row["description"], "atc_code": row["atc_code"]},
                ignore_index=True,
            )
            mixture_acts_df.at[index, "drop"] = True
    else:
        combination_atcs_df = combination_atcs_df.append(
            {"description": row["description"], "atc_code": row["atc_code"]},
            ignore_index=True,
        )
        mixture_acts_df.at[index, "drop"] = True
# -


final_mixture_acts_df = mixture_acts_df[~mixture_acts_df["drop"]]

final_mixture_acts_df.drop(
    ["combi drug 1", "combi drug 2", "combi drug 3", "drop", "split_desc"],
    axis=1,
    inplace=True,
)
final_mixture_acts_df

final_mixture_acts_df.to_csv("../data/atc/mixture_drugs_atc5.tsv", sep="\t", index=None)
len(final_mixture_acts_df)

# +
# TODO - create column of single drugs that are alternatives to combination drugs
combination_atcs_df["single ATC"] = ""
combination_atcs_df["single description"] = ""

for index2, row2 in combination_atcs_df.iterrows():
    for index, row in single_acts_df.iterrows():
        if row["description"] in row2["description"]:
            if combination_atcs_df.at[index2, "single ATC"] != "":
                # A better match will take precedence
                if row2["description"].startswith(row["description"]):
                    combination_atcs_df.at[index2, "single ATC"] = row["atc_code"]
                    combination_atcs_df.at[index2, "single description"] = row[
                        "description"
                    ]
            else:
                combination_atcs_df.at[index2, "single ATC"] = row["atc_code"]
                combination_atcs_df.at[index2, "single description"] = row[
                    "description"
                ]

# -

# And save all other combination/psycholeptics drugs
combination_atcs_df.to_csv(
    "../data/atc/combination_drugs_atc5.tsv", sep="\t", index=None
)
len(combination_atcs_df)

combination_atcs_df

single_acts_df[single_acts_df["description"] == "carbenoxolone"]
