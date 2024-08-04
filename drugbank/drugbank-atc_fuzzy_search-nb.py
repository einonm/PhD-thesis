# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.7.1
# ---

# # Drugbank fuzzy name search
#
# Lots of drugbank entries are missing ATC codes, where clearly an ATC code should exist (from a quick anecdotal check).
#
# So let's use a fuzzy text comparison algorithm and check how well the drug names from the ATC code list match the drugbank entries with missing ATC codes

# +
import os
import pandas as pd
from IPython.display import display

#!pip install fuzzywuzzy
#!pip install fuzzywuzzy[speedup]
import fuzzywuzzy

from fuzzywuzzy import fuzz
from fuzzywuzzy import process

# +
# Load up the ATC level 5 list (ATC code vs drug name)
data_path = "../data/drugbank/"


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


atc5 = get_atc_file(data_path + "../atc/level-5-atc.txt")
atc5.set_index("atc_code", inplace=True)

# +
# Load up the set of drugbank entries with protein targets but no ATC code

drugbank_noatc = pd.read_csv(
    data_path + "no_atc_code.csv",
    names=["db_id", "name", "atc_code", "targets"],
    sep="\t",
)

display(drugbank_noatc.head())
print(len(drugbank_noatc))

# +
# for each drug with no ATC code, find the top matches in the ATC list
match_list = pd.DataFrame(
    columns=["db_id", "db_name", "atc_match", "fuzz_score", "atc_code"]
)

for index, row in drugbank_noatc.iterrows():
    matches = process.extract(row["name"], atc5["description"])
    for df in matches:
        match_list = match_list.append(
            {
                "db_id": row["db_id"],
                "db_name": row["name"],
                "atc_match": df[0],
                "fuzz_score": df[1],
                "atc_code": df[2],
            },
            ignore_index=True,
        )
# -

match_list.set_index("db_id", inplace=True)
match_list.sort_values("fuzz_score", ascending=False, inplace=True)
display(match_list)
match_list.to_csv(data_path + "noatc_fuzzy-matches.txt", sep="\t")

# Testing
process.extract("ascorbic acid", atc5["description"])

# Testing
process.extract("vitamin c", atc5["description"])
