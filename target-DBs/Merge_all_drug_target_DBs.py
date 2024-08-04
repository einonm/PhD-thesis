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

import os
import numpy as np
import pandas as pd
from IPython.display import display

# #!pip install venn
# %matplotlib inline
import venn

# We want to see what the union of various drug target databases looks like.
#
# ## TTD

ttd_df = pd.read_csv("../data/target-dbs/ttd_ATC_targets.csv", sep="\t")
display(len(ttd_df))
ttd_df.head()

# Shouldn't be here... comes from a mistake in the TTD crossmatching file
ttd_df[ttd_df["ATC"] == "Telapre"]

# ## DrugBank

drugbank_df = pd.read_csv("../data/target-dbs/drugbank_atc_targets.csv", sep="\t")
drugbank_df = drugbank_df.rename(
    columns={
        "PubChem Compound": "PubChem CID",
        "atc_codes": "ATC",
        "genes": "targets_entrez",
    }
)
display(len(drugbank_df))
drugbank_df.head()

joint2_df = pd.merge(
    ttd_df, drugbank_df, how="outer", on="PubChem CID", suffixes=("_TTD", "_DrugBank")
)

display(len(joint2_df))
display(joint2_df.head())

# ## DGIdb

dgidb_df = pd.read_csv("../data/target-dbs/dgidb_targets.csv", sep="\t")
dgidb_df = dgidb_df.rename(
    columns={"pubchem_cid": "PubChem CID", "gene_list": "targets_entrez_DGIdb"}
)
display(dgidb_df.shape)
dgidb_df.head()

joint3_df = pd.merge(joint2_df, dgidb_df, how="outer", on="PubChem CID")

joint3_df.head(100)

joint3_df.shape

joint3_df[joint3_df['ATC_TTD'] == 'L01XA01']

# ## STITCH

stitch_df = pd.read_csv("../data/target-dbs/stitch_targets2.csv", sep="\t")
stitch_df = stitch_df.rename(
    columns={"PubChem_CID": "PubChem CID", "gene_list": "targets_entrez_STITCH"}
)
display(len(stitch_df))
stitch_df.head()

joint4_df = pd.merge(joint3_df, stitch_df, how="outer", on="PubChem CID")

# Now We have TTD, DrugBank, DGIdb and Stitch:

joint4_df.shape

joint4_df[joint4_df['ATC_TTD'] == 'L01XA01']

# ## Guide to Pharmacogenetics

gtopdb_df = pd.read_csv("../data/target-dbs/gtopdb_targets.csv", sep="\t")
gtopdb_df = gtopdb_df.rename(
    columns={"drug_pubchem_cid": "PubChem CID", "gene_list": "targets_entrez_GtoPdb"}
)
display(len(gtopdb_df))
gtopdb_df.head()

joint5_df = pd.merge(joint4_df, gtopdb_df, how="outer", on="PubChem CID")

joint5_df.fillna(value=np.nan, inplace=True)
joint5_df.shape

# ### T3DB

t3db_df = pd.read_csv("../data/target-dbs/T3DB/t3db_targets.tsv", sep="\t")
t3db_df = t3db_df.rename(
    columns={"pubchem_cid": "PubChem CID", "gene_list": "targets_entrez_T3DB"}
)
display(len(t3db_df))
t3db_df.head()

joint6_df = pd.merge(joint5_df, t3db_df, how="outer", on="PubChem CID")
joint6_df.fillna(value="nan", inplace=True)
joint6_df["PubChem CID"] = joint6_df["PubChem CID"].astype(str)
joint6_df.shape


joint6_df.head()

# Typical issue - same ATC code, different PubChem CIDs ??
joint6_df[joint6_df["ATC_TTD"] == "C01EB10"]

# ### DsigDB

dsigdb_df = pd.read_csv("../data/target-dbs/DsigDB/dsigdb_targets.tsv", sep="\t")
dsigdb_df = dsigdb_df.rename(
    columns={"atc_code": "ATC_DsigDB", "Entrez ID": "targets_entrez_DsigDB"}
)
dsigdb_df["targets_entrez_DsigDB"] = dsigdb_df["targets_entrez_DsigDB"].astype(str)
display(len(dsigdb_df))
dsigdb_df.head()

# TODO We should merge on something better than just TTD ATC codes
joint7_df = pd.merge(
    joint6_df, dsigdb_df, how="outer", left_on="ATC_TTD", right_on="ATC_DsigDB"
)
joint7_df.fillna(value="nan", inplace=True)
joint7_df.shape

joint7_df[joint7_df["ATC_DsigDB"] != "nan"]

# ### Drug Central

dc_df = pd.read_csv("../data/target-dbs/DrugCentral/drugcentral_targets.tsv", sep="\t")
dc_df = dc_df.rename(
    columns={"atc_code": "ATC_DrugCentral", "Entrez ID": "targets_entrez_DrugCentral"}
)
dc_df["targets_entrez_DrugCentral"] = dc_df["targets_entrez_DrugCentral"].astype(str)
display(len(dc_df))
dc_df.head()

# TODO We should merge on something better than just TTD ATC codes (but does doing an outer negate this?)
joint8_df = pd.merge(
    joint7_df, dc_df, how="outer", left_on="ATC_TTD", right_on="ATC_DrugCentral"
)
joint8_df.fillna(value="nan", inplace=True)
joint8_df.shape

joint8_df

joint8_df.to_csv("../data/target-dbs/jointall_targets-8dbs.tsv", sep="\t", index=None)
