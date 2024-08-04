# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# #!pip install psycopg2-binary
import psycopg2
import pandas as pd
import pandas.io.sql as sqlio

# A STITCH database copy is running on rocks, it has three schemas - 'items', 'evidence' and 'network'.
#
# Let's look at the drugs in items.chemicals table, where items are marked as being drugs, or not and are also listed in the chemicals_sources table, which lists their various identifiers.

item_conn = psycopg2.connect(
    user="postgres",
    password="treehouse",
    host="localhost",
    port="5432",
    database="evidence",
)

# TODO - add item_id_b for interacting pubchem IDs.

chem_sql = "SELECT flat_chemical_id,source_id,items.chemicals.preferred_name,mode,sources,species_id,item_id_b \
            FROM items.chemicals_sources, items.chemicals,evidence.actions_sets,evidence.sets_items \
            WHERE drug='t' \
            AND source_name='ATC' \
            AND items.chemicals_sources.flat_chemical_id=items.chemicals.chemical_id \
            AND evidence.actions_sets.item_id_a=items.chemicals.chemical_id \
            AND evidence.sets_items.item_id=items.chemicals.chemical_id \
            AND evidence.sets_items.set_id=evidence.actions_sets.sources;"

chem_df = sqlio.read_sql_query(chem_sql, item_conn)

display(chem_df.shape)
display(chem_df.head())

# De-duplicate this list, in order to reduce the number of chemical ID queries we need to run.

# The 'sources' column contains the set_id of the evidence for the link. We need to ensure that this is human (species_id=9606, or not specific to any species (-1))

# +
# for index,row in all_evidence_df.iterrows():
#    species_sql = "select species_id from evidence.sets_items where item_id='" \
#       + str(row['item_id_a']) + "' and set_id='"  + str(row['sources']) + "';"
#    species = sqlio.read_sql_query(species_sql, evidence_conn)
#    if len(species['species_id']) != 0:
#        all_evidence_df.at[index, 'species_id'] = species['species_id'][0];

# +
# sets_sql = "select * from evidence.sets_items where set_id='PRED292.21'"

# +
# sqlio.read_sql_query(sets_sql, evidence_conn)
# -

chem_df.to_csv("~/all_evidence_df.tsv", sep="\t")
