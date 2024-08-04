# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import os
import numpy as np
import pandas as pd
from IPython.display import display

# #!pip install neo4j
# #!pip install neo4j-driver
from neo4j import GraphDatabase

# -

# We're going to create a graph of Drug-Target Interaction (DTI) databases and their sources

# +
# Add Drug-Target databases
dti_dbs = "all_target_dbs.csv"

driver = GraphDatabase.driver("bolt://localhost", auth=("guest", "guest"))
result = driver.session().run(
    "LOAD csv FROM 'file:///"
    + dti_dbs
    + "' AS line FIELDTERMINATOR ',' \
                              MERGE (d:dti_db {name:line[0], site:line[1]})"
)
# +
# Add Drug-Gene Interaction data sources
dti_sources = "all_dti_sources.csv"

driver = GraphDatabase.driver("bolt://localhost", auth=("guest", "guest"))
result = driver.session().run(
    "LOAD csv FROM 'file:///"
    + dti_sources
    + "' AS line FIELDTERMINATOR ',' \
                              MERGE (s:dti_source {name:line[0]})"
)


# +
# Add DTI DB relationships to sources
dti_collects = "all_dti_src_rels.csv"

driver = GraphDatabase.driver("bolt://localhost", auth=("guest", "guest"))
result = driver.session().run(
    "LOAD csv FROM 'file:///"
    + dti_collects
    + "' AS line FIELDTERMINATOR ',' \
                              MATCH (d:dti_db {name:line[0]}) \
                              MATCH (s:dti_source {name:line[1]}) \
                              MERGE (s)-[:DTI_source]->(d)"
)

# +
# Add DTI DB relationships to other DBs
dti_includes = "all_dti_includes.csv"

driver = GraphDatabase.driver("bolt://localhost", auth=("guest", "guest"))
result = driver.session().run(
    "LOAD csv FROM 'file:///"
    + dti_includes
    + "' AS line FIELDTERMINATOR ',' \
                              MATCH (d1:dti_db {name:line[0]}) \
                              MATCH (d2:dti_db {name:line[1]}) \
                              MERGE (d2)-[:included_in]->(d1)"
)
# +
driver = GraphDatabase.driver("bolt://localhost", auth=("guest", "guest"))
result = driver.session().run("MERGE (d:duggie {name:'DUGGIE', site:'none'})")

# Add my DTI DB relationships to other DBs (duggie)
dti_includes = "duggie_includes.csv"

driver = GraphDatabase.driver("bolt://localhost", auth=("guest", "guest"))
result = driver.session().run(
    "LOAD csv FROM 'file:///"
    + dti_includes
    + "' AS line FIELDTERMINATOR ',' \
                              MATCH (d1:duggie {name:line[0]}) \
                              MATCH (d2:dti_db {name:line[1]}) \
                              MERGE (d2)-[:included_in]->(d1)"
)
# -
