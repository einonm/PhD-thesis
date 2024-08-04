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

# ## Get PubMed references
#
# Given a list of drug names, this searches PubMed and scrapes all publications found with the term '(disease) AND (drug)'.

# +
import os
import pandas as pd
from IPython.display import display

import numpy as np

# #!pip install matplotlib
# %matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns

# #!pip install pymed
from pymed import PubMed

# #!pip install metapub
from metapub import PubMedFetcher

import requests
import json

# #!pip install beautifulsoup4
from bs4 import BeautifulSoup

# +
diseases = ['schizophrenia', 'alzheimer', 'major depressive disorder', 'parkinson', 'huntington', 'bipolar']

# drugnames files have been generated from the thesis chapter5 code, re-run this to regenerate
sz_df = pd.read_csv("../data/wrangling/sz_drugnames.txt", header=None, names=["drug"])
al_df = pd.read_csv("../data/wrangling/al_drugnames.txt", header=None, names=["drug"])
mdd_df = pd.read_csv("../data/wrangling/mdd_drugnames.txt", header=None, names=["drug"])
pd_df = pd.read_csv("../data/wrangling/pd_drugnames.txt", header=None, names=["drug"])
hd_df = pd.read_csv("../data/wrangling/hd_drugnames.txt", header=None, names=["drug"])
bp_df = pd.read_csv("../data/wrangling/bp_drugnames.txt", header=None, names=["drug"])
disease_dfs = [sz_df, al_df, mdd_df, pd_df, hd_df, bp_df]


# -

def get_publications_pymed(disease, drug):
    pubmed = PubMed(tool="jupyter", email="ei@nn.com")
    query = disease + " " + drug
    results = pubmed.query(query, max_results=500)

    drug_refs_df = pd.DataFrame()

    for article in results:
        # Show information about the article
        myDict = article.toDict()
        drug_refs_df = drug_refs_df.append(pd.Series(myDict), ignore_index=True)
        
    return drug_refs_df


def get_publications_metapub(disease, drug):
    query = "(" + disease + ") AND (" + drug + ")"
    fetch = PubMedFetcher()

    # get the  PMID for first 3 articles with keyword sepsis
    pmids = fetch.pmids_for_query(query, retmax=50)

    # get  articles
    articles = {}
    for pmid in pmids:
        articles[pmid] = fetch.article_by_pmid(pmid)   

    return articles


def get_publications_url(disease, drug):
    domain = 'https://pubmed.ncbi.nlm.nih.gov'
    # standard query
    queryLinkSearch = f'{domain}/?term=%28{disease}%29+AND+%28{drug}%29&format=pmid&sort=date&size=50'
    response = requests.get(queryLinkSearch)
    
    parsed_html = BeautifulSoup(response.content)
    result = parsed_html.body.find('div', attrs = {"search-results"}).text
    
    return result[1:-1].split('\r\n')


# +
#get_publications_pymed("alzheimer", "drug")
# -

get_publications_metapub("alzheimer", "bisoprolol")

get_publications_url("alzheimer", "bisoprolol")


def get_all_publications(disease, drugs_df):
    results_path = '../data/wrangling/publications/' + disease + '/'
    os.makedirs(results_path, exist_ok=True)
    
    # crap manual hack if pubmed website starts refusing requests - restart at failed count
    count = 1
    print(disease, ": length ", len(drugs_df))
    for drug in drugs_df["drug"]:
        print(drug, count)
        count = count + 1
        # replace 1 with failed count here
        if count > 77:
            # search for ATC code drug, in and format (including combination drugs)
            result_df = get_publications_pymed(disease, drug)
            result_df.to_csv(results_path + drug + '.csv', index=None)
            if 'and' in drug:
                drug_names = []
                num_combination_drugs = 2
                first_other_drugs = drug.split(' and ')
                if ',' in first_other_drugs[0]:
                    num_combination_drugs = 3
                    first_second_drugs = first_other_drugs[0].split(', ')
                    if ',' in first_second_drugs[0]:
                        num_combination_drugs = 4
                        # Form of 'Drug1, Drug2, Drug3 and Drug4'
                        first_second_third_drugs = first_second_drugs[0].split(', ')
                        drug_names.append(first_second_third_drugs[0])
                        drug_names.append(first_second_third_drugs[1])
                    else:
                        # Form of 'Drug1, Drug2 and Drug3'
                        drug_names.append(first_second_drugs[0])
                    drug_names.append(first_second_drugs[1])
                else:
                    # Form of 'Drug1 and Drug2'
                    drug_names.append(first_other_drugs[0])
                # Get whatever drug is at the end
                drug_names.append(first_other_drugs[1])
                
                for comb_count in range(0,num_combination_drugs):
                    print(comb_count, drug_names[comb_count])
                    result_df = get_publications_pymed(disease, drug_names[comb_count])
                    result_df.to_csv(results_path + drug_names[comb_count] + '.csv', index=None)


# crap manual hack if pubmed website starts refusing requests - restart at failed count
for count in range(5,len(diseases)):
    get_all_publications(diseases[count], disease_dfs[count])


