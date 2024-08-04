#! /usr/bin/env bash
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .sh
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.6
#   kernelspec:
#     display_name: Bash
#     language: bash
#     name: bash
# ---


# Takes one parameter, the notebook file to view

export PATH=/opt/anaconda3/bin:$PATH

source activate ~/jupyter-notebooks/

jupyter nbconvert --to slides "$1" --reveal-prefix ../reveal.js --output index
mv index.slides.html index.html
