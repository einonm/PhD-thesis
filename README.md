# A Genomic Medicine Approach to Identifying Novel Drugs
--------------------------------------------------------

A staff PhD project - Mark Einon

## Abstract

Gene set enrichment analysis is an established method of identifying sub-sets
of related genes that are significantly associated with a disease whose
association with human genetic variation has been explored in a genome-wide
association study (GWAS). This method, alongside mining the ever increasing
wealth of experimentally obtained drug-target gene sets and functional genomic
datasets available, gives rise to the possibility of identifying new
repurposing opportunities for approved drugs. Such repurposing opportunities
are particularly sought after as therapeutic solutions for debilitating
neuropsychiatric and neurodegenerative diseases, for which available treatments
are scarce.

This thesis describes the creation, evaluation and use of an in-silico gene set
enrichment pipeline to discover novel evidence for repurposing already approved
drugs to treat six neuropsychiatric and neurodegenerative diseases -
schizophrenia, major depressive disorder, bipolar disorder, Alzheimer’s
disease, Parkinson’s disease and Huntington’s disease. Several strategies were
investigated with the aim of improving the yield and accuracy of pipeline
results - curating the most comprehensive dataset of drug-gene interactions
possible and exploring the use of more relevant functional QTL data with which
to annotate genes using available disease GWAS results. The application of
these strategies were evaluated on better understood diseases,
hypercholesterolemia and hypertension, for which drug treatments are more
plentiful.

The strategies investigated were found to be beneficial in identifying a
greater number of significantly disease-associated drugs, with nuanced results
leading to the adoption of a further strategy of running a battery of analyses
and ranking identified candidate repurposing drugs by the number of analyses in
which they appear. Executing this battery of analyses gave extensive results,
including a number of novel repurposing opportunities, for all six
neuropsychiatric and neurodegenerative diseases studied.

## Installation and usage

To use these scripts:

* Install Fedora 38 in a Virtual Machine (https://download.fedoraproject.org/pub/fedora/linux/releases/38/Workstation/x86_64/iso/Fedora-Workstation-Live-x86_64-38-1.6.iso)
* Open a bash terminal, and run the commands:
    * `$ sudo dnf install vim git-all meld R openldap-clients openldap-devel texlive-scheme-full`
    * `$ mkdir ~/source ~/source/code ~/source/code/data ~/source/magma-perms-analysis; cd ~/source`
* Copy SF7-results-data_archive.tar.gz to ~/source/code/data/, SF6-appendix-B_archive.zip to ~/source/magma-perms-analysis and SF5-code_archive.zip thesis supplimentary files to ~/source on the Fedora 38 VM
    * `$ cd ~/source/code/data; tar xvf SF7-results-data_archive.tar.gz`
    * `$ cd ~/source/code; unzip SF5-code_archive.zip`
    * `$ cd ~/source/appendix-b; unzip SF6-appendix-B_archive.zip`
    * `$ cd ~/source; wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh`
    * `$ sudo bash ./Anaconda3-2022.05-Linux-x86_64.sh # Install to /opt/anaconda3/ when prompted`
    * `$ git clone https://gitlab.com/einonm/jupyternb-setup.git`
    * `$ cd jupyternb-setup/gnome3-launcher/`
    * `$ ./install-jupyter-launcher-fedora.sh`
    * `$ run_jupyter_notebook_server.sh -n -e ~/source/code/environment.yml`
    * `$ sudo Rscript ~/source/code/install_qvalue.R`
* Copy libstdc++.so.6.0.30 or above from /usr/lib64/ to ~/jupyter-notebooks/lib (to prevent GLIBCXX errors):
    * `$ cp /usr/lib64/libstdc++.so.6.0.3x ~/jupyter-notebooks/lib/; cd ~/jupyter-notebooks/lib/`
    * `$ rm libstdc++.so.6 libstdc++.so`
    * `$ ln -s libstdc++.so.6.0.3x libstdc++.so.6`
    * `$ ln -s libstdc++.so.6.0.3x libstdc++.so`

### Generate thesis PDF

* In the Fedora VM after the above installation, run the 'Jupyter Notebook' program (search in 'Activities').
* In the Jupyter notebook browser window that appears, navigate to source/code/documentation/thesis/tex/.
* In the Jupyter browser window, navigate through each sub-directory, open each '.py' file as a jupyter notebook and run using the 'Cell -> Run All' menu option, and save each notebook when complete.
* Run:
    * `$ cd ~/source/code/documentation/thesis`
    * `$ ./scripts/generate_thesis.sh`

The finished thesis PDF appears as ~/source/code/documentation/thesis/tex/thesis/thesis.pdf

## Contributing
Contributing issues, feature proposals or code is actively welcomed - please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for more details.

## Code of Conduct
We want to create a welcoming environment for everyone who is interested in contributing. Please see the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md) file to learn more about our commitment to an open and welcoming environment.

## Copyright and Licence Information

See the [LICENSE file](LICENSE).
