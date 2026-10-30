Markdown
# UPF1 Analysis 

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-blue.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the code and data associated with the manuscript: **"UPF1 Analysis: [Integrative Degradome Profiling and Machine Learning Reveal 3′UTR Determinants of Nonsense-Mediated Decay in Arabidopsis ]"**. 

## 📂 Repository Structure

```text
.
├── colab/                      # Language Model & Feature Importance
│   ├── d300_training_importance__test_human.ipynb
│   ├── u300d300_train_importance.ipynb
│   └── file/                   # Processed sequence datasets
│       ├── input_human_seq300.txt
│       ├── input_human_TRAM1_kmers_12.txt
│       ├── downstream_300_seq_Nointron.txt
│       └── PN_sequence_up_300_downstream_300_seq_Nointron.txt
├── R/                          # Visualization
│   ├── All_Figure_rda.r        # Script for figure generation
│   ├── UPF1_figures.Rproj      # R Project configuration
│   ├── data/                   # Processed .rda objects (Fig 2e - S5)
│   └── output/                 # Generated PDFs
└── README.md                   # Project documentation
🔬 File Descriptions
1. LM Finetune (/colab)

d300_training_importance__test_human.ipynb: Training script for the language model focusing on the 300 nt downstream of the Termination Codon (TC). It includes human sequence test.

u300d300_train_importance.ipynb: Training script utilizing upstream (300 nt) and downstream (300 nt) regions of the TC.

file/ Inputs:

input_human_seq300.txt: human sequence test set.

input_human_TRAM1_kmers_12.txt: Serial deletion 12-mer data for TRAM1 importance analysis.

downstream_300_seq_Nointron.txt: Training data for TC-downstream regions.

PN_sequence_up_300_downstream_300_seq_Nointron.txt: Training data for upstream and downstream sequence analysis.

2. Figures Visualization (/R)
Scripts for manuscript figures.

All_Figure_rda.r: A R script that reads the pre-processed .rda objects and renders them into the PDF panels.

data/: R data objects for Figures.

output/: Stores the PDF outputs corresponding to the manuscript.

🚀 Usage Instructions
To Reproduce Manuscript Figures:

Open the UPF1_figures.Rproj file in RStudio.

Execute the All_Figure_rda.r script.

To Run NMD-LM Models:
Upload the .ipynb files and the /file directory to Google Colab or Jupyter server.

Ensure the data files are in the expected relative path (/file/).

Run the notebook cells to perform language model fine-tuning or calculate sequence importance.

📊 Data Availability
Raw Sequencing Data: Deposited at NCBI Gene Expression Omnibus (GEO) under accession number [GSEXXXXXX].

Processed Data: All intermediate .rda files and input .txt files required for the analysis are included in this repository.

🖋️ Citation
If you use the data provided here, please cite the following paper:

[Full Author List]. (2026). UPF1 Analysis: [Full Paper Title]. [Journal Name]. DOI: [Journal DOI]