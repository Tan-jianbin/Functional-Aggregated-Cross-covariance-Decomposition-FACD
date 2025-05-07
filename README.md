# Canonical Cross-Covariance Analysis of High-Dimensional Longitudinal Data

This repository contains the data and code associated with the paper:  **"Associating High-Dimensional Longitudinal Datasets through Canonical Cross-Covariance Analysis"**  by *Jianbin Tan* and *Pixu Shi*.

## 1. Data

The repository includes both simulated and real-world datasets used in our study. All data required to reproduce the results presented in the paper are provided.

### Data Description

The `Data/` directory contains raw data files, including blood samples from 36 participants who underwent a controlled exercise protocol following overnight fasting. Samples were collected at baseline, and at 2, 15, 30, and 60 minutes post-exercise. The plasma from these samples were profiled for untargeted proteomics by SWATH-MS, targeted proteomics by immunoassays, untargeted metabolomics using LC-MS, and semi-targeted lipidomics using Lipidyzer, capturing key dynamic molecular responses to acute physical
activity. The original dataset was obtained from [http://hmp2-data.stanford.edu/](http://hmp2-data.stanford.edu/).

## 2. Code

### Overview

This repository provides a complete implementation of the proposed Canonical Cross-Covariance Analysis (CCVA) framework. It includes scripts for simulation, data preprocessing, and analysis of real longitudinal datasets.

### Reproducibility

- **Preprocessing**  
  The script `data_preprocess.Rmd` performs preprocessing of the raw dataset.

- **Simulation Study**  
  Results in Section 5 of the paper can be reproduced by executing `Conduct_sim.R`.

- **Real Data Analysis**  
  The main analysis of the real dataset can be reproduced using `Data_analysis.Rmd`. Cross-validation results and comparative analyses are available in `Data_analysis_comp.R`.