# Microbial_Cohorts
This is the Github repository containing all necessary R-scripts to recreate the analyses from Milke et al. (2025)

Overview

This repository provides a complete R-based workflow for analyzing and comparing microbial community datasets across cohorts. The scripts implement data filtering, clustering, 
network analysis, compositional stability assessment, and various cohort-level comparisons used in the associated publication.

Contents

The repository is organized into individual scripts, each conducting an analysis shown in the original paper. Most scripts build on each other, so running them according to their numbers is mandatory.
The folder R/ contains scripts for certain functions. The main scripts depend on these functions. 

The scripts utilize so called datalist-formats, that wrap Count-Data and (environmental) Meta-Data into single objects. These objects can be manipulated using an analogue workflow as implemented by 
the well-known dplyr-package (see Github repository [ExCom](https://github.com/dermilke/ExCom)). The scripts download the necessary functions via the source() command directly from GitHub.

Prerequisites

Before running the scripts, ensure you have the following installed:

• R (≥ 4.0)

• Required R packages (tidyverse, vegan, igraph, ggplot2, SpiecEasi, Cluster, mclust)

• FastSpar installation via Conda-Environment

Additional specialized packages may be needed depending on script-specific dependencies.

Citation

If you use this repository for your research, please cite:

Milke et al. (2025) (in review).  ￼
