## This script runs the clustering pipeline for all 25 different datasets
# In the original publication we used two different correlation-coefficient thresholds (r = 0.25 and r = 0.4)
#
# In this script we use the automated version of the clustering-pipeline and
# use pre-defined clustering-parameters, such as expected number of clusters,
# consensus-approach, and clustering-resolution interval used for consensus analysis.
#
# These values were manually inferred by running the semi-automatic clustering pipeline-mode.
# These values are the same that we used in the original publication.
#
# We implemented a parallelization to improve computation runspeed. Since this
# script is heavy on I/O operations, adjust the number of workers carefully.
#
# Author: Felix Milke
# Date: 26.02.2025

library(tidyverse)
library(igraph)
library(cluster)
library(mclust)
library(furrr)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Import_SparCC_Network.R")
source("R/Cluster_Pipeline.R")

# Define number of threads for parallelization here
plan(multisession, workers = 6)

#### Params ####

# Define r-threshold (in publication we used 0.25 or 0.4 but all values 0-1 are possible)
r_threshold = 0.25

# Iterate clustering over all FastSpar-correlation matrices in the output-folder 
# (one for each dataset)
datasets <- list.files("output/fastspar/output_files/", pattern = "cor")

clustering_datasets_0.25 <- furrr::future_map(datasets, function(datasets) {
  
  if (grepl(pattern = "Comparison", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 2)            # for r = 0.25
      cluster_num <- 8                   # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.7, 1.6)            # for r = 0.4
      cluster_num <-  10                   # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Comparison_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Comparison_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Comparison")
    
  } else if (grepl(pattern = "CalCOFI", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.7, 1.6)            # for r = 0.25
      cluster_num <- 6                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1.5, 2.5)            # for r = 0.4
      cluster_num <- 5                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_CalCOFI_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_CalCOFI_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "CalCOFI")
    
  } else if (grepl(pattern = "SPOT", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 1.8)            # for r = 0.25
      cluster_num <- 6                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.8, 2)            # for r = 0.4
      cluster_num <- 5                    # for r = 0.4
      cluster_method <- "Cluster_consensus" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_SPOT_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_SPOT_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "SPOT")
  
  } else if (grepl(pattern = "Soil_Australia", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.9, 1.7)            # for r = 0.25
      cluster_num <- 7                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1.2, 2)            # for r = 0.4
      cluster_num <- 9                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Soil_Australia_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Soil_Australia_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Soil_Australia")
    
  } else if (grepl(pattern = "Lake_Mendota", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1.3, 2)            # for r = 0.25
      cluster_num <- 6                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 1.8)            # for r = 0.4
      cluster_num <- 10                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Lake_Mendota_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Lake_Mendota_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Lake_Mendota")
    
    
  } else if (grepl(pattern = "Distalgut_Human", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.7, 2)            # for r = 0.25
      cluster_num <- 10                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 3)            # for r = 0.4
      cluster_num <- 6                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Distalgut_Human_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Distalgut_Human_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)

    cluster_result <- c(cluster_result, Type = "Distalgut_Human")
    
  } else if (grepl(pattern = "Distalgut_Monkey", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 1.5)            # for r = 0.25
      cluster_num <- 5                    # for r = 0.25
      cluster_method <- "Cluster_consensus" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 1.8)            # for r = 0.4
      cluster_num <- 5                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Distalgut_Monkey_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Distalgut_Monkey_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Distalgut_Monkey")
    
  } else if (grepl(pattern = "Distalgut_Kangaroo", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.9, 2)            # for r = 0.25
      cluster_num <- 7                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 2)            # for r = 0.4
      cluster_num <- 10                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Distalgut_Kangaroo_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Distalgut_Kangaroo_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Distalgut_Kangaroo")
    
  } else if (grepl(pattern = "Distalgut_Rabbit", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.9, 2)            # for r = 0.25
      cluster_num <- 6                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1.1, 1.9)            # for r = 0.4
      cluster_num <- 8                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Distalgut_Rabbit_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Distalgut_Rabbit_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Distalgut_Rabbit")
    
  } else if (grepl(pattern = "Distalgut_Deer", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1.1, 2)            # for r = 0.25
      cluster_num <- 8                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 1.6)            # for r = 0.4
      cluster_num <- 6                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Distalgut_Deer_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Distalgut_Deer_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Distalgut_Deer")
    
  } else if (grepl(pattern = "Human_Saliva", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.8, 1.9)            # for r = 0.25
      cluster_num <- 7                    # for r = 0.25
      cluster_method <- "Cluster_consensus" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 1.5)            # for r = 0.4
      cluster_num <- 5                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Human_Saliva_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Human_Saliva_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Human_Saliva")
    
  } else if (grepl(pattern = "Nose_Secretion", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.7, 1.4)            # for r = 0.25
      cluster_num <- 4                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 1.8)            # for r = 0.4
      cluster_num <- 6                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Nose_Secretion_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Nose_Secretion_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Nose_Secretion")
    
  } else if (grepl(pattern = "Plant_Rhizosphere", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.6, 1.5)            # for r = 0.25
      cluster_num <- 6                    # for r = 0.25
      cluster_method <- "Cluster_consensus" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 1.5)            # for r = 0.4
      cluster_num <- 8                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Plant_Rhizosphere_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Plant_Rhizosphere_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Plant_Rhizosphere")
    
  } else if (grepl(pattern = "Plant_Roots", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 2)            # for r = 0.25
      cluster_num <- 5                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 2)            # for r = 0.4
      cluster_num <- 6                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Plant_Roots_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Plant_Roots_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Plant_Roots")
    
  } else if (grepl(pattern = "Plant_Surface", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.9, 1.8)            # for r = 0.25
      cluster_num <- 7                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 1.8)            # for r = 0.4
      cluster_num <- 5                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Plant_Surface_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Plant_Surface_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Plant_Surface")
    
  } else if (grepl(pattern = "Sediment_Nonsaline", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.9, 1.5)            # for r = 0.25
      cluster_num <- 5                    # for r = 0.25
      cluster_method <- "Cluster_consensus" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 1.7)            # for r = 0.4
      cluster_num <- 4                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Sediment_Nonsaline_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Sediment_Nonsaline_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Sediment_Nonsaline")
    
  } else if (grepl(pattern = "Sediment_Saline", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.7, 1.5)            # for r = 0.25
      cluster_num <- 3                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.8, 1.8)            # for r = 0.4
      cluster_num <- 4                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Sediment_Saline_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Sediment_Saline_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Sediment_Saline")
    
  } else if (grepl(pattern = "Surface_Foot", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 2)            # for r = 0.25
      cluster_num <- 9                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 2)            # for r = 0.4
      cluster_num <- 9                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Surface_Foot_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Surface_Foot_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Surface_Foot")
    
  } else if (grepl(pattern = "Surface_Hand", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 2.7)            # for r = 0.25
      cluster_num <- 7                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.8, 1.5)            # for r = 0.4
      cluster_num <- 9                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Surface_Hand_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Surface_Hand_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Surface_Hand")
    
  } else if (grepl(pattern = "Surface_Nonsaline", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.9, 1.8)            # for r = 0.25
      cluster_num <- 7                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1, 2)            # for r = 0.4
      cluster_num <- 10                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Surface_Nonsaline_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Surface_Nonsaline_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Surface_Nonsaline")
    
  } else if (grepl(pattern = "Soil_Cultivated", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.5, 1.5)            # for r = 0.25
      cluster_num <- 3                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 1.6)            # for r = 0.4
      cluster_num <- 3                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Soil_Cultivated_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Soil_Cultivated_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Soil_Cultivated")
    
  } else if (grepl(pattern = "Soil_Field", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1.5, 2.2)            # for r = 0.25
      cluster_num <- 7                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1.2, 2)            # for r = 0.4
      cluster_num <- 8                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Soil_Field_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Soil_Field_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Soil_Field")
    
  } else if (grepl(pattern = "Soil_Garden", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.7, 1.3)            # for r = 0.25
      cluster_num <- 3                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.8, 1.5)            # for r = 0.4
      cluster_num <- 4                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Soil_Garden_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Soil_Garden_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Soil_Garden")
    
  } else if (grepl(pattern = "Soil_Sandfilter", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(0.9, 1.6)            # for r = 0.25
      cluster_num <- 6                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(1.2, 1.8)            # for r = 0.4
      cluster_num <- 9                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Soil_Sandfilter_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Soil_Sandfilter_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Soil_Sandfilter")
    
  } else if (grepl(pattern = "Water_Nonsaline_Germany", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 1.2)            # for r = 0.25
      cluster_num <- 6                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.8, 1.2)            # for r = 0.4
      cluster_num <- 5                    # for r = 0.4
      cluster_method <- "Cluster_consensus" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Water_Nonsaline_Germany_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Water_Nonsaline_Germany_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Water_Nonsaline_Germany")
    
  } else if (grepl(pattern = "Water_Nonsaline_Timeseries", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1, 1.5)            # for r = 0.25
      cluster_num <- 5                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 1.7)            # for r = 0.4
      cluster_num <- 6                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Water_Nonsaline_Timeseries_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Water_Nonsaline_Timeseries_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Water_Nonsaline_Timeseries")
    
  } else if (grepl(pattern = "Water_Nonsaline_USA", x = datasets)) {
    
    if (r_threshold == 0.25) {
      res_param <- c(1.1, 1.6)            # for r = 0.25
      cluster_num <- 8                    # for r = 0.25
      cluster_method <- "Cluster_network" # for r = 0.25
    } else if (r_threshold == 0.4) {
      res_param <- c(0.9, 1.8)            # for r = 0.4
      cluster_num <- 5                    # for r = 0.4
      cluster_method <- "Cluster_network" # for r = 0.4
    }
    
    cluster_result <- import_sparcc_network(cor_file = "output/fastspar/output_files/fastspar_cor_Water_Nonsaline_USA_table.tsv", 
                                            pval_file = "output/fastspar/output_files/fastspar_pvalues_Water_Nonsaline_USA_table.tsv", 
                                            min_r = 0, min_p = 0.05) %>%
      filter(weight >= r_threshold) %>%
      complete_cluster_pipeline(find_parameters = FALSE, res_param = res_param, 
                                cluster_method = cluster_method, cluster_num = cluster_num)
    
    cluster_result <- c(cluster_result, Type = "Water_Nonsaline_USA")
    
  }
  
  return(cluster_result)
  
  }, .options = furrr_options(seed = 123))

# We save results as R-objects depending on the applied r-threshold (here r = 0.25)
clustering_datasets_0.25 <- purrr::map(clustering_datasets_0.25, function(x) c(x, r_threshold = 0.25))
save(clustering_datasets_0.25, file = "output/clustering_datasets_025.RData")

# (here r = 0.4)
clustering_datasets_0.4 <- purrr::map(clustering_datasets_0.4, function(x) c(x, r_threshold = 0.4))
save(clustering_datasets_0.4, file = "output/clustering_datasets_04.RData")
