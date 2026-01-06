## This script runs the cluster-verification including network-attack analysis
# Here, we first load the clustering-results for all datasets and both r-thresholds into the environment
# and then iterate over all elements of the list of datasets to infer cluster-significance and stability.
# See the verify_network_cluster() script for further details.
#
# This script uses parallelization. Adjust the number of threads manually.
#
# Author: Felix Milke
# Date: 26.02.2025

library(tidyverse)
library(igraph)
library(mclust)
library(furrr)

source("R/Cluster_Verification.R")

# This script contains new names for datasets and color-vectors for plotting
source("R/Names_Colors_Improvement.R")

load(file = "output/clustering_datasets_025.RData")
load(file = "output/clustering_datasets_04.RData")

# Adjust number of threads here as "workers"
plan(multisession, workers = 6)

# Run the verification and stability analysis pipeline on each dataset (here: r_threshold = 0.4)
cluster_verification_0.4_noderemoval <- furrr::future_map(clustering_datasets_0.4, function(cluster_data) {
  
  verify_network_cluster(cluster_data$Network, cluster_data$Cluster_Tab, remove_nodes = F,
                         res_param = mean(cluster_data$Resolutionrange), bootstrap_num = 50, random_expectation = T,
                         niter_rewire = V(cluster_data$Network)*10, attack_bootstrap = 50, dataset_name = cluster_data$Type)
  
}, .progress = T)

# Run the verification and stability analysis pipeline on each dataset (here: r_threshold = 0.25)
cluster_verification_0.25_noderemoval <- furrr::future_map(clustering_datasets_0.25, function(cluster_data) {
  
  verify_network_cluster(cluster_data$Network, cluster_data$Cluster_Tab, remove_nodes = F,
                         res_param = mean(cluster_data$Resolutionrange), bootstrap_num = 50, random_expectation = T,
                         niter_rewire = V(cluster_data$Network)*10, attack_bootstrap = 50, dataset_name = cluster_data$Type)
  
}, .progress = T)

plan(sequential)

saveRDS(cluster_verification_0.4_noderemoval, file = "output/Cluster_Verification_r04.rds")
saveRDS(cluster_verification_0.25_noderemoval, file = "output/Cluster_Verification_r025.rds")

cluster_verification_0.4_noderemoval <- readRDS("output/Cluster_Verification_r04.rds")
cluster_verification_0.25_noderemoval <- readRDS("output/Cluster_Verification_r025.rds")

# Combine both verification-datasets for r-treshold 0.25 and 0.4.
# Compute differences between cluster stability of observed and rewired networks 
# Use average ARI-values instead of individual bootstrap-iterations for each attack-strength and both r-thresholds
# Add the dataset-categories and new names for the datasets to improve visualization
cluster_verification_combined <- purrr::map(seq(1, length(cluster_verification_0.25_noderemoval)), function(x) { 
  cluster_verification_0.25_noderemoval[[x]]$Cluster_Attack_Table %>%
    group_by(Attack_Strength) %>%
    summarize(ARI_difference = (mean(ARI[Type == "Observed"]) - mean(ARI[Type == "Rewired"]))) %>%
    mutate(Type = clustering_datasets_0.25[[x]]$Type)
}) %>%
  bind_rows() %>%
  mutate(r_threshold = "r-threshold: 0.25") %>%
  bind_rows(., 
            purrr::map(seq(1, length(cluster_verification_0.4_noderemoval)), function(x) { 
              cluster_verification_0.4_noderemoval[[x]]$Cluster_Attack_Table %>%
                group_by(Attack_Strength) %>%
                summarize(ARI_difference = (mean(ARI[Type == "Observed"]) - mean(ARI[Type == "Rewired"]))) %>%
                mutate(Type = clustering_datasets_0.4[[x]]$Type)
            }) %>%
              bind_rows() %>%
              mutate(r_threshold = "r-threshold: 0.4")) %>%
  left_join(., biome_type, by = "Type") %>% 
  mutate(Type = ordered(Type, levels = arrange(., Biome) %>% .$Type %>% unique())) %>%
  left_join(., select(new_dataset_names, -New_Type), by = "Type") %>%
  mutate(New_Type = ordered(New_Type, levels = arrange(., Biome) %>% .$New_Type %>% unique()))

write_csv(cluster_verification_combined, "output/Source_Data_Figure_S4.csv")

# Visualize cluster-stability-significance as difference between ARI values in a 
# smoothed scatter-plot
cluster_verification_combined %>%
  ggplot(aes(x = (1-Attack_Strength), y = ARI_difference, col = New_Type)) +
  geom_point() +
  geom_smooth(se = F) +
  scale_color_manual(values = comb_palette) +
  facet_grid(Biome~r_threshold) +
  theme_bw() +
  labs(x = "Attack strength (proportion of removed nodes)", 
       y = "ARI difference\n(Difference between observed and randomized network-cluster)",
       col = "Dataset") +
  guides(colour = guide_legend(ncol = 1))

ggsave("figs/Cluster_Attack_Stability.png", width = 6, height = 10, dpi = 300)

# Test if higher r-threshold yields significantly higher cluster-stability-significance
wilcox.test(filter(cluster_verification_combined, r_threshold == "r-threshold: 0.25") %>% .$ARI_difference,
            filter(cluster_verification_combined, r_threshold == "r-threshold: 0.4") %>% .$ARI_difference)
