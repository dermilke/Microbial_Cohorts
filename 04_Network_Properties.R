## This script computes network properties for all datasets under varying r-thresholds
# To do that, the script loads all FastSpar correlation matrices for each dataset
# and then iterates over an interval of correlation-thresholds (0.1 to 0.7 with total 30 steps between)
# For each correlation threshold a network is built from the FastSpar correlation matrix.
# From the network four metrics are inferred: 
# a) Modularity - based on cluster_walktrap() memberships
# b) Clustering coefficient - using the transitivity() function
# c) avg. number of ASVs per cluster - based on the cluster_walktrap() memberships
# d) total clusters - based on the cluster_walktrap() memberships
#
# Afterwards, all metrics are combined and visualized.
#
# The script uses parallelization, adjust to your likening.
#
# Author: Felix Milke
# Date: 26.02.2025

library(tidyverse)
library(igraph)
library(furrr)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Import_SparCC_Network.R")
source("R/Names_Colors_Improvement.R")

# Number of workers should be set to your preferred number of threads.
plan(multisession, workers = 2)

# Identify all datasets
cor_mat <- list.files("output/10_percent_threshold", pattern = "cor", full.names = T, recursive = T)

types <- str_extract(cor_mat, pattern = "cor_.*_table") %>%
  str_replace_all(pattern = "cor_", replacement = "") %>%
  str_replace_all(pattern = "_table", replacement = "")

# Iterate over each dataset
network_props_biomes <- purrr::map(cor_mat, function(cor_mat) {
  
  # extract dataset-Type from file-name
  dataset_type <- str_extract(cor_mat, pattern = "cor_.*_table") %>%
    str_replace_all(pattern = "cor_", replacement = "") %>%
    str_replace_all(pattern = "_table", replacement = "")
  
  cat(paste0(dataset_type, "\n"))
  
  # load correlation- and pvalue-matrices and create an edgelist from them.
  # Don't filter for r-thresholds (except r >= 0)
  network_data <- import_sparcc_network(cor_file = cor_mat, 
                          pval_file = str_replace_all(cor_mat, pattern = "cor", replacement = "pvalues"),
                          min_r = 0, min_p = 0.05) %>%
    filter(Type == "Positive")
  
  # Iterate over each r-threshold in the tested interval
  network_props <- purrr::map(seq(0.1, 0.7, length.out = 30), function(r_thresh) {
      
      # filter edgelist according to selected r-threshold and create a network from it
      network_tmp <- network_data %>%
        filter(weight >= r_thresh) %>%
        igraph::graph_from_data_frame(d = ., directed = F) %>%
        igraph::set_vertex_attr(graph = ., name = "label", value = NA) %>%
        igraph::simplify()
      
      # infer all four network properties from the network
      cluster_tmp <- cluster_walktrap(network_tmp)
      mod <- modularity(network_tmp, cluster_tmp$membership)
      cluster_coef <- transitivity(network_tmp)
      
      # consider only clusters that have at least 11 members, all others are discarded since
      # they are likely spurious and not relevant here.
      cluster_table_tmp <-  table(cluster_tmp$membership) %>%
        tibble::enframe() %>%
        filter(value > 10)
      
      avg_members <- mean(cluster_table_tmp$value)
      total_cluster <- length(cluster_table_tmp$value)
      
      return(tibble(r_val = r_thresh, modularity = mod, cluster_coef = cluster_coef,
                    avg_members = avg_members, total_cluster = total_cluster))
      
    }) %>%
      bind_rows() %>%
      mutate(Type = dataset_type)
  
  
  return(list(network_data = network_data, network_coef = network_props))
  
})

# Combine all correlation-values for all datasets
network_data_combined <- purrr::map(seq(1, length(network_props_biomes)), function(x) {
  network_props_biomes[[x]]$network_data %>%
    mutate(Type = network_props_biomes[[x]]$network_coef$Type[1])
}) %>%
  bind_rows() %>%
  mutate(Type = str_replace_all(Type, pattern = "Freshwater_", replacement = "Water_Nonsaline_") %>%
           str_replace_all(pattern = "Indoor_Surface", replacement = "Surface_Nonsaline") %>%
           str_replace_all(pattern = "Skin_Foot", replacement = "Surface_Foot") %>%
           str_replace_all(pattern = "Plant_Root", replacement = "Plant_Roots") %>%
           str_replace_all(pattern = "Secretion_Nose", replacement = "Nose_Secretion") %>%
           str_replace_all(pattern = "Secretion_Saliva", replacement = "Human_Saliva") %>%
           str_replace_all(pattern = "Skin_Hand", replacement = "Surface_Hand")) %>%
  left_join(., biome_type, by = "Type")

# Combine all network-properties for every dataset and every r-threshold
network_props_combined <- purrr::map(network_props_biomes, "network_coef") %>%
  bind_rows() %>%
  mutate(Type = str_replace_all(Type, pattern = "Freshwater_", replacement = "Water_Nonsaline_") %>%
           str_replace_all(pattern = "Indoor_Surface", replacement = "Surface_Nonsaline") %>%
           str_replace_all(pattern = "Skin_Foot", replacement = "Surface_Foot") %>%
           str_replace_all(pattern = "Plant_Root", replacement = "Plant_Roots") %>%
           str_replace_all(pattern = "Secretion_Nose", replacement = "Nose_Secretion") %>%
           str_replace_all(pattern = "Secretion_Saliva", replacement = "Human_Saliva") %>%
           str_replace_all(pattern = "Skin_Hand", replacement = "Surface_Hand")) %>%
  left_join(., biome_type, by = "Type") %>%
  filter(avg_members > 10) %>%
  mutate(New_Type = ordered(New_Type, levels = arrange(., Biome) %>% .$New_Type %>% unique()))

write_csv(network_props_combined, "output/Source_Data_Figure_2b.csv")

# Visualize modularity scores
network_props_combined %>%
  
  ggplot(aes(x = r_val, y = modularity, col = New_Type)) +
    geom_point() +
    geom_smooth(se = T) +
    scale_color_manual(values = comb_palette) +
    facet_wrap(~Biome, ncol = 1) +
    theme_bw() +
    labs(y = "Network modularity", x = "SparCC r-value") +
    theme(legend.position = "none")

ggsave("figs/corco_modularity_fixy_soil_marine.png", width = 3, height = 14, dpi = 300)

# Visualize clustering coefficient scores
network_props_combined %>%
  ggplot(aes(x = r_val, y = cluster_coef, col = New_Type)) +
    geom_point() +
    geom_smooth(se = F) +
    scale_color_manual(values = comb_palette) +
    facet_wrap(~Biome, ncol = 1)  + 
    theme_bw() +
    labs(y = "Network clustering coefficient", x = "SparCC r-value") +
    theme(legend.position = "none")

ggsave("figs/corco_clustercoef_fixy_soil_marine.png", width = 3, height = 14, dpi = 300)

# Visualize avg. number of cluster members
network_props_combined %>%
  ggplot(aes(x = r_val, y = avg_members, col = New_Type)) +
    geom_point() +
    geom_smooth(se = F) +
    scale_y_log10() +
    scale_color_manual(values = comb_palette) +
    facet_wrap(~Biome, ncol = 1) +
    theme_bw() +
    labs(y = "Avg. members of Cl", x = "SparCC r-value") +
    theme(legend.position = "none")

ggsave("figs/corco_avgmember_fixy_soil_marine.png", width = 3, height = 14, dpi = 300)

# Visualize total number of clusters 
# (drop one instance where number of clusters was extremely high and which skewed y-axis scale)
network_props_combined %>%
  filter(total_cluster < 20) %>%
  ggplot(aes(x = r_val, y = total_cluster, col = New_Type)) +
    geom_point() +
    geom_smooth(se = F) +
    scale_color_manual(values = comb_palette) +
    facet_wrap(~Biome, ncol = 1) + 
    theme_bw() +
    labs(y = "Number of cluster", x = "SparCC r-value", col = "Dataset") +
    guides(col = guide_legend(ncol = 1))

ggsave("figs/corco_numberofcluster_fixy_soil_marine.png", width = 4.8, height = 14, dpi = 300)

# Visualize distribution density of correlation-values for each dataset
network_data_combined %>%
  ggplot(aes(x = weight, fill = New_Type)) +
    geom_density(alpha = 0.6) +
    scale_x_log10() +
    scale_fill_manual(values = comb_palette) +
    facet_wrap(~Biome, ncol = 1) + 
    theme_bw() +
    theme(legend.position = "none") +
    labs(y = "Density distribution", x = "SparCC r-values")

ggsave("figs/corco_distribution_fixy_soil_marine.png", width = 3, height = 14, dpi = 300)

