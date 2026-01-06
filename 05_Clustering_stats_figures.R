## This script computes cluster properties for all datasets under two different r-thresholds
# To do that, the script loads all datasets and clusters count-data according to the network-clusters
# inferred in script 02_Run_Clustering_Studio.R for both r-thresholds (r = 0.25 & r = 0.4).
# Then, based on the cluster-information, the following metrics are inferred:
# - ASV per cohort
# - Richness of cohorts
# - Effective Number of Species (Cohorts)
# - Proportion of total counts associated to cohorts
#
# Afterwards, all metrics are combined and visualized.
#
# The script uses parallelization, adjust to your likening.
#
# Author: Felix Milke
# Date: 26.06.2025

library(tidyverse)
library(furrr)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

# Here I store tables for plotting / improving dataset-names for visualization
source("R/Names_Colors_Improvement.R")

load("output/clustering_datasets_025.RData")
load("output/clustering_datasets_04.RData")

# Initialising parallelisation. Adjust number of workers according to your CPUs.
plan(multisession, workers = 2)

# Iterate over all datasets, calculate the chosen metrics for both r-thresholds
cluster_stats <- furrr::future_map2(clustering_datasets_0.25, clustering_datasets_0.4, function(cluster_object_0.25, cluster_object_0.4) {
  
  cat(cluster_object_0.25$Type)
  
  datalist_tmp <- import_data(paste0("data/Count_Data_Filtered/", cluster_object_0.25$Type, "/"), kingdom = "Prok", min_counts = 2000) %>%
    mutate_count_datalist(function(x) x/sum(x))
  
  # Combine cluster-information of both r-thresholds
  cluster_tab_combined <- bind_rows(
    cluster_object_0.25$Cluster_Tab %>%
      mutate(r_threshold = 0.25),
    cluster_object_0.4$Cluster_Tab %>%
      mutate(r_threshold = 0.4)
  )
  
  # Iterate over both r-thresholds:
  # - Filter out all ASVs from the chosen-ASVs that were assigned to a cluster/cohort
  # - Sum up their total rel. abundance of a sample
  # -> Get the metric: "Proportion of total counts associated to cohorts"
  mapped_counts <- purrr::map(c(0.25, 0.4), function(r_thresh) {
    datalist_tmp$Count_Data %>%
    left_join(., dplyr::rename(cluster_tab_combined, "OTU_ID" = "ASV") %>% filter(r_threshold == r_thresh), by = "OTU_ID") %>%
    filter(!is.na(Cluster)) %>%
    select(-Cluster, -r_threshold) %>%
    select_if(is.numeric) %>%
    colSums() %>%
    tibble::enframe(name = "Sample_ID") %>%
    dplyr::rename("Mapped_Counts" = "value") %>%
    mutate(r_threshold = r_thresh)
  }) %>%
    bind_rows()
  
  # Identify number of ASVs associate to a metric for both r-thresholds
  asv_per_cohort <- cluster_tab_combined %>%
    group_by(Cluster, r_threshold) %>%
    summarize(N = n()) %>%
    ungroup()
  
  # Iterate over both r-thresholds:
  # - Assign to each ASV in the filtered dataset their cluster-information
  # - Summmarize count-data to get the total rel. abundance for each Cluster in each sample 
  # - Use this cluster-abundance-data to infer richness, Evenness, Shannon and ESN using 
  #   the "diversity_datatable()" function from the ExCom function-suite.
  cluster_stats <- purrr::map(c(0.25, 0.4), function(r_thresh) {
    
    datalist_tmp_cluster <- datalist_tmp
    
    datalist_tmp_cluster$Count_Data <- datalist_tmp_cluster$Count_Data %>%
      left_join(., dplyr::rename(cluster_tab_combined, "OTU_ID" = "ASV") %>% filter(r_threshold == r_thresh), by = "OTU_ID") %>%
      filter(!is.na(Cluster)) %>%
      mutate(Cluster = as.character(Cluster)) %>%
      select(-r_threshold) %>%
      group_by(Cluster) %>%
      summarize_if(is.numeric, sum)
    
    datalist_tmp_cluster %>%
      diversity_datatable() %>%
      left_join(., filter(mapped_counts, r_threshold == r_thresh), by = "Sample_ID") %>%
      select(Sample_ID, Shannon, Richness, Evenness, ESN, Mapped_Counts, r_threshold) %>%
      mutate(Type = cluster_object_0.25$Type) 
  }) %>%
    bind_rows()

  # Total number of samples in the dataset
  sample_num <- datalist_tmp$Meta_Data %>% nrow()
  
  return(list(ASV_per_cohort = asv_per_cohort,
              cluster_stats = cluster_stats,
              sample_num = sample_num,
              Type = cluster_object_0.25$Type))
})

# Improve names of Biome-Types for plotting
cluster_stats_tab_filtered <- cluster_stats %>%
  purrr::map("cluster_stats") %>%
  bind_rows() %>%
  mutate(Biome = left_join(., biome_type, by = "Type")$Biome)

# Add number of samples to the names of Biome-Types
biome_type_with_num <- cluster_stats_tab_filtered %>%
  group_by(Biome) %>%
  summarize(N = n()) %>%
  left_join(., biome_type, by = "Biome") %>%
  mutate(Biome = paste0(Biome, " (N = ", N/2, ")")) %>%
  select(-N) %>%
  right_join(., select(cluster_stats_tab_filtered, -Biome), by = "Type")

# Assign an order to the different Biome-Types based on their total number of clusters/cohorts
Biome_Order <- biome_type_with_num %>%
  group_by(Biome) %>%
  summarize(Richness = median(Richness)) %>%
  arrange(desc(Richness)) %>%
  distinct() %>%
  .$Biome

# Plot the number of cohorts per sample metric
biome_type_with_num %>%
  mutate(Biome = ordered(Biome, levels = Biome_Order)) %>%
  ggplot(aes(x = Biome, y = Richness, col = as.factor(r_threshold), fill = as.factor(r_threshold))) +
  geom_boxplot(col = "black", width = .2, position = position_dodge(width = .3), linewidth = 0.25, outlier.size = .25) +
  geom_jitter(alpha = .05, position = position_jitterdodge(dodge.width = 1.2, jitter.width = .3,
                                                           jitter.height = .4)) +
  coord_flip() +
  labs(x = "", y = "Number of cohorts per sample\n ", col = "r-threshold", fill = "r-threshold") +
  scale_color_manual(values = cbbPalette[c(7, 4)]) +
  scale_fill_manual(values = cbbPalette[c(7, 4)]) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))

ggsave("figs/Sample_Cohort_Richness_jittered2.png", width = 5, height = 3.5, dpi = 300)

# Plot the effective number of cohorts metric
biome_type_with_num %>%
  mutate(Biome = ordered(Biome, levels = Biome_Order)) %>%
  ggplot(aes(x = Biome, y = ESN, fill = as.factor(r_threshold))) +
    geom_violin(scale = "width", linewidth = 0.25) +
    coord_flip() +
    scale_fill_manual(values = cbbPalette[c(7, 4)]) +
    labs(x = "", y = "Effective Number of\ncohorts per sample", fill = "r-threshold") +
    theme_bw() +
  theme(axis.text = element_text(colour = "black"))

ggsave("figs/Sample_Cohort_ENS.png", width = 5, height = 3.5, dpi = 300)

# Plot the proportion of counts associated to a cohort-metric
biome_type_with_num %>%
  mutate(Biome = ordered(Biome, levels = Biome_Order)) %>%
  ggplot(aes(x = Biome, y = Mapped_Counts*100, fill = as.factor(r_threshold))) +
  geom_violin(scale = "width", linewidth = 0.25) +
  scale_fill_manual(values = cbbPalette[c(7, 4)]) +
  coord_flip() +
  labs(x = "", y = "Mapped counts (%)\n ", fill = "r-threshold") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))

ggsave("figs/Sample_Cohort_Mapped.png", width = 5, height = 3.5, dpi = 300)

# Plot the Number of ASVs associated to a cohort metric
cluster_stats %>%
  purrr::map(function(x) {
    mutate(x$ASV_per_cohort, Type = unique(x$cluster_stats$Type)) %>%
      return()
    }) %>%
  bind_rows() %>%
  left_join(., select(biome_type_with_num, Type, Biome) %>% distinct(), by = "Type") %>%
  filter(Type != "HMP") %>%
  filter(!(Type %in% c("Sediment_Saline", "Sediment_Nonsaline"))) %>%
  mutate(Biome = ordered(Biome, levels = Biome_Order)) %>%
  ggplot(aes(x = Biome, y = N, fill = as.factor(r_threshold))) +
    geom_boxplot(linewidth = 0.25, outlier.size = .25) +
    coord_flip() +
    scale_fill_manual(values = cbbPalette[c(7, 4)]) +
    scale_y_log10() +
    labs(x = "", y = "Number of ASVs per cohort\n ", fill = "r-threshold") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"))

ggsave("figs/ASVs_per_Cluster.png", width = 5, height = 3.5, dpi = 300)

#### Create Source-Data File ####

biome_type_with_num %>%
  select(Biome, ESN, Mapped_Counts, Richness, r_threshold) %>%
  write_csv("output/Source_Data_Figure_1a_1.csv")

cluster_stats %>%
  purrr::map(function(x) {
    mutate(x$ASV_per_cohort, Type = unique(x$cluster_stats$Type)) %>%
      return()
  }) %>%
  bind_rows() %>%
  left_join(., select(biome_type_with_num, Type, Biome) %>% distinct(), by = "Type") %>%
  filter(Type != "HMP") %>%
  filter(!(Type %in% c("Sediment_Saline", "Sediment_Nonsaline"))) %>%
  select(Biome, N, r_threshold) %>%
  dplyr::rename("ASV_per_Cohort" = "N") %>%
  write_csv("output/Source_Data_Figure_1a_2.csv")
