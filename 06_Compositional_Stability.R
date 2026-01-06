## This script computes average cluster-composition and compares them using Bray-Curtis
# To do that, the script loads all datasets and clusters count-data according to the network-clusters
# inferred in script 02_Run_Clustering_Studio.R for both r-thresholds (r = 0.25 & r = 0.4).
# From this clustered abundance-data, the average abundance of each ASV in each dataset and cluster
# is saved to compute the compositional similarity (here Bray-Curtis) between differernt cluster from different datasets.
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

source("R/Names_Colors_Improvement.R")

load("output/clustering_datasets_025.RData")
load("output/clustering_datasets_04.RData")

# Initialising parallelisation. Adjust number of workers according to your CPUs.
plan(multisession, workers = 2)

# This function iterates over all datasets and assigns the cluster-information to all
# ASVs. Then, it computes the mean-abundance for each ASV and saves it in a table
# together with the information of which ASV, which Cluster, which dataset, and which r-threshold.
cluster_composition <- furrr::future_map2(clustering_datasets_0.25, clustering_datasets_0.4, function(cluster_object_0.25, cluster_object_0.4) {
  
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
  
  # Iterate over both r-thresholds and calculate the average abundance of each ASV,
  # then assign cluster-information to each ASV and only keep those ASVs that could be assigned
  # to a cluster.
  avg_composition <- purrr::map(c(0.25, 0.4), function(r_thresh) {
    
    datalist_tmp$Count_Data %>%
      mutate(total = rowSums(select_if(., is.numeric))/ncol(select_if(., is.numeric))) %>%
      select(OTU_ID, total) %>%
      left_join(., dplyr::rename(cluster_tab_combined, "OTU_ID" = "ASV") %>% filter(r_threshold == r_thresh), by = "OTU_ID") %>%
      filter(!is.na(Cluster)) %>%
      with(., reshape2::dcast(data = ., formula = OTU_ID ~ Cluster, value.var = "total")) %>%
      tibble::as_tibble() %>%
      mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
      mutate(r_threshold = r_thresh)
    
  })
  
  return(list(Type = cluster_object_0.25$Type,
              avg_composition = avg_composition))
})

# Save the results into an RDS-file
saveRDS(cluster_composition, "output/Cluster_Composition_Complete.rds")

cluster_composition <- readRDS("output/Cluster_Composition_Complete.rds")

# Iterate over both r-thresholds and calculate bray-curtis similarity matrix for each data-subset
# (here we use Distalgut-samples, Secretion, Skin and Indoor-Surface samples, since we identified
# overlapping clusters in these datasets)
bc_mat_list <- purrr::map(c(0.25, 0.4), function(r_thresh) {
  
  if (r_thresh == 0.25) ind = 1 else ind = 2
  
  # Calculate Bray-Curtis similarity of the dataset
  bc_mat <- purrr::map(cluster_composition, function(comp_tab) {
    tmp <- comp_tab$avg_composition[[ind]] %>%
      mutate(Type = comp_tab$Type) %>%
      left_join(., biome_type, by = "Type") %>%
      filter(Biome %in% c("Distalgut", "Human Secretion", "Human skin", "Indoor surface")) %>%
      filter(Dataset == "EMP") %>%
      magrittr::set_colnames(c("OTU_ID", paste0(seq(1, ncol(.)-(ncol(biome_type)+2)), "_", unique(.$Type)), "r_threshold", names(biome_type))) %>% 
      select(-r_threshold, -Type) %>%
      mutate_if(is.numeric, function(x) x/sum(x))
    
    if (nrow(tmp) == 0) return(NULL) else return(tmp)
  }) %>%
    bind_rows() %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
    group_by(OTU_ID) %>%
    summarize_if(is.numeric, sum) %>%
    select_if(is.numeric) %>%
    t() %>%
    vegan::vegdist()
})
  
# Make annotation-tables for the pheatmap-plot (here: r-threshold = 0.25)
annotation_0.25 <- tibble(Type = colnames(as.matrix(bc_mat_list[[1]])) %>%
        str_remove(pattern = "^[0-9]*_")) %>%
  left_join(., biome_type) %>%
  filter(Biome %in% c("Distalgut", "Human Secretion", "Human skin", "Indoor surface")) %>%
  filter(Dataset == "EMP") %>%
  select(Biome, Host) %>%
  as.data.frame() %>%
  magrittr::set_rownames(colnames(as.matrix(bc_mat_list[[1]])))
  
# Make annotation-tables for the pheatmap-plot (here: r-threshold = 0.4)
annotation_0.4 <- tibble(Type = colnames(as.matrix(bc_mat_list[[2]])) %>%
                            str_remove(pattern = "^[0-9]*_")) %>%
  left_join(., biome_type) %>%
  filter(Biome %in% c("Distalgut", "Human Secretion", "Human skin", "Indoor surface")) %>%
  filter(Dataset == "EMP") %>%
  select(Biome, Host) %>%
  as.data.frame() %>%
  magrittr::set_rownames(colnames(as.matrix(bc_mat_list[[2]])))

# Set colors for the annotation-tables for the pheatmap-plot
annotation_colors = list(Biome = magrittr::set_names(Palette3[-4], unique(annotation_0.4$Biome)),
                         Host = magrittr::set_names(cbbPalette[1:length(unique(annotation_0.4$Host))], unique(annotation_0.4$Host)))

# Use pheatmap to visualize results (r-threshold = 0.4)
pheat <- pheatmap::pheatmap(1-as.matrix(bc_mat_list[[2]]), annotation_col = annotation_0.4, 
                   show_rownames = F, show_colnames = F, border_color = NA,
                   color = c(hcl.colors(n = 10, palette = "Sunset", rev = F), 
                             rep(hcl.colors(n = 10, palette = "Sunset", rev = F)[10], 1)),
                   annotation_row = annotation_0.4, annotation_colors = annotation_colors,
                   treeheight_row = 0, treeheight_col = 0)

ggsave(filename = "figs/Heatmap_Cluster_BC_0.4.png", plot = pheat, width = 6, height = 4.8, dpi = 300)

# Use pheatmap to visualize results (r-threshold = 0.25)
pheat <- pheatmap::pheatmap(1-as.matrix(bc_mat_list[[1]]), annotation_col = annotation_0.25, 
                            show_rownames = F, show_colnames = F, border_color = NA,
                            color = c(hcl.colors(n = 10, palette = "Sunset", rev = F), 
                                      rep(hcl.colors(n = 10, palette = "Sunset", rev = F)[10], 1)),
                            annotation_row = annotation_0.25, annotation_colors = annotation_colors,
                            treeheight_row = 0, treeheight_col = 0)

ggsave(filename = "figs/Heatmap_Cluster_BC_0.25.png", plot = pheat, width = 6, height = 4.8, dpi = 300)

#### Create Source-Data File ####

bind_rows(
  {1-as.matrix(bc_mat_list[[2]])} %>%
    reshape2::melt() %>%
    magrittr::set_colnames(c("From_Cluster", "To_Cluster", "BC_Similarity")) %>%
    mutate(r_threshold = 0.4) %>%
    as_tibble() %>%
    left_join(., as_tibble(annotation_0.4, rownames = "Cluster") %>%
                dplyr::rename_all(function(x) paste0("From_", x)), by = "From_Cluster") %>%
    left_join(., as_tibble(annotation_0.4, rownames = "Cluster") %>%
                dplyr::rename_all(function(x) paste0("To_", x)), by = "To_Cluster"),
  
  {1-as.matrix(bc_mat_list[[1]])} %>%
    reshape2::melt() %>%
    magrittr::set_colnames(c("From_Cluster", "To_Cluster", "BC_Similarity")) %>%
    mutate(r_threshold = 0.25) %>%
    as_tibble() %>%
    left_join(., as_tibble(annotation_0.25, rownames = "Cluster") %>%
                dplyr::rename_all(function(x) paste0("From_", x)), by = "From_Cluster") %>%
    left_join(., as_tibble(annotation_0.25, rownames = "Cluster") %>%
                dplyr::rename_all(function(x) paste0("To_", x)), by = "To_Cluster")
) %>%
  write_csv("output/Source_Data_Figure_2a.csv")
