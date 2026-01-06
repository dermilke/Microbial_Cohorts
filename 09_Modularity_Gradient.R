## This script visualizes environmental preferences of cohorts along the temperature gradient
# together with the  modularity of individual samples along the temperature gradient
# 
# Author: Felix Milke
# Date: 26.06.2025

library(tidyverse)
library(vegan)
library(igraph)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Names_Colors_Improvement.R")

# Load abundance-data network data
datalist_Combined <- import_data("data/Count_Data_Filtered/Comparison/", kingdom = "Prok")

network_Combined <- read_graph("../Archive/Comparison_Paper/output/SparCC_Network.txt", format = "graphml") %>%
  igraph::set_vertex_attr(name = "Cluster", index = read_csv("../Archive/Comparison_Paper/output/SparCC_Cluster/SparCC_Cluster.csv")$OTU_ID,
                          value = read_csv("../Archive/Comparison_Paper/output/SparCC_Cluster/SparCC_Cluster.csv")$Cluster)

datalist_Mendota <- import_data("data/Count_Data_Filtered/Lake_Mendota/", kingdom = "Prok")

datalist_Australia <- import_data("data/Count_Data_Filtered/Soil_Australia/", kingdom = "Prok")

load(file = "output/clustering_datasets_04.RData")

# Create datalists where rows are summarized into cohort abundances
datalist_cluster_Mendota <- datalist_Mendota
datalist_cluster_Combined <- datalist_Combined
datalist_cluster_Australia <- datalist_Australia

datalist_cluster_Mendota$Count_Data <- datalist_cluster_Mendota %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., dplyr::rename(clustering_datasets_0.4[[12]]$Cluster_Tab, "OTU_ID" = "ASV"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster))

datalist_cluster_Combined$Count_Data <- datalist_cluster_Combined %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., read_csv("../Archive/Comparison_Paper/output/SparCC_Cluster/SparCC_Cluster.csv"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster))

datalist_cluster_Australia$Count_Data <- datalist_cluster_Australia %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., dplyr::rename(clustering_datasets_0.4[[20]]$Cluster_Tab, "OTU_ID" = "ASV"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster))

# Transform datalist into flat datatable format
datatable_cluster_Mendota <- datalist_cluster_Mendota %>%
  create_datatable(., grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0) 

datatable_cluster_Combined <- datalist_cluster_Combined %>%
  create_datatable(., grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0) 

datatable_cluster_Australia <- datalist_cluster_Australia %>%
  create_datatable(., grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0) 

# Visualize cohort abundances along temperature gradients
p1_Mendota <- datatable_cluster_Mendota %>%
  filter(!is.na(Temperature)) %>%
  filter(Group != "Other Cluster") %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Temperature)$Sample_ID))) %>%
  group_by(Sample_ID, Group) %>%
  summarize(Abundance = mean(Abundance)) %>%
  ungroup() %>%
  mutate(Group = ordered(Group, levels = seq(1, max(Group)))) %>%
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cbbPalette[1:max(datatable_cluster_Mendota$Group)]) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Rel. Abundance (%)", fill = "Cluster") +
  theme(axis.text.x = element_blank())

p1_Combined <- datatable_cluster_Combined %>%
  filter(!is.na(Pot_Temperature)) %>%
  filter(Group != "Other Cluster") %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Pot_Temperature)$Sample_ID))) %>%
  group_by(Sample_ID, Group) %>%
  summarize(Abundance = mean(Abundance)) %>%
  ungroup() %>%
  mutate(Group = ordered(Group, levels = seq(1, max(Group)))) %>%
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(cbbPalette, "darkgreen", "purple")) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Rel. Abundance (%)", fill = "Cluster") +
  theme(axis.text.x = element_blank())

p1_Australia <- datatable_cluster_Australia %>%
  filter(!is.na(temp_max)) %>%
  filter(Group != "Other Cluster") %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., temp_max)$Sample_ID))) %>%
  group_by(Sample_ID, Group) %>%
  summarize(Abundance = mean(Abundance)) %>%
  ungroup() %>%
  mutate(Group = ordered(Group, levels = seq(1, max(Group)))) %>%
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cbbPalette[1:max(datatable_cluster_Mendota$Group)]) +
  theme_bw() +
  labs(x = "Max. annual temperature (°C)", y = "Rel. Abundance (%)", fill = "Cluster") +
  theme(axis.text.x = element_blank())

## Step 2: Compute local network modularity in each sample

# Compute local network modularity by filtering out all nodes in the network
# that are present in an individual sample and calculate the resulting modularity
# for that particular sample. Iterate over all samples in a dataset.
network_properties_Mendota <- purrr::map_df(datalist_Mendota$Meta_Data$Sample_ID, function(x) {
  OTUs_tmp <- datalist_Mendota$Count_Data$OTU_ID[pull(datalist_Mendota$Count_Data, !!x) > 0]
  
  network_tmp <- induced_subgraph(clustering_datasets_0.4[[12]]$Network, 
                                  vids = names(V(clustering_datasets_0.4[[12]]$Network)) %in% OTUs_tmp)
  modularity_tmp <- modularity(network_tmp, membership = vertex_attr(network_tmp, name = "Cluster"))
  
  return(tibble(Sample_ID = x, Modularity = modularity_tmp))
})

network_properties_Combined <- purrr::map_df(datalist_Combined$Meta_Data$Sample_ID, function(x) {
  OTUs_tmp <- datalist_Combined$Count_Data$OTU_ID[pull(datalist_Combined$Count_Data, !!x) > 0]
  
  network_tmp <- induced_subgraph(network_Combined, 
                                  vids = names(V(network_Combined)) %in% OTUs_tmp)
  modularity_tmp <- modularity(network_tmp, membership = vertex_attr(network_tmp, name = "Cluster"))
  
  return(tibble(Sample_ID = x, Modularity = modularity_tmp))
})

network_properties_Australia <- purrr::map_df(datalist_Australia$Meta_Data$Sample_ID, function(x) {
  OTUs_tmp <- datalist_Australia$Count_Data$OTU_ID[pull(datalist_Australia$Count_Data, !!x) > 0]
  
  network_tmp <- induced_subgraph(clustering_datasets_0.4[[20]]$Network, 
                                  vids = names(V(clustering_datasets_0.4[[20]]$Network)) %in% OTUs_tmp)
  modularity_tmp <- modularity(network_tmp, membership = vertex_attr(network_tmp, name = "Cluster"))
  
  return(tibble(Sample_ID = x, Modularity = modularity_tmp))
})

# Visualize the local modularity along the Temperature gradient
p2_Mendota <- datatable_cluster_Mendota %>%
  left_join(., network_properties_Mendota, by = "Sample_ID") %>%
  ggplot(aes(x = Temperature, y = Modularity)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1, se = F, col = "darkred", lwd = 1.3) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Modularity")

p2_Combined <- datatable_cluster_Combined %>%
  left_join(., network_properties_Combined, by = "Sample_ID") %>%
  ggplot(aes(x = Pot_Temperature, y = Modularity)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1, se = F, col = "darkred", lwd = 1.3) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Modularity")

p2_Australia <- datatable_cluster_Australia %>%
  left_join(., network_properties_Australia, by = "Sample_ID") %>%
  ggplot(aes(x = temp_max, y = Modularity)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1, se = F, col = "darkred", lwd = 1.3) +
  theme_bw() +
  labs(x = "Max. annual temperature (°C)", y = "Modularity")

cowplot::plot_grid(p1_Mendota + theme(legend.position = "none"), 
                   p1_Combined + theme(legend.position = "none"), 
                   p1_Australia + theme(legend.position = "none"), 
                   p2_Mendota, p2_Combined, p2_Australia,
                   ncol = 3, nrow = 2)

ggsave("figs/Ratio_Modularity.png", width = 14, height = 6, dpi = 300)

ggsave(plot = cowplot::get_legend(p1_Combined), filename = "figs/Legend_Modularity_Barplot.png")

### Create Source Data File ###

bind_rows(
  datatable_cluster_Mendota %>%
    left_join(., network_properties_Mendota, by = "Sample_ID") %>%
    select(Temperature, Modularity, Abundance, Group, Sample_ID) %>%
    mutate(Group = as.character(paste0("Cluster ", Group))) %>%
    mutate(Dataset = "Lake Mendota"),
  
  datatable_cluster_Combined %>%
    left_join(., network_properties_Combined, by = "Sample_ID") %>%
    select(Pot_Temperature, Modularity, Abundance, Group, Sample_ID) %>%
    mutate(Group = as.character(paste0("Cluster ", Group))) %>%
    dplyr::rename("Temperature" = "Pot_Temperature") %>%
    mutate(Dataset = "Ocean Latitude Transects"),
  
  datatable_cluster_Australia %>%
    left_join(., network_properties_Australia, by = "Sample_ID") %>%
    select(temp_max, Modularity, Abundance, Group, Sample_ID) %>%
    mutate(Group = as.character(paste0("Cluster ", Group))) %>%
    dplyr::rename("Temperature" = "temp_max") %>%
    mutate(Dataset = "Soil Australia")
) %>%
  filter(!is.na(Temperature)) %>%
  write_csv("output/Source_Data_Figure_4.csv")
