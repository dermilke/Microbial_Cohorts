## This script calculates environmental preferences of cohorts from the two datasets
# Soil: Australia and Marine: Comparison
# Further, this script calculates the amount of abundance-variance that can be explained by 
# environmental parameter using PERMANOVA (adonis2() function).
#
# The Soil: Australia dataset uses the co-occurrence clustering used in this study,
# whereas the Marine: Comparison dataset uses the previously characterized open ocean
# cohorts from Milke et al. (2023). We specifically used these clusters here since they
# have already been ecologically described and validated. We further used these 
# cohorts in the biogeographic analysis of the CalCOFI dataset (see script "07_Abundance_Ratios.R").
#
# We only analyse here cohorts with at least 10 ASV members, since clusters with too few members
# are likely undersampled. 
#
# First we calculate the environmental preferences for each cohort as their abundance-weighted 
# average of environmental parameter values. Then, the env. preferences are visualized
# using a clustered heatmap.
#
# For the PERMANOVA data we subset the dataset to only include samples where the env.
# parameters were measured. Then, it creates a Bray-Curtis matrix for the cohort-abundance-data
# and then it measures the explained variance of all env. parameter using adonis2() R^2 value.
#
# Author: Felix Milke
# Date: 26.06.2025

library(tidyverse)
library(vegan)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

# This function calculates the abundance-weighted average of all organisms in a datalist using the selected "param" 
# Parameter that is found in the datalist-Meta-Data.
weighted_mean_datalist <- function(datalist, param) {
  param <- enquo(param)
  otu_id <- pull(datalist$Count_Data, 1)
  param_val <- pull(datalist$Meta_Data, !!param)
  
  final <- purrr::map(otu_id, function(x) {
    tmp <- weighted.mean(x = param_val, w = filter(datalist$Count_Data, pull(datalist$Count_Data, 1) == !!x) %>%
                           select_if(is.numeric) %>%
                           as.numeric(), na.rm = T)
    return(tibble(Cluster = rlang::as_name(x), !!rlang::quo_name(param) := tmp))
  }) %>%
    bind_rows()
  
  return(final)
}

# Load cluster data for the Soil: Australia dataset
load("output/clustering_datasets_04.RData")

# Load abundance data. The Marine: Comparison dataset consists of two latitudinal transects,
# one conducted in the Atlantic Ocean and one in the Pacific Ocean (see Milke et al. (2023)).
# For our analyses here we combine both datasets into a single datalist.
datalist_Atlantic <- import_data("../Archive/Comparison_Paper/data/Atlantic/", kingdom = "Prok") %>%
  mutate_meta_datalist(Station = as.character(Station))
datalist_Pacific <- import_data("../Archive/Comparison_Paper/data/Pacific/", kingdom = "Prok")

datalist_Combined <- combine_data(datalist_Atlantic, datalist_Pacific)

datalist_Australia <- import_data("../Australia_Soils/data/Count_Data_Filtered/", kingdom = "Prok")

# Add WordlClim data to the Soil: Australia Meta-Data
datalist_Australia$Meta_Data <- left_join(datalist_Australia$Meta_Data, 
                                          read_tsv("../Australia_Soils/data/worldclim_data/bioclimatic.tsv"), by = "Sample_ID")

# Filter out those clusters with fewer than 10 members
cluster_to_choose_Australia <- clustering_datasets_0.4[[14]]$Cluster_Tab %>%
  group_by(Cluster) %>%
  summarize(N = n()) %>%
  filter(N > 10) %>%
  .$Cluster

cluster_to_choose_Combined <- read_csv("../Archive/Comparison_Paper/output/SparCC_Cluster/SparCC_Cluster.csv") %>%
  group_by(Cluster) %>%
  summarize(N = n()) %>%
  filter(N > 10) %>%
  .$Cluster

# Create new datalists where count-data is summarized into cluster-data
datalist_cluster_Australia <- datalist_Australia
datalist_cluster_Combined <- datalist_Combined

datalist_cluster_Australia$Count_Data <- datalist_cluster_Australia %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., dplyr::rename(clustering_datasets_0.4[[14]]$Cluster_Tab, "OTU_ID" = "ASV"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  slice(cluster_to_choose_Australia)

datalist_cluster_Combined$Count_Data <- datalist_cluster_Combined %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., read_csv("../Archive/Comparison_Paper/output/SparCC_Cluster/SparCC_Cluster.csv"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  slice(cluster_to_choose_Combined)

# Selected environmental parameter
params_Australia <- c("ammonium_nitrogen_wt", "clay", "fine_sand",  "sand", "silt",
                      "exc_calcium", "exc_potassium", "nitrate_nitrogen", "organic_carbon", "ph", "phosphorus_colwell", 
                      "sulphur", "temp_annual", "temp_range", "precipitation_annual", "precip_seasonality")

params_Combined <- c("Si", "NO3", "Depth",  "Bacterial_Generation_Time", "Fluorescence",
                      "Pot_Temperature", "Bacteria_Abundance_Flow", "Salinity", "Oxygen")

# Compute env. preferences using the weighted_mean_datalist() function from above.
env_preference_cluster_Australia <- purrr::map(params_Australia, function(x) {
  weighted_mean_datalist(datalist_cluster_Australia, x) %>%
    dplyr::rename(!!x := "x") %>%
    select(-Cluster)
}) %>%
  bind_cols() %>%
  mutate(Cluster = datalist_cluster_Australia$Count_Data$Cluster, .before = 1)

env_preference_cluster_Combined <- purrr::map(params_Combined, function(x) {
  weighted_mean_datalist(datalist_cluster_Combined, x) %>%
    dplyr::rename(!!x := "x") %>%
    select(-Cluster)
}) %>%
  bind_cols() %>%
  mutate(Cluster = datalist_cluster_Combined$Count_Data$Cluster, .before = 1)

# Visualize env. preferences. For that, we first normalize the env. preferences using
# TSS (total sum scaling) and then standardize the values using z-scoring.
p1 <- env_preference_cluster_Australia %>%
  select_if(is.numeric) %>%
  mutate_all(function(x) x/sum(x)) %>%
  with(., apply(., 1, function(x) (x-mean(x))/sd(x))) %>%
  magrittr::set_rownames(c("Ammonium-nitrogen", "clay", "Fine sand", "Sand", "Silt",
                           "Exc. calcium", "Exc. potassium", "Nitrate-nitrogen", "Organic carbon", 
                           "PH", "Phosphorous", "Sulphur", "Annual avg. temperature",
                           "Temperature range", "Annual avg. precipitation", "Precipitation seasonality")) %>%
  pheatmap::pheatmap(labels_col = as.character(1:nrow(datalist_cluster_Australia$Count_Data)))

ggsave(p1, filename = "fig/Cohort_Env_Preferences_Australia.png", width = 7, height = 5.5, dpi = 300)

p2 <- env_preference_cluster_Combined %>%
  select_if(is.numeric) %>%
  mutate_all(function(x) x/sum(x)) %>%
  with(., apply(., 1, function(x) (x-mean(x))/sd(x))) %>%
  magrittr::set_rownames(c("Silicate", "Nitrate", "Depth", "Prokaryotic Generation Time", "Fluorescence",
                           "Potential Temperature", "Prokaryotic Cell Numbers", "Salinity", "Oxygen")) %>%
  pheatmap::pheatmap(labels_col = as.character(1:nrow(datalist_cluster_Combined$Count_Data)))

ggsave(p2, filename = "fig/Cohort_Env_Preferences_Combined.png", width = 7, height = 5.5, dpi = 300)

# Subset each dataset to only include samples where the considered environmental parameter
# are not NA
select_station_Australia <- datalist_cluster_Australia$Meta_Data %>%
  select(params_Australia) %>%
  mutate(index = seq(1, nrow(.))) %>%
  filter(!is.na(rowSums(.))) %>%
  mutate(Sample_ID = datalist_cluster_Australia$Meta_Data$Sample_ID[index], .before = 1) %>%
  .$Sample_ID

select_station_Combined <- datalist_cluster_Combined$Meta_Data %>%
  select(params_Combined) %>%
  mutate(index = seq(1, nrow(.))) %>%
  filter(!is.na(rowSums(.))) %>%
  mutate(Sample_ID = datalist_cluster_Combined$Meta_Data$Sample_ID[index], .before = 1) %>%
  .$Sample_ID

# Calculate Bray-Curtis matrices for the data-subsets
bc_count_Australia <- datalist_cluster_Australia %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_station_datalist(Sample_ID %in% select_station_Australia) %>%
  select_meta_datalist(params_Australia) %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  t() %>% vegan::vegdist()

bc_count_Combined <- datalist_cluster_Combined %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_station_datalist(Sample_ID %in% select_station_Combined) %>%
  select_meta_datalist(params_Combined) %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  t() %>% vegan::vegdist()

# Run the adonis2() function using all considered env. parameter with standard-parameters
model_Australia <- vegan::adonis2(as.formula(paste0("bc_count_Australia ~ ", paste0(params_Australia, collapse = " + "))), 
                                  data = filter_station_datalist(datalist_cluster_Australia, 
                                                                 Sample_ID %in% select_station_Australia)$Meta_Data)

model_Comparison <- vegan::adonis2(as.formula(paste0("bc_count_Combined ~ ", paste0(params_Combined, collapse = " + "))), 
                                   data = filter_station_datalist(datalist_cluster_Combined, 
                                                                  Sample_ID %in% select_station_Combined)$Meta_Data)
