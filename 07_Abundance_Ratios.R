## This script computes abundance-ratios of cohorts in the two datasets
# Soil: Australia and Marine: CalCOFI
# For that, count-data is transformed into a flat datatable and
# their log-abundance-ratio is computed for every sample in the dataset.
# Resulting abundance-ratios are plotted against either temperature or
# max. annual temperature based on the WorldClim-data.
#
# Author: Felix Milke
# Date: 26.06.2025

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

# Load information about clusters
load(file = "output/clustering_datasets_025.RData")
load(file = "output/clustering_datasets_04.RData")

# Define max. number of organisms used in prior analyses
max_organisms <- 5000

# Load CalCOFI abundance data
datalist_CalCOFI <- import_data("data/Count_Data_Filtered/CalCOFI/", kingdom = "Prok") %>%
  mutate_meta_datalist(Date = lubridate::mdy(Date))

# Load Soil:Australia abundance data
datalist_Australia <- import_data("data/Count_Data_Filtered/Soil_Australia/", kingdom = "Prok")

# Create new datalist with abundance data clustered according to the co-occurrence clusters
datalist_cluster_CalCOFI <- datalist_CalCOFI
datalist_cluster_Australia <- datalist_Australia

datalist_cluster_CalCOFI$Count_Data <- datalist_CalCOFI %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., read_csv("../Archive/Comparison_Paper/output/SparCC_Cluster/SparCC_Cluster.csv"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster))

datalist_cluster_Australia$Count_Data <- datalist_Australia %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., dplyr::rename(clustering_datasets_0.4[[20]]$Cluster_Tab, "OTU_ID" = "ASV"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster))

# Transform cluster-datalist into flat datatable with abundance- and meta-data combined
datatable_cluster_CalCOFI <- datalist_cluster_CalCOFI %>%
  create_datatable(., grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0)

datatable_cluster_Australia <- datalist_cluster_Australia %>%
  create_datatable(., grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0)

# Compute abundance-ratios for each sample and visualize the results in 
# scatter-plots together with linear models
datatable_cluster_CalCOFI %>%
  filter(Depth <= 20) %>%
  mutate(Season = ordered(Season, levels = c("Spring", "Summer", "Autumn", "Winter"))) %>%
  group_by(Sample_ID) %>%
  mutate(Abundance_ratio = log(Abundance[Group == "2"] / Abundance[Group == "6"], base = 2)) %>%
  ungroup() %>%
  filter(is.finite(Abundance_ratio)) %>%
  ggplot(aes(y = Abundance_ratio, x = T_degC)) +
    geom_point() +
    geom_smooth(method = "lm", se = F, col = "darkred") +
    facet_wrap(~Season) +
    theme_bw() +
    labs(y = "Log. ratio (cohort 2 : cohort 6)", x = "Temperature (°C)")

ggsave("figs/CalCOFI_Abundance_ratio_years.png", width = 6, height = 5, dpi = 300)

datatable_cluster_Australia %>%
  mutate(collection_date = str_replace_all(string = collection_date, pattern = " UTC", replacement = "") %>%
           str_replace_all(pattern = " [0-9]*:[0-9]*:[0-9]*", replacement = "")) %>%
  mutate(Date = lubridate::ymd(collection_date)) %>%
  mutate(Month = lubridate::month(Date)) %>%
  mutate(Season = ifelse(Month %in% c(12, 1, 2), "Summer",
                         ifelse(Month %in% c(3, 4, 5), "Autumn",
                                ifelse(Month %in% c(6, 7, 8), "Winter", "Spring")))) %>%
  group_by(Sample_ID) %>%
  mutate(Abundance_ratio = log(Abundance[Group == "3"] / Abundance[Group == "1"], base = 2)) %>%
  ungroup() %>%
  filter(is.finite(Abundance_ratio)) %>%
  ggplot(aes(y = Abundance_ratio, x = temp_max)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, col = "darkred") +
  facet_wrap(~Season) +
  theme_bw() +
  labs(y = "Log. ratio (cohort 3 : cohort 1)", x = "Max. annual temperature (°C)")

ggsave("figs/Australia_Abundance_ratio_years.png", width = 6, height = 5, dpi = 300)

#### Create Source-Data File ####

datatable_cluster_Australia %>%
  mutate(collection_date = str_replace_all(string = collection_date, pattern = " UTC", replacement = "") %>%
           str_replace_all(pattern = " [0-9]*:[0-9]*:[0-9]*", replacement = "")) %>%
  mutate(Date = lubridate::ymd(collection_date)) %>%
  mutate(Month = lubridate::month(Date)) %>%
  mutate(Season = ifelse(Month %in% c(12, 1, 2), "Summer",
                         ifelse(Month %in% c(3, 4, 5), "Autumn",
                                ifelse(Month %in% c(6, 7, 8), "Winter", "Spring")))) %>%
  group_by(Sample_ID) %>%
  mutate(Abundance_ratio = log(Abundance[Group == "3"] / Abundance[Group == "1"], base = 2)) %>%
  ungroup() %>%
  filter(is.finite(Abundance_ratio)) %>%
  select(Abundance_ratio, Sample_ID, temp_max, Season) %>%
  write_csv("output/Source_Data_Figure_3b.csv")

datatable_cluster_CalCOFI %>%
  filter(Depth <= 20) %>%
  mutate(Season = ordered(Season, levels = c("Spring", "Summer", "Autumn", "Winter"))) %>%
  group_by(Sample_ID) %>%
  mutate(Abundance_ratio = log(Abundance[Group == "2"] / Abundance[Group == "6"], base = 2)) %>%
  ungroup() %>%
  filter(is.finite(Abundance_ratio)) %>%
  select(Abundance_ratio, Sample_ID, T_degC, Season) %>%
  write_csv("output/Source_Data_Figure_3e.csv")
