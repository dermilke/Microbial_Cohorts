library(tidyverse)

filter_wrapper <- function(datalist, prop_samples = 0.1, min_organisms) {
  
  set.seed(123)
  
  datalist_singleton <- datalist %>%
    singleton_filter(min_count = 0, min_station = {.$Count_Data %>% select_if(is.numeric) %>% ncol()} * prop_samples)
  
  if (nrow(datalist_singleton$Count_Data) < min_organisms) {
    
    chosen_asvs <- datalist %>%
      mutate_count_datalist(function(x) x/sum(x)) %>%
      .$Count_Data %>%
      mutate(total = rowMeans(select_if(., is.numeric))) %>%
      select(OTU_ID, total) %>%
      arrange(desc(total)) %>%
      dplyr::slice(1:min_organisms) %>%
      .$OTU_ID
    
    datalist_filtered <- datalist %>%
      filter_taxa_datalist(OTU_ID %in% chosen_asvs)
  } else {
    datalist_filtered <- datalist_singleton
  }
  
  return(datalist_filtered)
}

# Define maximum number of ASVs for FastSpar calculations
min_organisms <- 5000

load("output/clustering_datasets_025.RData")
load("output/clustering_datasets_04.RData")

count_data_files <- list.files("data/Count_Data_/", pattern = "Full_Prok_Count.tsv", recursive = T, full.names = T) %>%
  str_replace_all(pattern = "\\/\\/", replacement = "\\/")

types_vec <- purrr::map_chr(clustering_datasets_0.25, "Type")

info_table <- purrr::map(count_data_files, function(count_data_tmp) {
  
  get_infos <- function(datalist_input) {
    tmp <- datalist_input$Count_Data %>%
      select_if(is.numeric)
    
    tibble(Samples = ncol(tmp), ASVs = nrow(tmp), Proportion_avg = mean(colSums(tmp))) %>%
      return()
  }
  
  type_tmp <- str_replace_all(count_data_tmp, pattern = "\\/Processed.*$", replacement = "") %>%
    str_replace_all(pattern = "^data/Count_Data_/", replacement = "")
  
  cluster_data_tmp_025 <- clustering_datasets_0.25[[which(types_vec == type_tmp)]]
  cluster_data_tmp_04 <- clustering_datasets_0.4[[which(types_vec == type_tmp)]]
  
  datalist_tmp <- import_data(paste0("data/Count_Data_/", type_tmp, "/"), kingdom = "Prok") %>%
    mutate_count_datalist(function(x) x/sum(x))
  
  datalist_tmp_filtered <- datalist_tmp %>%
    singleton_filter(min_count = 0, min_station = {.$Count_Data %>% select_if(is.numeric) %>% ncol()} * 0.1)
  
  datalist_tmp_inNetwork <- datalist_tmp %>%
    filter_wrapper(min_organisms = min_organisms)
  
  bind_cols(get_infos(datalist_tmp) %>% dplyr::rename_all(function(x) paste0("Complete_", x)),
            get_infos(datalist_tmp_filtered)[,-1] %>% dplyr::rename_all(function(x) paste0("Filtered_10_", x)),
            get_infos(datalist_tmp_inNetwork)[,-1] %>% dplyr::rename_all(function(x) paste0("Final_Filter", x))) %>%
    mutate(Type = type_tmp, .before = 1) %>%
    mutate(InNetwork_025_ASVs = igraph::vcount(cluster_data_tmp_025$Network)) %>%
    mutate(InNetwork_025_Edges = igraph::ecount(cluster_data_tmp_025$Network)) %>%
    mutate(InNetwork_025_Proportion_Avg = datalist_tmp %>%
             filter_taxa_datalist(OTU_ID %in% cluster_data_tmp_025$Cluster_Tab$ASV) %>%
             .$Count_Data %>% select_if(is.numeric) %>% colSums() %>% mean()) %>%
    mutate(InNetwork_025_Proportion_Avg_FinalFilter = datalist_tmp_inNetwork %>%
             mutate_count_datalist(function(x) x/sum(x)) %>%
             filter_taxa_datalist(OTU_ID %in% cluster_data_tmp_025$Cluster_Tab$ASV) %>%
             .$Count_Data %>% select_if(is.numeric) %>% colSums() %>% mean()) %>%
    mutate(InNetwork_04_ASVs = igraph::vcount(cluster_data_tmp_04$Network)) %>%
    mutate(InNetwork_04_Edges = igraph::ecount(cluster_data_tmp_04$Network)) %>%
    mutate(InNetwork_04_Proportion_Avg = datalist_tmp %>%
             filter_taxa_datalist(OTU_ID %in% cluster_data_tmp_04$Cluster_Tab$ASV) %>%
             .$Count_Data %>% select_if(is.numeric) %>% colSums() %>% mean()) %>%
    mutate(InNetwork_04_Proportion_Avg_FinalFilter = datalist_tmp_inNetwork %>%
             mutate_count_datalist(function(x) x/sum(x)) %>%
             filter_taxa_datalist(OTU_ID %in% cluster_data_tmp_04$Cluster_Tab$ASV) %>%
             .$Count_Data %>% select_if(is.numeric) %>% colSums() %>% mean()) %>%
    return()
  
}) %>%
  bind_rows()

write_csv(info_table, "output/Complete_Info_Table.csv")
