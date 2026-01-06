#### Validation Marine: Comparison - Systematic Ocean Basin Subsets ####

## This script runs the validation of clustering results using data subsets and three different network inference
## algorithms (SparCC, SpiecEasi, Pearson)
# First create systematic data subsets based on individual ocean basins (Atlantic vs. Pacific Ocean), 
# and find all shared ASVs (Comparison between Validation runs is more reliable and fair when 
# considering only shared ASVs).
# Then compute:
# 1. SparCC
# 2. Pearson correlation using CLR transformation
# 3. SpiecEasi
# Since SpiecEasi does not compute correlation weights for the network, we use a network
# clustering algorithm that does not consider weights (cluster_fast_greedy()).
# The other networks are clustered using the in-house clustering pipeline.
# Next, we use five different comparison analyses:
# 1. Cluster Composition: The ASV composition of the resulting network cluster are compared using Bray-Curtis 
# Similarities (Composition is based on average abundance of ASVs).
# 2. Biogeography: We visualize the abundance of clusters along major environmental axes (here max. Temperature)
# using stacked barplots. To identify which cluster belongs to which color, we cluster
# the clusters based on compositional similarity with hierarchical clustering.
# For that, we use their shilouette-width to define the optimal cluster number.
# 3. Environmental Prediction: We compare prediction performance of the resulting cluster
# abundance-ratios to predict the major environmental axes (here max. Temperature).
# 4. Mantel Test Comparison: Compute matrix-correlation between cluster abundance tables.
# 5. Gold Standard Comparison: Finally, we compare the composition of each sample against
# its "gold standard" (defined as the SparCC results using the complete dataset).
# For that, we use Bray-Curtis similarity between each sample and its respective Gold Standard sample.

# This script uses parallelization. Adjust the number of threads manually.
# Especially the FastSpar implementation takes a lot of time.
#
# Author: Felix Milke
# Date: 28.10.2025

library(tidyverse)
library(SpiecEasi)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Fastspar_new.R")
source("R/Cluster_Pipeline.R")
source("R/Import_SparCC_Network.R")

source("R/Names_Colors_Improvement.R")

#### Construct Data Subsets ####

datalist_Combined <- import_data("../Archive/Comparison_Paper/data/Combined/", kingdom = "Prok") %>%
  filter_station_datalist(Size_Fraction == "0.22")

datalist_Atlantic <- datalist_Combined %>%
  filter_station_datalist(Ocean == "Atlantic")

datalist_Pacific <- datalist_Combined %>%
  filter_station_datalist(Ocean == "Pacific")

min_organisms = 5000

chosen_asvs <- datalist_Combined %>%
  filter_taxa_datalist(OTU_ID %in% datalist_Atlantic$Count_Data$OTU_ID & OTU_ID %in% datalist_Pacific$Count_Data$OTU_ID) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  mutate(total = rowMeans(select_if(., is.numeric))) %>%
  select(OTU_ID, total) %>%
  arrange(desc(total)) %>%
  dplyr::slice(1:min_organisms) %>%
  .$OTU_ID

datalist_Atlantic$Count_Data <- datalist_Atlantic$Count_Data %>%
  slice(match(chosen_asvs, OTU_ID))

datalist_Pacific$Count_Data <- datalist_Pacific$Count_Data %>%
  slice(match(chosen_asvs, OTU_ID))

datalist_Combined$Count_Data <- datalist_Combined$Count_Data %>%
  slice(match(chosen_asvs, OTU_ID))

# -> Datalist Atlantic and Pacific only shared 3126 ASVs -> less than 5000!
# Median total abundance per sample covered: 95.8% in Atlantic Ocean and 91% in Pacific Ocean

#### Run SparCC Network Inference ####

## SparCC Data Preparation
datalist_Combined %>%
  .$Count_Data %>%
  select(-c(2:8)) %>%
  rename(`#OTU ID` = "OTU_ID") %>%
  write_tsv("output/Verification/Combined_FL_table.tsv")

datalist_Atlantic %>%
  .$Count_Data %>%
  select(-c(2:8)) %>%
  rename(`#OTU ID` = "OTU_ID") %>%
  write_tsv("output/Verification/Atlantic_FL_table.tsv")

datalist_Pacific %>%
  .$Count_Data %>%
  select(-c(2:8)) %>%
  rename(`#OTU ID` = "OTU_ID") %>%
  write_tsv("output/Verification/Pacific_FL_table.tsv")

## SparCC inference
message("Running FastSpar. This may take some time (hours)...")

FastSparCC_Server_Function(otu_table_input = "output/verification/Combined_FL_table.tsv", 
                           output_folder = "output/verification", 
                           fastSpar_iterations = 20, bootstrap_num = 200, threads = 4)

FastSparCC_Server_Function(otu_table_input = "output/verification/Atlantic_FL_table.tsv", 
                           output_folder = "output/verification", 
                           fastSpar_iterations = 20, bootstrap_num = 200, threads = 4)

FastSparCC_Server_Function(otu_table_input = "output/verification/Pacific_FL_table.tsv", 
                           output_folder = "output/verification", 
                           fastSpar_iterations = 20, bootstrap_num = 200, threads = 4)

# The values defined for the complete_cluster_pipeline() have been selected based on the
# "find_parameters = T" setting (see publication for details)
cluster_Combined <- import_sparcc_network(cor_file = "output/Verification/output_files/fastspar_cor_Combined_FL_table.tsv", 
                                        pval_file = "output/Verification/output_files/fastspar_pvalues_Combined_FL_table.tsv", 
                                        min_r = 0.4, min_p = 0.05) %>%
  filter(Type == "Positive") %>%
#  complete_cluster_pipeline(find_parameters = T)
  complete_cluster_pipeline(find_parameters = F, res_param = c(0.8, 2), cluster_num = 8,
                            cluster_method = "Cluster_network", vertex_num_drop = 10)

cluster_Pacific <- import_sparcc_network(cor_file = "output/Verification/output_files/fastspar_cor_Pacific_FL_table.tsv", 
                                          pval_file = "output/Verification/output_files/fastspar_pvalues_Pacific_FL_table.tsv", 
                                          min_r = 0.4, min_p = 0.05) %>%
  filter(Type == "Positive") %>%
  #complete_cluster_pipeline(find_parameters = T)
  complete_cluster_pipeline(find_parameters = F, res_param = c(0.8, 2), cluster_num = 8,
                            cluster_method = "Cluster_network", vertex_num_drop = 10)

cluster_Atlantic <- import_sparcc_network(cor_file = "output/Verification/output_files/fastspar_cor_Atlantic_FL_table.tsv", 
                                          pval_file = "output/Verification/output_files/fastspar_pvalues_Atlantic_FL_table.tsv", 
                                          min_r = 0.4, min_p = 0.05) %>%
  filter(Type == "Positive") %>%
  #complete_cluster_pipeline(find_parameters = T)
  complete_cluster_pipeline(find_parameters = F, res_param = c(0.8, 2), cluster_num = 8,
                            cluster_method = "Cluster_network", vertex_num_drop = 10)

cluster_tab_sparcc <- full_join(dplyr::rename(cluster_Combined$Cluster_Tab, "Combined_SparCC" = "Cluster"),
                                dplyr::rename(cluster_Atlantic$Cluster_Tab, "Atlantic_SparCC" = "Cluster")) %>%
  full_join(., dplyr::rename(cluster_Pacific$Cluster_Tab, "Pacific_SparCC" = "Cluster")) %>%
  dplyr::rename("OTU_ID" = "ASV") %>%
  mutate_if(is.numeric, function(x) paste0("Cluster_", x)) %>%
  filter(Combined_SparCC != "Cluster_NA" & Atlantic_SparCC != "Cluster_NA" & Pacific_SparCC != "Cluster_NA")

write_tsv(cluster_tab_sparcc, "output/Verification/Cluster_Tab_SparCC.tsv")

#### Run Pearson-Correlation Network Inference ####

count_matrix_Combined <- datalist_Combined %>%
  .$Count_Data %>%
  select(-c(1:8)) %>%
  as.matrix() 

count_matrix_Atlantic <- datalist_Atlantic %>%
  .$Count_Data %>%
  select(-c(1:8)) %>%
  as.matrix() 

count_matrix_Pacific <- datalist_Pacific %>%
  .$Count_Data %>%
  select(-c(1:8)) %>%
  as.matrix() 

## CLR-Transform count tables
clr_matrix_Combined <- SpiecEasi::clr(count_matrix_Combined+1) %>%
  t() %>%
  magrittr::set_colnames(datalist_Combined$Count_Data$OTU_ID)

clr_matrix_Atlantic <- SpiecEasi::clr(count_matrix_Atlantic+1) %>%
  t() %>%
  magrittr::set_colnames(datalist_Atlantic$Count_Data$OTU_ID)

clr_matrix_Pacific <- SpiecEasi::clr(count_matrix_Pacific+1) %>%
  t() %>%
  magrittr::set_colnames(datalist_Pacific$Count_Data$OTU_ID)

## Run Pearson Correlation
corr_matrix_Combined <- Hmisc::rcorr(clr_matrix_Combined, type = "pearson")
corr_matrix_Atlantic <- Hmisc::rcorr(clr_matrix_Atlantic, type = "pearson")
corr_matrix_Pacific <- Hmisc::rcorr(clr_matrix_Pacific, type = "pearson")

## Filter out all correlations with p >= 0.05 & r < 0.5
pearson_cooccurrence_Combined <- corr_matrix_Combined$r
pearson_cooccurrence_Combined[corr_matrix_Combined$P > 0.05 | corr_matrix_Combined$r < 0.5] <- 0

pearson_cooccurrence_Atlantic <- corr_matrix_Atlantic$r
pearson_cooccurrence_Atlantic[corr_matrix_Atlantic$P > 0.05 | corr_matrix_Atlantic$r < 0.5] <- 0

pearson_cooccurrence_Pacific <- corr_matrix_Pacific$r
pearson_cooccurrence_Pacific[corr_matrix_Pacific$P > 0.05 | corr_matrix_Pacific$r < 0.5] <- 0

# The values defined for the complete_cluster_pipeline() have been selected based on the
# "find_parameters = T" setting (see publication for details)
cluster_combined <- pearson_cooccurrence_Combined %>%
  reshape2::melt() %>%
  magrittr::set_colnames(c("From", "To", "weight")) %>%
  as_tibble() %>%
  filter(weight > 0 & weight < 1) %>%
  complete_cluster_pipeline(find_parameters = F, res_param = c(0.5, 1), cluster_method = "Cluster_network",
                            cluster_num = 10, vertex_num_drop = 10) 

cluster_Atlantic <- pearson_cooccurrence_Atlantic %>%
  reshape2::melt() %>%
  magrittr::set_colnames(c("From", "To", "weight")) %>%
  as_tibble() %>%
  filter(weight > 0 & weight < 1) %>%
  complete_cluster_pipeline(find_parameters = F, res_param = c(0.6, 1.2), cluster_method = "Cluster_network",
                            cluster_num = 10, vertex_num_drop = 10) 

cluster_Pacific <- pearson_cooccurrence_Pacific %>%
  reshape2::melt() %>%
  magrittr::set_colnames(c("From", "To", "weight")) %>%
  as_tibble() %>%
  filter(weight > 0 & weight < 1) %>%
  complete_cluster_pipeline(find_parameters = F, res_param = c(0.6, 1.2), cluster_method = "Cluster_network",
                            cluster_num = 10, vertex_num_drop = 10)

cluster_tab_spear <- full_join(dplyr::rename(cluster_combined$Cluster_Tab, "Combined_Spear" = "Cluster"), 
                               dplyr::rename(cluster_Atlantic$Cluster_Tab, "Atlantic_Spear" = "Cluster"), by = "ASV") %>%
  full_join(dplyr::rename(cluster_Pacific$Cluster_Tab, "Pacific_Spear" = "Cluster"), by = "ASV") %>%
  dplyr::rename("OTU_ID" = "ASV") %>%
  filter(!is.na(Atlantic_Spear) & !is.na(Pacific_Spear)) %>%
  mutate_if(is.numeric, function(x) paste0("Cluster_", x))  %>%
  filter(Combined_Spear != "Cluster_NA" & Atlantic_Spear != "Cluster_NA" & Pacific_Spear != "Cluster_NA")

write_tsv(cluster_tab_spear, "output/Verification/Cluster_Tab_Pearson.tsv")

#### Run SpiecEasi Network Inference ####

se_Combined <- spiec.easi(count_matrix_Combined, method='mb', lambda.min.ratio=1e-1,
                          nlambda=20, pulsar.params=list(rep.num=50, ncores = 8))

se_Atlantic <- spiec.easi(count_matrix_Atlantic, method='mb', lambda.min.ratio=1e-1,
                          nlambda=20, pulsar.params=list(rep.num=50, ncores = 8))

se_Pacific <- spiec.easi(count_matrix_Pacific, method='mb', lambda.min.ratio=1e-1,
                         nlambda=20, pulsar.params=list(rep.num=50, ncores = 8))

# Extract adjacency matrix
network_Combined_se <- se_Combined %>%
  getRefit() %>%
  adj2igraph() %>%
  set_vertex_attr(name = "name", value = colnames(count_matrix_Combined))

network_Atlantc_se <- se_Atlantic %>%
  getRefit() %>%
  adj2igraph() %>%
  set_vertex_attr(name = "name", value = colnames(count_matrix_Atlantic))

network_Pacific_se <- se_Pacific %>%
  getRefit() %>%
  adj2igraph() %>%
  set_vertex_attr(name = "name", value = colnames(count_matrix_Pacific))

## Use cluster_fast_greedy() for network clustering, since SpiecEasi doesnt yield correlation weights
cluster_tab_se <- bind_cols(OTU_ID = names(V(network_Combined_se)),
                         Combined_SE = paste0("Cluster_", cluster_fast_greedy(network_Combined_se)$membership)) %>%
  full_join(., bind_cols(OTU_ID = names(V(network_Atlantc_se)),
                         Atlantic_SE = paste0("Cluster_", cluster_fast_greedy(network_Atlantc_se)$membership)), by = "OTU_ID") %>%
  full_join(., bind_cols(OTU_ID = names(V(network_Pacific_se)),
                         Pacific_SE = paste0("Cluster_", cluster_fast_greedy(network_Pacific_se)$membership)), by = "OTU_ID") %>%
  filter(!is.na(Atlantic_SE) & !is.na(Pacific_SE))

write_tsv(cluster_tab_se, "output/Verification/Cluster_Tab_SpiecEasi.tsv")

#### ARI Comparison of clustering results ####

cluster_tab_sparcc <- read_tsv("output/Verification/Cluster_Tab_SparCC.tsv")
cluster_tab_spear <- read_tsv("output/Verification/Cluster_Tab_Pearson.tsv")
cluster_tab_se <- read_tsv("output/Verification/Cluster_Tab_SpiecEasi.tsv")

cluster_tab_all <- inner_join(cluster_tab_sparcc, cluster_tab_spear, by = "OTU_ID") %>%
  inner_join(., cluster_tab_se, by = "OTU_ID")

combinations <- expand.grid(expand.grid(c("Combined", "Atlantic", "Pacific"), c("SparCC", "Spear", "SE")) %>% 
  apply(., 1, function(x) paste(x, collapse = "_")),
  expand.grid(c("Combined", "Atlantic", "Pacific"), c("SparCC", "Spear", "SE")) %>% 
    apply(., 1, function(x) paste(x, collapse = "_")))

ARI_combinations <- purrr::map2(combinations$Var1, combinations$Var2, function(x, y) {
  ARI_tmp <- mclust::adjustedRandIndex(pull(cluster_tab_all, x), pull(cluster_tab_all, y))
  
  tibble(From = x, To = y, ARI = ARI_tmp)
}) %>%
  bind_rows() %>%
  filter(ARI != 1) %>%
  mutate(From = as.character(From),
         To = as.character(To)) %>%
  filter(From > To)

ARI_combinations %>%
  filter(To == "Combined_SparCC" | From == "Combined_SparCC") %>%
  arrange(desc(ARI))

#### 1. Cluster Composition ####

cluster_tmp <- cluster_tab_all

cluster_compositions <- purrr::map(unique(combinations$Var1), function(subset_tmp) {
  avg_composition <- datalist_Combined %>%
      mutate_count_datalist(function(x) x/sum(x)) %>%
      .$Count_Data %>%
      mutate(total = rowSums(select_if(., is.numeric))/ncol(select_if(., is.numeric))) %>%
      select(OTU_ID, total) %>%
      left_join(., select(cluster_tmp, OTU_ID, !!subset_tmp), by = "OTU_ID") %>%
      dplyr::rename("Cluster" = !!subset_tmp) %>%
      filter(!is.na(Cluster)) %>%
      with(., reshape2::dcast(data = ., formula = OTU_ID ~ Cluster, value.var = "total")) %>%
      tibble::as_tibble() %>%
      mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
      dplyr::rename_at(-1, function(x) paste0(subset_tmp, "_", x))
  return(avg_composition)
}) %>%
  reduce(full_join, by = "OTU_ID")

annotation <- tibble(DB_Type = colnames(cluster_compositions)[-1] %>%
                       str_replace_all(pattern = "_.*$", replacement = ""),
                     Method_Type = colnames(cluster_compositions)[-1] %>%
                       str_replace_all(pattern = "_Cluster.*$", replacement = "") %>%
                       str_replace_all(pattern = "^.*_", replacement = "")) %>%
  as.data.frame() %>%
  magrittr::set_rownames(colnames(cluster_compositions)[-1]) %>%
  mutate(Method_Type = str_replace_all(Method_Type, pattern = "Spear", replacement = "Pearson") %>%
           str_replace_all(pattern = "SE", replacement = "SpiecEasi"))

annotation_colors = list(DB_Type = magrittr::set_names(Palette3[1:3], unique(annotation$DB_Type)),
                         Method_Type = magrittr::set_names(cbbPalette[2:4], unique(annotation$Method_Type)))

BC_mat <- cluster_compositions %>%
  select(-1) %>%
  t() %>%
  vegan::vegdist() %>%
  as.matrix() %>%
  magrittr::set_colnames(colnames(cluster_compositions)[-1]) %>%
  magrittr::set_rownames(colnames(cluster_compositions)[-1])

pheat <- pheatmap::pheatmap(1-as.matrix(BC_mat), annotation_col = annotation, 
                   show_rownames = F, show_colnames = F, border_color = NA,
                   color = c(hcl.colors(n = 10, palette = "RdYlBu", rev = T), 
                             rep(hcl.colors(n = 10, palette = "RdYlBu", rev = T)[10], 1)),
                   annotation_row = annotation, annotation_colors = annotation_colors,
                   treeheight_row = 0, treeheight_col = 0)

ggsave(filename = "figs/Verification_Cluster_Stability_Ocean.png", plot = pheat, width = 8, height = 7, dpi = 300)

write.csv((1-as.matrix(BC_mat)), file = "output/Source_Data_Figure_5d.csv")

#### 2. Biogeographic Cluster Comparison ####

composition_dist <- cluster_compositions %>%
  mutate_if(is.numeric, function(x) x/sum(x)) %>%
  select(-1) %>%
  t() %>%
  vegan::vegdist() 

silhouette_analysis <- purrr::map_dbl(seq(2, 20), function(k) {
  membership_tmp <- hclust(composition_dist, method = "ward.D2") %>%
    cutree(k = k)
  mean(cluster::silhouette(membership_tmp, composition_dist)[,3])
})
# Based on visualization k = 7 yields optimal results
plot(2:20, silhouette_analysis)
  
cluster_cluster <- composition_dist %>%
  hclust(method = "ward.D2") %>% cutree(k = 7) %>%
  tibble::enframe() %>%
  magrittr::set_colnames(c("Group", "Cluster_Cluster")) %>%
  mutate(Dataset = str_replace_all(Group, pattern = "_.*$", replacement = "")) %>%
  mutate(Method = str_replace_all(Group, pattern = "_Cluster.*$", replacement = "") %>%
           str_replace_all(pattern = "^.*_", replacement = "")) %>%
  mutate(Group = paste0(Method, "_", Dataset, "_", str_replace_all(Group, pattern = "^.*_", replacement = "")))

datatable_final <- purrr::map(c("SparCC", "Spear", "SE"), function(method_tmp) {
  
  datatable_combined <- purrr::map(c("Combined", "Atlantic", "Pacific"), function(test_tmp) {
    
    data_subset_tmp <- paste0(test_tmp, "_", method_tmp)
    
    datalist_cluster <- datalist_Combined
    
    datalist_cluster$Count_Data <- datalist_cluster %>%
      mutate_count_datalist(function(x) x/sum(x)) %>%
      .$Count_Data %>%
      left_join(., select(cluster_tab_all, OTU_ID, data_subset_tmp), by = "OTU_ID") %>%
      dplyr::rename("Cluster" = !!data_subset_tmp) %>%
      filter(!is.na(Cluster)) %>%
      group_by(Cluster) %>%
      summarize_if(is.numeric, sum) %>%
      mutate(Cluster = as.character(Cluster))
    
    datatable_cluster <- datalist_cluster %>%
      create_datatable(., grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0) %>%
      mutate(Type = paste0(test_tmp)) %>%
      mutate(Group = paste0(method_tmp, "_", test_tmp, "_", str_replace_all(as.character(Group), pattern = "Cluster_", replacement = "")))
    
    return(datatable_cluster)
    
  }) %>%
    bind_rows()
  
  datatable_combined %>%
    mutate(Method = method_tmp) %>%
    return()
  
}) %>%
  bind_rows() %>%
  mutate(Method = str_replace_all(Method, pattern = "Spear", replacement = "Pearson") %>%
           str_replace_all(pattern = "SE", replacement = "SpiecEasi"))

p1 <- datatable_final %>%
  left_join(., select(cluster_cluster, Group, Cluster_Cluster), by = "Group") %>%
  filter(Ocean == "Pacific" & Depth < 100) %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster)) %>%
  group_by(Sample_ID, Method, Type) %>%
  mutate(Abundance = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Latitude)$Sample_ID))) %>% 
  
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Cluster_Cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 1.2) +
    scale_fill_manual(values = read_csv("../GRUMP_Model_Comparison/data/Cluster_Colour.csv")$Colour[c(2, 6, 10, 4, 5, 1, 3)]) +
    theme_bw() +
    facet_grid(Type ~ Method, scales = "free_x") +
    labs(x = "Latitude", y = "Rel. Abundance (%)", fill = "Cluster", title = "Pacific - <100m Depth") +
    theme(axis.text.x = element_blank())

p2 <- datatable_final %>%
  left_join(., select(cluster_cluster, Group, Cluster_Cluster), by = "Group") %>%
  filter(Ocean == "Atlantic" & Depth < 100) %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster)) %>%
  group_by(Sample_ID, Method, Type) %>%
  mutate(Abundance = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Latitude)$Sample_ID))) %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster, levels = c(2, 1, 3, 4, 5, 6, 7))) %>%
  
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Cluster_Cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 1.2) +
    scale_fill_manual(values = read_csv("../GRUMP_Model_Comparison/data/Cluster_Colour.csv")$Colour[c(2, 6, 10, 4, 5, 1, 3)]) +
    theme_bw() +
    facet_grid(Type ~ Method, scales = "free_x") +
    labs(x = "Latitude", y = "Rel. Abundance (%)", fill = "Cluster", title = "Atlantic - <100m Depth") +
    theme(axis.text.x = element_blank())

p3 <- datatable_final %>%
  left_join(., select(cluster_cluster, Group, Cluster_Cluster), by = "Group") %>%
  filter(Ocean == "Pacific" & Depth >= 100) %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster)) %>%
  group_by(Sample_ID, Method, Type) %>%
  mutate(Abundance = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Latitude)$Sample_ID))) %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster, levels = c(2, 1, 3, 4, 5, 6, 7))) %>%
  
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Cluster_Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 1.2) +
  scale_fill_manual(values = read_csv("../GRUMP_Model_Comparison/data/Cluster_Colour.csv")$Colour[c(2, 6, 10, 4, 5, 1, 3)]) +
  theme_bw() +
  facet_grid(Type ~ Method, scales = "free_x") +
  labs(x = "Latitude", y = "Rel. Abundance (%)", fill = "Cluster", title = "Pacific - >=100m Depth") +
  theme(axis.text.x = element_blank())

p4 <- datatable_final %>%
  left_join(., select(cluster_cluster, Group, Cluster_Cluster), by = "Group") %>%
  filter(Ocean == "Atlantic" & Depth >= 100) %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster)) %>%
  group_by(Sample_ID, Method, Type) %>%
  mutate(Abundance = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Latitude)$Sample_ID))) %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster, levels = c(2, 1, 3, 4, 5, 6, 7))) %>%
  
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Cluster_Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 1.2) +
  scale_fill_manual(values = read_csv("../GRUMP_Model_Comparison/data/Cluster_Colour.csv")$Colour[c(2, 6, 10, 4, 5, 1, 3)]) +
  theme_bw() +
  facet_grid(Type ~ Method, scales = "free_x") +
  labs(x = "Latitude", y = "Rel. Abundance (%)", fill = "Cluster", title = "Atlantic - >=100m Depth") +
  theme(axis.text.x = element_blank())

cowplot::plot_grid(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), 
                   p3 + theme(legend.position = "none"), p4 + theme(legend.position = "none"), ncol = 2)

ggsave("figs/Verification_Cluster_Biogeography_Ocean.png", width = 12, height = 12, dpi = 300)

ggsave("figs/Verification_Cluster_Biogeography_Ocean_Legend.png", cowplot::get_legend(p1))

#### 3. Environmental Prediction based on log-abundance-ratios ####

ratio_table <- datatable_final %>%
  left_join(., select(cluster_cluster, Group, Cluster_Cluster), by = "Group") %>%
  group_by(Sample_ID, Type, Method, Cluster_Cluster) %>%
  summarize(Abundance = sum(Abundance),
            Pot_Temperature = mean(Pot_Temperature)) %>%
  
  group_by(Sample_ID, Type, Method) %>% 
  mutate(ratio_1_2 = log(Abundance[Cluster_Cluster == 2] /Abundance[Cluster_Cluster == 1])) %>%
  filter(is.finite(ratio_1_2))

write_csv(ratio_table, "output/Source_Data_Figure_S7.csv")
  
ratio_table %>%
  ggplot(aes(x = Pot_Temperature, y = ratio_1_2)) +
    geom_point() +
    geom_smooth(method = "lm", col = "darkred") +
    facet_grid(Type ~ Method) +
    theme_bw() +
    labs(x = "Temperature (°C)", y = "Log. Abundance Ratio (Cohort 2 : Cohort 1)")

ggsave("figs/Verification_Ratio_Correlation_Ocean.png", width = 6.5, height = 6, dpi = 300)

ratio_table %>%
  group_by(Type, Method) %>%
  summarize(Adj_R2 = (function(x) summary(lm(x ~ Pot_Temperature))$adj.r.squared)(ratio_1_2),
            P_val = (function(x) {
              fstat <- summary(lm(x ~ Pot_Temperature))$fstatistic
              pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
              })(ratio_1_2),
            Slope = (function(x) lm(x ~ Pot_Temperature)$coefficients[2])(ratio_1_2),
            Intercept = (function(x) lm(x ~ Pot_Temperature)$coefficients[1])(ratio_1_2)) %>%
  arrange(Method)


#### 4. Mantel Test Comparison ####

method_vec <- rep(c("SparCC", "Spear", "SE"), 3)
subset_vec <- c(rep("Combined", 3), rep("Atlantic", 3), rep("Pacific", 3))

datalist_list <- purrr::map2(method_vec, subset_vec, function(method_tmp, test_tmp) {
  
  data_subset_tmp <- paste0(test_tmp, "_", method_tmp)
    
  datalist_cluster <- datalist_Combined
    
  datalist_cluster$Count_Data <- datalist_cluster %>%
    mutate_count_datalist(function(x) x/sum(x)) %>%
    .$Count_Data %>%
    left_join(., select(cluster_tab_all, OTU_ID, data_subset_tmp), by = "OTU_ID") %>%
    dplyr::rename("Cluster" = !!data_subset_tmp) %>%
    filter(!is.na(Cluster)) %>%
    group_by(Cluster) %>%
    summarize_if(is.numeric, sum) %>%
    mutate(Cluster = as.character(Cluster))
    
  return(datalist_cluster)
  
})

mat_names <- paste0(method_vec, "_", subset_vec)

gold_standard <- datalist_list[[1]]$Count_Data %>%
  select_if(is.numeric) %>%
  t() %>%
  vegan::vegdist()

mantel_comparison <- purrr::map(seq(2, length(mat_names)), function(ind) {
  
  test_tmp <- datalist_list[[ind]]$Count_Data %>%
    select_if(is.numeric) %>%
    t() %>%
    vegan::vegdist()
  
  mantel_tmp <- vegan::mantel(gold_standard, test_tmp, method = "pearson", permutations = 999)
  
  result <- tibble(Type = mat_names[ind], Method = method_vec[ind], Subset = subset_vec[ind],
                   Mantel_R = mantel_tmp$statistic, Mantel_pval = mantel_tmp$signif)
  
  return(result)
  
}) %>%
  bind_rows()

mantel_comparison %>%
  mutate(Method = str_replace_all(Method, pattern = "Spear", replacement = "Pearson") %>%
           str_replace_all(pattern = "SE", replacement = "SpiecEasi"))

#### 5. Gold Standard Comparison ####

cluster_table_tmp <- cluster_cluster %>%
  mutate(Cluster = str_replace_all(string = Group, pattern = "^.*_", replacement = "Cluster_")) %>%
  mutate(Cluster_Cluster = as.character(Cluster_Cluster))

gold_standard <- datalist_list[[1]]$Count_Data %>%
  with(., full_join(filter(cluster_table_tmp, Dataset == "Combined" & Method == "SparCC") %>% select(Cluster_Cluster, Cluster),
                     ., by = "Cluster")) %>%
  arrange(Cluster_Cluster) %>%
  group_by(Cluster_Cluster) %>%
  summarize_if(is.numeric, sum)

BC_comparison_result <- purrr::map(seq(2, length(mat_names)), function(ind) {
  
  test_tmp <- datalist_list[[ind]]$Count_Data %>%
    with(., right_join(filter(cluster_table_tmp, Dataset == subset_vec[ind] & Method == method_vec[ind]) %>% select(Cluster_Cluster, Cluster),
                       ., by = "Cluster")) %>%
    arrange(Cluster_Cluster) %>%
    group_by(Cluster_Cluster) %>%
    summarize_if(is.numeric, sum) 
  
  gold_standard_tmp <- full_join(gold_standard, tibble(Cluster_Cluster = test_tmp$Cluster_Cluster), by = "Cluster_Cluster") %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
    arrange(Cluster_Cluster) %>%
    select(-1)
  
  test_tmp <- full_join(test_tmp, tibble(Cluster_Cluster = gold_standard$Cluster_Cluster), by = "Cluster_Cluster") %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
    arrange(Cluster_Cluster) %>%
    select(-1)
  
  bc_comparison <- purrr::map_dbl(seq(1, ncol(gold_standard_tmp)), function(sample_tmp) {
    
    x <- pull(gold_standard_tmp, sample_tmp)
    y <- pull(test_tmp, sample_tmp)
    1- sum(abs(x - y)) / sum(x + y)
    
  })
  
  result <- tibble(Type = mat_names[ind], Method = method_vec[ind], Subset = subset_vec[ind],
                   BC = bc_comparison)
  
  return(result)
  
}) %>%
  bind_rows()

library(ggridges)

BC_comparison_result %>%
  mutate(Method = str_replace_all(Method, pattern = "Spear", replacement = "Pearson") %>%
           str_replace_all(pattern = "SE", replacement = "SpiecEasi")) %>% 
  ggplot(aes(x = BC, y = Type, fill = Method)) +
    geom_density_ridges() +
    scale_fill_manual(values = cbbPalette[2:4]) +
    scale_y_discrete(labels = c("Pacific", "Combined", "Atlantic", "Pacific", "Atlantic", "Pacific", "Combined", "Atlantic")) +
    theme_ridges() +
    lims(x = c(0, 1)) +
    labs(x = "Bray-Curtis Similarity", y = "")

ggsave("figs/Verification_Ocean_BC_single.png", width = 5, height = 4, dpi = 300)  

            