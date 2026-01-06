#### Validation Soil: Australia - Random Subsets ####

## This script runs the validation of clustering results using data subsets and three different network inference
## algorithms (SparCC, SpiecEasi, Pearson)
# First create random data subsets, store them for future use, and find all shared ASVs (Comparison between 
# Validation runs is more reliable and fair when considering only shared ASVs).
# Then compute:
# 1. SpiecEasi
# 2. SparCC
# 3. Pearson correlation using CLR transformation
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

source("R/Fastspar_Studio.R")
source("R/Cluster_Pipeline.R")
source("R/Import_SparCC_Network.R")

source("R/Names_Colors_Improvement.R")

threads <- 4

#### Construct Test-Data Subsets ####
# only run once to keep the sampleset
number_of_test_splits <- 3
size_of_test_splits <- 0.7
min_organisms = 5000

datalist_Australia <- import_data("../Australia_Soils/data/Count_Data_Filtered/", kingdom = "Prok")

datalist_Australia$Meta_Data <- datalist_Australia$Meta_Data %>%
  left_join(., read_tsv("../Australia_Soils/data/worldclim_data/bioclimatic.tsv"))

set.seed(123)
test_set <- map(seq(1, number_of_test_splits), function(x) {
  sample(datalist_Australia$Meta_Data$Sample_ID, size = ceiling(length(datalist_Australia$Meta_Data$Sample_ID) * size_of_test_splits))
})

test_set[[(number_of_test_splits+1)]] <- datalist_Australia$Meta_Data$Sample_ID
number_of_test_splits <- number_of_test_splits + 1

test_set <- write_rds(test_set, file = "data/Verification_Australia_Test_Split.rds")

test_set <- read_rds("data/Verification_Australia_Test_Split.rds")

available_in_all <- purrr::map(seq(1, number_of_test_splits), function(test_split_num) {

  count_data_selected <- datalist_Australia %>%
    filter_station_datalist(Sample_ID %in% test_set[[!!test_split_num]]) %>%
    .$Count_Data %>%
    .$OTU_ID
  
}) %>%
  reduce(intersect)

chosen_asvs <- datalist_Australia$Count_Data %>%
  filter(OTU_ID %in% available_in_all) %>%
  mutate_if(is.numeric, function(x) x/sum(x)) %>%
  mutate(total = rowMeans(select_if(., is.numeric))) %>%
  arrange(desc(total)) %>%
  dplyr::slice(1:min_organisms) %>%
  .$OTU_ID

#### Run SpiecEasi Network Inference ####
purrr::map(seq(1, number_of_test_splits), function(test_split_num) {
  
  datalist_test <- datalist_Australia %>%
    filter_station_datalist(Sample_ID %in% test_set[[!!test_split_num]]) %>%
    filter_taxa_datalist(OTU_ID %in% chosen_asvs)
  
  count_matrix <- datalist_test %>%
    .$Count_Data %>%
    select(-c(1:8)) %>%
    as.matrix() %>%
    t() %>%
    magrittr::set_colnames(datalist_test$Count_Data$OTU_ID)
  
  se_test <- spiec.easi(count_matrix, method='mb', lambda.min.ratio=1e-1,
                        nlambda=20, pulsar.params=list(rep.num=50, ncores = threads))
  
  # Extract adjacency matrix
  network_se <- se_test %>%
    getRefit() %>%
    adj2igraph() %>%
    igraph::set_vertex_attr(name = "name", value = colnames(count_matrix))
  
  bind_cols(OTU_ID = names(igraph::V(network_se)),
            !!paste0("SpiecEasi_Australia_", test_split_num) := paste0("Cluster_", igraph::cluster_fast_greedy(network_se)$membership)) %>%
    write_tsv(paste0("output/Verification/Cluster_Tab_Australia_", test_split_num, "_SpiecEasi.tsv"))
  
})

#### Run SparCC and Pearson Correlation Network Inference ####
purrr::map(seq(1, number_of_test_splits), function(test_split_num) {

  datalist_test <- datalist_Australia %>%
    filter_station_datalist(Sample_ID %in% test_set[[!!test_split_num]]) %>%
    filter_taxa_datalist(OTU_ID %in% chosen_asvs)
  
  datalist_test %>%
    .$Count_Data %>%
    select(-c(2:8)) %>%
    rename(`#OTU ID` = "OTU_ID") %>%
    write_tsv(paste0("output/Verification/Australia_", test_split_num, "_table.tsv"))
  
  #### SparCC ####
  
  message("Running FastSpar. This may take some time (hours)...")
  
  FastSparCC_Server_Function(otu_table_input = paste0("output/verification/Australia_", test_split_num, "_table.tsv"), 
                             output_folder = "output/verification", 
                             fastSpar_iterations = 20, bootstrap_num = 200, threads = threads)
  
  # The values defined for the complete_cluster_pipeline() have been selected based on the
  # "find_parameters = T" setting (see publication for details)
  cluster_sparcc <- import_sparcc_network(cor_file = paste0("output/Verification/output_files/fastspar_cor_Australia_", test_split_num, "_table.tsv"),
                                          pval_file = paste0("output/Verification/output_files/fastspar_pvalues_Australia_", test_split_num, "_table.tsv"),
                                          min_r = 0.4, min_p = 0.05) %>%
    filter(Type == "Positive") %>%
    #complete_cluster_pipeline(find_parameters = T)
    complete_cluster_pipeline(find_parameters = F, res_param = c(1, 2), cluster_num = 8,
                              cluster_method = "Cluster_network", vertex_num_drop = 10)
    
  cluster_sparcc$Cluster_Tab %>%
    dplyr::rename("OTU_ID" = "ASV") %>%
    mutate(Cluster = paste0("Cluster_", Cluster)) %>%
    dplyr::rename(!!paste0("SparCC_Australia_", test_split_num) := "Cluster") %>%
    write_tsv(paste0("output/Verification/Cluster_Tab_Australia_", test_split_num, "_SparCC.tsv"))
  
  #### CLR-transformed Pearson Correlation ####
  
  count_matrix <- datalist_test %>%
    .$Count_Data %>%
    select(-c(1:8)) %>%
    as.matrix() 
  
  clr_matrix <- clr(count_matrix+1) %>%
    t() %>%
    magrittr::set_colnames(datalist_test$Count_Data$OTU_ID)
  
  corr_matrix <- Hmisc::rcorr(clr_matrix, type = "pearson")
  
  pearson_cooccurrence <- corr_matrix$r
  pearson_cooccurrence[corr_matrix$P > 0.05 | corr_matrix$r < 0.6] <- 0
  
  # The values defined for the complete_cluster_pipeline() have been selected based on the
  # "find_parameters = T" setting (see publication for details)
  cluster_spear <- pearson_cooccurrence %>%
    reshape2::melt() %>%
    magrittr::set_colnames(c("From", "To", "weight")) %>%
    as_tibble() %>%
    filter(weight > 0 & weight < 1) %>%
    complete_cluster_pipeline(find_parameters = F, res_param = c(.5, 1.5), cluster_method = "Cluster_network",
                              cluster_num = 8, vertex_num_drop = 10) 
  
  cluster_spear$Cluster_Tab %>%
    dplyr::rename("OTU_ID" = "ASV") %>%
    mutate(Cluster = paste0("Cluster_", Cluster)) %>%
    dplyr::rename(!!paste0("Pearson_Australia_", test_split_num) := "Cluster") %>%
    write_tsv(paste0("output/Verification/Cluster_Tab_Australia_", test_split_num, "_Pearson.tsv"))
  
})

#### ARI Comparison of clustering results ####

cluster_tab_sparcc <- purrr::map(seq(1, number_of_test_splits), function(x) read_tsv(paste0("output/Verification/Cluster_Tab_Australia_", x, "_SparCC.tsv"))) %>%
  reduce(full_join, by = "OTU_ID")
cluster_tab_spear <- purrr::map(seq(1, number_of_test_splits), function(x) read_tsv(paste0("output/Verification/Cluster_Tab_Australia_", x, "_Pearson.tsv"))) %>%
  reduce(full_join, by = "OTU_ID")
cluster_tab_se <- purrr::map(seq(1, number_of_test_splits), function(x) read_tsv(paste0("output/Verification/Cluster_Tab_Australia_", x, "_SpiecEasi.tsv"))) %>%
  reduce(full_join, by = "OTU_ID")

cluster_tab_all <- inner_join(cluster_tab_sparcc, cluster_tab_spear, by = "OTU_ID") %>%
  inner_join(., cluster_tab_se, by = "OTU_ID")

combinations <- expand.grid(expand.grid(c("SparCC", "Pearson", "SpiecEasi"), paste0("Australia_", seq(1, number_of_test_splits))) %>% 
                              apply(., 1, function(x) paste(x, collapse = "_")),
                            expand.grid(c("SparCC", "Pearson", "SpiecEasi"), paste0("Australia_", seq(1, number_of_test_splits))) %>% 
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

ARI_combinations %>% view()

#### 1. Cluster Composition ####

cluster_tmp <- cluster_tab_all

cluster_compositions <- purrr::map(unique(combinations$Var1), function(subset_tmp) {
  avg_composition <- datalist_Australia %>%
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
  reduce(full_join, by = "OTU_ID") %>%
  mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x))

annotation <- tibble(Method_Type = colnames(cluster_compositions)[-1] %>%
                       str_replace_all(pattern = "_.*$", replacement = ""),
                     DB_Type = colnames(cluster_compositions)[-1] %>%
                       str_replace_all(pattern = "_Cluster.*$", replacement = "") %>%
                       str_replace_all(pattern = "^.*_", replacement = "")) %>%
  as.data.frame() %>%
  magrittr::set_rownames(colnames(cluster_compositions)[-1])

annotation_colors = list(DB_Type = magrittr::set_names(Palette3[1:4], unique(annotation$DB_Type)),
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

ggsave(filename = "figs/Verification_Cluster_Stability_Australia.png", plot = pheat, width = 8, height = 7, dpi = 300)

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
# Based on visualization k = 8 yields optimal results
plot(2:20, silhouette_analysis)

cluster_cluster <- composition_dist %>%
  hclust(method = "ward.D2") %>% cutree(k = 8) %>%
  tibble::enframe() %>%
  magrittr::set_colnames(c("Group", "Cluster_Cluster")) %>%
  mutate(Method = str_replace_all(Group, pattern = "_.*$", replacement = "")) %>%
  mutate(Dataset = str_replace_all(Group, pattern = "_Cluster.*$", replacement = "") %>%
           str_replace_all(pattern = "^.*_", replacement = "")) %>%
  mutate(Group = paste0(Method, "_", Dataset, "_", str_replace_all(Group, pattern = "^.*_", replacement = "")))

datatable_final <- purrr::map(c("SparCC", "Pearson", "SpiecEasi"), function(method_tmp) {
  
  datatable_combined <- purrr::map(seq(1, number_of_test_splits), function(test_tmp) {
    
    data_subset_tmp <- paste0(method_tmp, "_Australia_", test_tmp)
    
    datalist_cluster_Australia <- datalist_Australia
    
    datalist_cluster_Australia$Count_Data <- datalist_cluster_Australia %>%
      mutate_count_datalist(function(x) x/sum(x)) %>%
      .$Count_Data %>%
      left_join(., select(cluster_tab_all, OTU_ID, data_subset_tmp), by = "OTU_ID") %>%
      dplyr::rename("Cluster" = !!data_subset_tmp) %>%
      filter(!is.na(Cluster)) %>%
      group_by(Cluster) %>%
      summarize_if(is.numeric, sum) %>%
      mutate(Cluster = as.character(Cluster))
    
    datatable_cluster_Australia <- datalist_cluster_Australia %>%
      create_datatable(., grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0) %>%
      mutate(Type = paste0("Test_Split_", test_tmp)) %>%
      mutate(Group = paste0(method_tmp, "_", test_tmp, "_", str_replace_all(as.character(Group), pattern = "Cluster_", replacement = "")))
    
    return(datatable_cluster_Australia)
    
  }) %>%
    bind_rows()
  
  datatable_combined %>%
    mutate(Method = method_tmp) %>%
    return()
  
}) %>%
  bind_rows()

datatable_final %>%
  left_join(., select(cluster_cluster, Group, Cluster_Cluster), by = "Group") %>%
  mutate(Cluster_Cluster = ordered(Cluster_Cluster)) %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., temp_max, ph)$Sample_ID))) %>%
  mutate(Type = str_replace_all(Type, pattern = "Test_Split_4", replacement = "Complete Dataset")) %>%
  group_by(Sample_ID, Method, Type) %>%
  mutate(Abundance = Abundance/sum(Abundance)) %>%
  ungroup() %>% 
  mutate(Type = str_replace_all(Type, pattern = "Test_Split_", replacement = "Subset ")) %>%
  
  ggplot(., aes(x = Sample_ID, y = Abundance*100, fill = Cluster_Cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 1.2) +
    scale_fill_manual(values = c(cbbPalette[c(2, 3, 7, 4, 5, 6, 8)], "red")) +
    theme_bw() +
    facet_grid(Type ~ Method, scales = "free_x") +
    labs(x = "Max. annual temperature (°C)", y = "Rel. Abundance (%)", fill = "Cluster") +
    theme(axis.text.x = element_blank())

ggsave("figs/Verification_Cluster_Biogeography_Australia.png", width = 8, height = 7, dpi = 300)

#### 3. Environmental Prediction based on log-abundance-ratios ####

ratio_table <- datatable_final %>%
  left_join(., select(cluster_cluster, Group, Cluster_Cluster), by = "Group") %>%
  group_by(Sample_ID, Type, Method, Cluster_Cluster) %>%
  summarize(Abundance = sum(Abundance),
            Temp_Max = mean(temp_max)) %>%
  
  group_by(Sample_ID, Type, Method) %>% 
  mutate(ratio_1_2 = log(Abundance[Cluster_Cluster == 2] /Abundance[Cluster_Cluster == 5])) %>%
  filter(is.finite(ratio_1_2))

ratio_table %>%
  mutate(Type = str_replace_all(Type, pattern = "Test_Split_", replacement = "Subset ")) %>%
  ggplot(aes(x = Temp_Max, y = ratio_1_2)) +
  geom_point(size = .75) +
  geom_smooth(method = "lm", col = "darkred") +
  facet_grid(Type ~ Method) +
  theme_bw() +
  labs(x = "Max. annual temperature (°C)", y = "Log. Abundance Ratio (Cohort 2 : Cohort 5)")

ggsave("figs/Verification_Ratio_Correlation_Soil.png", width = 6.5, height = 6, dpi = 300)

ratio_table %>%
  group_by(Type, Method) %>%
  summarize(Adj_R2 = (function(x) summary(lm(x ~ Temp_Max))$adj.r.squared)(ratio_1_2),
            P_val = (function(x) {
              fstat <- summary(lm(x ~ Temp_Max))$fstatistic
              pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
            })(ratio_1_2),
            Slope = (function(x) lm(x ~ Temp_Max)$coefficients[2])(ratio_1_2),
            Intercept = (function(x) lm(x ~ Temp_Max)$coefficients[1])(ratio_1_2)) %>%
  arrange(Method)

#### 4. Mantel Test Comparison ####

method_vec <- rep(c("SparCC", "Pearson", "SpiecEasi"), 4)
subset_vec <- c(rep("Australia_1", 3), rep("Australia_2", 3), rep("Australia_3", 3), rep("Australia_4", 3))

datalist_list <- purrr::map2(method_vec, seq(1, length(subset_vec)), function(method_tmp, test_tmp_ind) {
  
  data_subset_tmp <- paste0(method_tmp, "_", subset_vec[test_tmp_ind])
  
  datalist_cluster <- datalist_Australia
  
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

gold_standard <- datalist_list[[10]]$Count_Data %>%
  select_if(is.numeric) %>%
  t() %>%
  vegan::vegdist()

mantel_comparison <- purrr::map(seq(1, length(mat_names))[-10], function(ind) {
  
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
  ggplot(aes(x = Method, y = Mantel_R)) +
    geom_boxplot() +
    labs(y = 'Pearson Mantel Correlation with "Gold Standard"')

ggsave("figs/Validation_Australia_Mantel_Method.png", width = 4.5, height = 5, dpi = 300)

#### 5. Gold Standard Comparison ####

cluster_table_tmp <- cluster_cluster %>%
  mutate(Dataset = str_replace_all(Group, pattern = "_[0-9]$", replacement = "") %>%
           str_replace_all(pattern = "^.*_", replacement = "Australia_")) %>%
  mutate(Method = str_replace_all(Group, pattern = "_.*$", replacement = "")) %>%
  mutate(Cluster = str_replace_all(string = Group, pattern = "^.*_", replacement = "Cluster_")) %>%
  mutate(Cluster_Cluster = as.character(Cluster_Cluster))

gold_standard <- datalist_list[[10]]$Count_Data %>%
  with(., full_join(filter(cluster_table_tmp, Dataset == "Australia_4" & Method == "SparCC") %>% select(Cluster_Cluster, Cluster),
                    ., by = "Cluster")) %>%
  arrange(Cluster_Cluster) %>%
  group_by(Cluster_Cluster) %>%
  summarize_if(is.numeric, sum)

BC_comparison_result <- purrr::map(seq(1, length(mat_names))[-10], function(ind) {
  
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
  ggplot(aes(x = BC, y = Type, fill = Method)) +
  geom_density_ridges() +
  scale_y_discrete(labels = rev(c(paste0("Subset ", 1:4), paste0("Subset ", 1:4), paste0("Subset ", 1:3)))) +
  scale_fill_manual(values = cbbPalette[2:4]) +
  theme_ridges() +
  lims(x = c(0, 1))

ggsave("figs/Verification_Soil_BC_single.png", width = 5.5, height = 5, dpi = 300)
