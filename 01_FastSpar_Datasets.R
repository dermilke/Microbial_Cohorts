#### This script runs FastSpar on all datasets of the manuscript ####
#
# Author: Felix Milke
# Date: 26.02.2025

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Fastspar_new.R")

## Wrapper-function to run on each dataset.
# It selects the maximum number of ASVs for each table using average abundances
# And formats the count-data table to be readable for the FastSpar program.
# Then it runs the FastSpar script, which automatically runs all FastSpar
# functions. Because FastSpar is written in C++ it is executed via the shell,
# using the local conda-environment in which it is installed.
# (for further information, see the FastSpar.R script)
fastSpar_wrapper <- function(datalist, max_organisms, dataset_name,
                                      fastspar_bootstrap = 500, fastspar_iter = 20, threads = 2) {
  
  set.seed(123)
  
  chosen_asvs <- datalist %>%
    mutate_count_datalist(function(x) x/sum(x)) %>%
    .$Count_Data %>%
    mutate(total = rowMeans(select_if(., is.numeric))) %>%
    select(OTU_ID, total) %>%
    arrange(desc(total)) %>%
    dplyr::slice(1:max_organisms) %>%
    .$OTU_ID
  
  datalist_filtered <- datalist %>%
    filter_taxa_datalist(OTU_ID %in% chosen_asvs)
  
  datalist_filtered %>%
    .$Count_Data %>%
    select(-c(2:8)) %>%
    rename(`#OTU ID` = "OTU_ID") %>%
    write_tsv(paste0("output/fastspar/", dataset_name, "_table.tsv"))
    
  message("Running FastSpar. This may take some time (hours)...")
    
  FastSparCC_Server_Function(otu_table_input = paste0("output/fastspar/", dataset_name, "_table.tsv"), 
                             output_folder = "output/fastspar", 
                             fastSpar_iterations = fastspar_iter, bootstrap_num = fastspar_bootstrap, threads = threads)
    
}

# Define maximum number of ASVs for FastSpar calculations
max_organisms <- 5000

### Soil: Australia ####

datalist_australia <- list(Count_Data = read_tsv("../Australia_Soils/data/Count_Data_Filtered/Processed/Prok/Full_Prok_Count.tsv"),
                           Meta_Data = read_tsv("../Australia_Soils/data/Count_Data_Filtered/Meta_Data/Prok/Meta_Data.tsv") %>%
                             left_join(., read_tsv("../Australia_Soils/data/worldclim_data/bioclimatic.tsv"), by = "Sample_ID"))

fastSpar_wrapper(datalist_australia, max_organisms, dataset_name = "Soil_Australia",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("datalist_australia")

### Freshwater: Lake Mendota ####

datalist_lake <- import_data("../Lake_Mendota/data/Count_Data/Complete/", min_counts = 2000, kingdom = "Prok", abundance_filter = F)

fastSpar_wrapper(datalist_lake, max_organisms, dataset_name = "Lake_Mendota",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("datalist_lake")

### Marine: Comparison ####

datalist_comparison <- import_data("../Archive/Comparison_Paper/data/Combined_Reduced/", min_counts = 2000, kingdom = "Prok", abundance_filter = F)

fastSpar_wrapper(datalist_comparison, max_organisms, dataset_name = "Comparison",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 4)

rm("datalist_comparison")

### Marine: SPOT ####

datalist_spot <- import_data("../SPOT/data/Count_Data/Complete/", min_counts = 2000, kingdom = "Prok", abundance_filter = F)

fastSpar_wrapper(datalist_spot, max_organisms, dataset_name = "SPOT",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 4)

rm("datalist_spot")

### Marine: CalCOFI ####

datalist_calcofi <- import_data("../CalCOFI/data/Count_Data/Complete/", min_counts = 2000, kingdom = "Prok", abundance_filter = F)

fastSpar_wrapper(datalist_calcofi, max_organisms, dataset_name = "CalCOFI",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 4)

rm("datalist_calcofi")

#### EMP datasets ####

data_secretion <- import_data("../EMP/data/Count_Data/Animal_secretion/", kingdom = "Prok", min_counts = 2000) %>%
  mutate_meta_datalist(Type = ifelse(grepl(pattern = "saliva", x = Description), "Saliva",
                                     ifelse(grepl(pattern = "Nose", x = Description), "Nose",
                                            ifelse(grepl(pattern = "oral", x = Description) | grepl(pattern = "mouth", x = Description), "Mouth", "Various")))) 
    
#### Secretion: Nose ####

data_secretion_human_nose <- data_secretion %>%
  filter_station_datalist(Type == "Nose" & host_common_name_provided == "human")
    
fastSpar_wrapper(data_secretion_human_nose, max_organisms, dataset_name = "Nose_Secretion",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_secretion_human_nose")

#### Secretion: Saliva ####

data_secretion_human_saliva <- data_secretion %>%
  filter_station_datalist(Type == "Saliva" & host_common_name_provided == "human")

fastSpar_wrapper(data_secretion_human_saliva, max_organisms, dataset_name = "Human_Saliva",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)
    
rm("data_secretion_human_saliva", "data_secretion")

#### Distalgut: Human ####

data_distal_gut <- import_data("../EMP/data/Count_Data/Animal_distal_gut/", kingdom = "Prok", min_counts = 2000)
    
data_distal_gut_human <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "human")

fastSpar_wrapper(data_distal_gut_human, max_organisms, dataset_name = "Distalgut_Human",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_distal_gut_human")
    
#### Distalgut: Kangaroo ####

data_distal_gut_kangaroo <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "Kangaroo")

fastSpar_wrapper(data_distal_gut_kangaroo, max_organisms, dataset_name = "Distalgut_Kangaroo",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_distal_gut_kangaroo")

#### Distalgut: Rabbit ####
    
data_distal_gut_rabbit <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "Rabbit") 

fastSpar_wrapper(data_distal_gut_rabbit, max_organisms, dataset_name = "Distalgut_Rabbit",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_distal_gut_rabbit")

#### Distalgut: Deer ####

data_distal_gut_deer <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "Sambar deer")

fastSpar_wrapper(data_distal_gut_deer, max_organisms, dataset_name = "Distalgut_Deer",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_distal_gut_deer")

#### Distalgut: Monkey ####

data_distal_gut_monkey <- data_distal_gut %>%
      filter_station_datalist(host_common_name_provided == "Spider monkey")

fastSpar_wrapper(data_distal_gut_monkey, max_organisms, dataset_name = "Distalgut_Monkey",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_distal_gut_monkey", "data_distal_gut")

#### Skin: Foot ####

data_animal_surface <- import_data("../EMP/data/Count_Data/Animal_surface/", kingdom = "Prok", min_counts = 2000) %>%
  mutate_meta_datalist(Bodypart = str_replace_all(Description, pattern = "^.*\\.", replacement = "") %>%
                         str_replace_all(pattern = "sample_[0-9]* ", replacement = ""))
    
data_surface_human_foot <- data_animal_surface %>%
  filter_station_datalist(host_common_name_provided == "human" & Bodypart == "Foot")

fastSpar_wrapper(data_surface_human_foot, max_organisms, dataset_name = "Surface_Foot",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_surface_human_foot")

#### Skin: Hand ####

data_surface_human_hand <- data_animal_surface %>%
  filter_station_datalist(host_common_name_provided == "human" & Bodypart == "Hand")

fastSpar_wrapper(data_surface_human_hand, max_organisms, dataset_name = "Surface_Hand",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_surface_human_hand", "data_animal_surface")

#### Plant: Rhizosphere ####

data_plant_rhizo_root <- import_data("../EMP/data/Count_Data/Plant_rhizosphere/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria")
    
data_plant_rhizosphere <- data_plant_rhizo_root %>%
  filter_station_datalist(Description == "Rhizosphere")

fastSpar_wrapper(data_plant_rhizosphere, max_organisms, dataset_name = "Plant_Rhizosphere",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_plant_rhizosphere")

#### Plant: Root ####

data_plant_roots <- data_plant_rhizo_root %>%
  filter_station_datalist(Description == "Roots")

fastSpar_wrapper(data_plant_roots, max_organisms, dataset_name = "Plant_Roots",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_plant_roots", "data_plant_rhizo_root")

#### Plant: Surface ####

data_plant_surface <- import_data("../EMP/data/Count_Data/Plant_surface/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria")

fastSpar_wrapper(data_plant_surface, max_organisms, dataset_name = "Plant_Surface",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_plant_surface")

#### Soil: Sandfilter ####

data_soil <- import_data("../EMP/data/Count_Data/Soil_non_saline/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria")
    
data_soil_sandfilter <- data_soil %>%
  filter_station_datalist(Description == "samples of sand from slow_sand filter water purification system" & depth_m <= 0.4)

fastSpar_wrapper(data_soil_sandfilter, max_organisms, dataset_name = "Soil_Sandfilter",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_soil_sandfilter")

#### Soil: Cultivated ####

data_soil_cultivated_habitat <- data_soil %>%
  filter_station_datalist(study_id == 1721 & collection_timestamp == "16.11.11")

fastSpar_wrapper(data_soil_cultivated_habitat, max_organisms, dataset_name = "Soil_Cultivated",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_soil_cultivated_habitat")

#### Soil: Field ####

data_soil_field <- data_soil %>%
  filter_station_datalist(study_id == 990)

fastSpar_wrapper(data_soil_field, max_organisms, dataset_name = "Soil_Field",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_soil_field")

#### Soil: Garden ####

data_soil_garden <- data_soil %>%
  filter_station_datalist(study_id == 1674)

fastSpar_wrapper(data_soil_garden, max_organisms, dataset_name = "Soil_Garden",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_soil_garden", "data_soil")

#### Indoor surface ####

data_surface_nonsaline_indoor <- import_data("../EMP/data/Count_Data/Surface_non_saline/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria") %>%
  filter_station_datalist(study_id == 2192)

fastSpar_wrapper(data_surface_nonsaline_indoor, max_organisms, dataset_name = "Surface_Nonsaline",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_surface_nonsaline_indoor")

#### Freshwater: Timeseries ####

data_water_nonsaline <- import_data("../EMP/data/Count_Data/Water_non_saline/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria")
    
data_water_nonsaline_timeseries <- data_water_nonsaline %>%
  filter_station_datalist(study_id == 1288)

fastSpar_wrapper(data_water_nonsaline_timeseries, max_organisms, dataset_name = "Water_Nonsaline_Timeseries",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_water_nonsaline_timeseries")

#### Freshwater: USA ####

data_water_nonsaline_usa <- data_water_nonsaline %>%
  filter_station_datalist(study_id == 1883 & env_feature == "freshwater habitat" & country == "GAZ:United States of America")

fastSpar_wrapper(data_water_nonsaline_usa, max_organisms, dataset_name = "Water_Nonsaline_USA",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_water_nonsaline_usa")

#### Freshwater: Germany ####

data_water_nonsaline_germany <- data_water_nonsaline %>%
  filter_station_datalist(study_id == 945)

fastSpar_wrapper(data_water_nonsaline_germany, max_organisms, dataset_name = "Water_Nonsaline_Germany",
                 fastspar_bootstrap = 500, fastspar_iter = 20, threads = 8)

rm("data_water_nonsaline_germany", "data_water_nonsaline")
