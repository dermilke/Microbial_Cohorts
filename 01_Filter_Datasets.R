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

### Soil: Australia ####

datalist_australia <- list(Count_Data = read_tsv("../Australia_Soils/data/Count_Data_Filtered/Processed/Prok/Full_Prok_Count.tsv"),
                           Meta_Data = read_tsv("../Australia_Soils/data/Count_Data_Filtered/Meta_Data/Prok/Meta_Data.tsv") %>%
                             left_join(., read_tsv("../Australia_Soils/data/worldclim_data/bioclimatic.tsv"), by = "Sample_ID")) %>%
  filter_wrapper(min_organisms = min_organisms)

write_tsv(datalist_australia$Count_Data, "data/Count_Data_Filtered/Soil_Australia/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(datalist_australia$Meta_Data, "data/Count_Data_Filtered/Soil_Australia/Meta_Data/Prok/Meta_Data.tsv")

rm("datalist_australia")

### Freshwater: Lake Mendota ####

datalist_lake <- import_data("../Lake_Mendota/data/Count_Data/Complete/", min_counts = 2000, kingdom = "Prok", abundance_filter = F) %>%
  filter_wrapper(min_organisms = min_organisms)

write_tsv(datalist_lake$Count_Data, "data/Count_Data_Filtered/Lake_Mendota/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(datalist_lake$Meta_Data, "data/Count_Data_Filtered/Lake_Mendota/Meta_Data/Prok/Meta_Data.tsv")

rm("datalist_lake")

### Marine: Comparison ####

datalist_comparison <- import_data("../Archive/Comparison_Paper/data/Combined/", min_counts = 2000, kingdom = "Prok", abundance_filter = F)  %>%
  filter_wrapper(min_organisms = min_organisms)

write_tsv(datalist_comparison$Count_Data, "data/Count_Data_Filtered/Comparison/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(datalist_comparison$Meta_Data, "data/Count_Data_Filtered/Comparison/Meta_Data/Prok/Meta_Data.tsv")

rm("datalist_comparison")

### Marine: SPOT ####

datalist_spot <- import_data("../SPOT/data/Count_Data/Complete/", min_counts = 2000, kingdom = "Prok", abundance_filter = F) %>%
  filter_wrapper(min_organisms = min_organisms)

write_tsv(datalist_spot$Count_Data, "data/Count_Data_Filtered/SPOT/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(datalist_spot$Meta_Data, "data/Count_Data_Filtered/SPOT/Meta_Data/Prok/Meta_Data.tsv")

rm("datalist_spot")

### Marine: CalCOFI ####

datalist_calcofi <- import_data("../CalCOFI/data/Count_Data/Complete/", min_counts = 2000, kingdom = "Prok", abundance_filter = F) %>%
  filter_wrapper(min_organisms = min_organisms)

write_tsv(datalist_calcofi$Count_Data, "data/Count_Data_Filtered/CalCOFI/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(datalist_calcofi$Meta_Data, "data/Count_Data_Filtered/CalCOFI/Meta_Data/Prok/Meta_Data.tsv")

rm("datalist_calcofi")

#### EMP datasets ####

data_secretion <- import_data("../EMP/data/Count_Data/Animal_secretion/", kingdom = "Prok", min_counts = 2000) %>%
  mutate_meta_datalist(Type = ifelse(grepl(pattern = "saliva", x = Description), "Saliva",
                                     ifelse(grepl(pattern = "Nose", x = Description), "Nose",
                                            ifelse(grepl(pattern = "oral", x = Description) | grepl(pattern = "mouth", x = Description), "Mouth", "Various")))) 
    
#### Secretion: Nose ####

data_secretion_human_nose <- data_secretion %>%
  filter_station_datalist(Type == "Nose" & host_common_name_provided == "human") %>%
  filter_wrapper(min_organisms = min_organisms)
    
data_secretion_human_nose %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_secretion_human_nose$Count_Data, "data/Count_Data_Filtered/Nose_Secretion/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_secretion_human_nose$Meta_Data, "data/Count_Data_Filtered/Nose_Secretion/Meta_Data/Prok/Meta_Data.tsv")

rm("data_secretion_human_nose")

#### Secretion: Saliva ####

data_secretion_human_saliva <- data_secretion %>%
  filter_station_datalist(Type == "Saliva" & host_common_name_provided == "human") %>%
  filter_wrapper(min_organisms = min_organisms)

data_secretion_human_saliva %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_secretion_human_saliva$Count_Data, "data/Count_Data_Filtered/Human_Saliva/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_secretion_human_saliva$Meta_Data, "data/Count_Data_Filtered/Human_Saliva/Meta_Data/Prok/Meta_Data.tsv")
    
rm("data_secretion_human_saliva", "data_secretion")

#### Distalgut: Human ####

data_distal_gut <- import_data("../EMP/data/Count_Data/Animal_distal_gut/", kingdom = "Prok", min_counts = 2000)
    
data_distal_gut_human <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "human") %>%
  filter_wrapper(min_organisms = min_organisms)

data_distal_gut_human %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_distal_gut_human$Count_Data, "data/Count_Data_Filtered/Distalgut_Human/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_distal_gut_human$Meta_Data, "data/Count_Data_Filtered/Distalgut_Human/Meta_Data/Prok/Meta_Data.tsv")

rm("data_distal_gut_human")
    
#### Distalgut: Kangaroo ####

data_distal_gut_kangaroo <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "Kangaroo") %>%
  filter_wrapper(min_organisms = min_organisms)

data_distal_gut_kangaroo %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_distal_gut_kangaroo$Count_Data, "data/Count_Data_Filtered/Distalgut_Kangaroo/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_distal_gut_kangaroo$Meta_Data, "data/Count_Data_Filtered/Distalgut_Kangaroo/Meta_Data/Prok/Meta_Data.tsv")

rm("data_distal_gut_kangaroo")

#### Distalgut: Rabbit ####
    
data_distal_gut_rabbit <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "Rabbit")  %>%
  filter_wrapper(min_organisms = min_organisms)

data_distal_gut_rabbit %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_distal_gut_rabbit$Count_Data, "data/Count_Data_Filtered/Distalgut_Rabbit/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_distal_gut_rabbit$Meta_Data, "data/Count_Data_Filtered/Distalgut_Rabbit/Meta_Data/Prok/Meta_Data.tsv")

rm("data_distal_gut_rabbit")

#### Distalgut: Deer ####

data_distal_gut_deer <- data_distal_gut %>%
  filter_station_datalist(host_common_name_provided == "Sambar deer") %>%
  filter_wrapper(min_organisms = min_organisms)

data_distal_gut_deer %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_distal_gut_deer$Count_Data, "data/Count_Data_Filtered/Distalgut_Deer/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_distal_gut_deer$Meta_Data, "data/Count_Data_Filtered/Distalgut_Deer/Meta_Data/Prok/Meta_Data.tsv")

rm("data_distal_gut_deer")

#### Distalgut: Monkey ####

data_distal_gut_monkey <- data_distal_gut %>%
      filter_station_datalist(host_common_name_provided == "Spider monkey") %>%
  filter_wrapper(min_organisms = min_organisms)

data_distal_gut_monkey %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_distal_gut_monkey$Count_Data, "data/Count_Data_Filtered/Distalgut_Monkey/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_distal_gut_monkey$Meta_Data, "data/Count_Data_Filtered/Distalgut_Monkey/Meta_Data/Prok/Meta_Data.tsv")

rm("data_distal_gut_monkey", "data_distal_gut")

#### Skin: Foot ####

data_animal_surface <- import_data("../EMP/data/Count_Data/Animal_surface/", kingdom = "Prok", min_counts = 2000) %>%
  mutate_meta_datalist(Bodypart = str_replace_all(Description, pattern = "^.*\\.", replacement = "") %>%
                         str_replace_all(pattern = "sample_[0-9]* ", replacement = ""))
    
data_surface_human_foot <- data_animal_surface %>%
  filter_station_datalist(host_common_name_provided == "human" & Bodypart == "Foot") %>%
  filter_wrapper(min_organisms = min_organisms)

data_surface_human_foot %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_surface_human_foot$Count_Data, "data/Count_Data_Filtered/Surface_Foot/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_surface_human_foot$Meta_Data, "data/Count_Data_Filtered/Surface_Foot/Meta_Data/Prok/Meta_Data.tsv")

rm("data_surface_human_foot")

#### Skin: Hand ####

data_surface_human_hand <- data_animal_surface %>%
  filter_station_datalist(host_common_name_provided == "human" & Bodypart == "Hand") %>%
  filter_wrapper(min_organisms = min_organisms)

data_surface_human_hand %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_surface_human_hand$Count_Data, "data/Count_Data_Filtered/Surface_Hand/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_surface_human_hand$Meta_Data, "data/Count_Data_Filtered/Surface_Hand/Meta_Data/Prok/Meta_Data.tsv")

rm("data_surface_human_hand", "data_animal_surface")

#### Plant: Rhizosphere ####

data_plant_rhizo_root <- import_data("../EMP/data/Count_Data/Plant_rhizosphere/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria")
    
data_plant_rhizosphere <- data_plant_rhizo_root %>%
  filter_station_datalist(Description == "Rhizosphere") %>%
  filter_wrapper(min_organisms = min_organisms)

data_plant_rhizosphere %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_plant_rhizosphere$Count_Data, "data/Count_Data_Filtered/Plant_Rhizosphere/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_plant_rhizosphere$Meta_Data, "data/Count_Data_Filtered/Plant_Rhizosphere/Meta_Data/Prok/Meta_Data.tsv")

rm("data_plant_rhizosphere")

#### Plant: Root ####

data_plant_roots <- data_plant_rhizo_root %>%
  filter_station_datalist(Description == "Roots") %>%
  filter_wrapper(min_organisms = min_organisms)

data_plant_roots %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_plant_roots$Count_Data, "data/Count_Data_Filtered/Plant_Roots/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_plant_roots$Meta_Data, "data/Count_Data_Filtered/Plant_Roots/Meta_Data/Prok/Meta_Data.tsv")

rm("data_plant_roots", "data_plant_rhizo_root")

#### Plant: Surface ####

data_plant_surface <- import_data("../EMP/data/Count_Data/Plant_surface/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria") %>%
  filter_wrapper(min_organisms = min_organisms)

data_plant_surface %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_plant_surface$Count_Data, "data/Count_Data_Filtered/Plant_Surface/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_plant_surface$Meta_Data, "data/Count_Data_Filtered/Plant_Surface/Meta_Data/Prok/Meta_Data.tsv")

rm("data_plant_surface")

#### Soil: Sandfilter ####

data_soil <- import_data("../EMP/data/Count_Data/Soil_non_saline/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria")
    
data_soil_sandfilter <- data_soil %>%
  filter_station_datalist(Description == "samples of sand from slow_sand filter water purification system" & depth_m <= 0.4) %>%
  filter_wrapper(min_organisms = min_organisms)

data_soil_sandfilter %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_soil_sandfilter$Count_Data, "data/Count_Data_Filtered/Soil_Sandfilter/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_soil_sandfilter$Meta_Data, "data/Count_Data_Filtered/Soil_Sandfilter/Meta_Data/Prok/Meta_Data.tsv")

rm("data_soil_sandfilter")

#### Soil: Cultivated ####

data_soil_cultivated_habitat <- data_soil %>%
  filter_station_datalist(study_id == 1721 & collection_timestamp == "16.11.11") %>%
  filter_wrapper(min_organisms = min_organisms)

data_soil_cultivated_habitat %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_soil_cultivated_habitat$Count_Data, "data/Count_Data_Filtered/Soil_Cultivated/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_soil_cultivated_habitat$Meta_Data, "data/Count_Data_Filtered/Soil_Cultivated/Meta_Data/Prok/Meta_Data.tsv")

rm("data_soil_cultivated_habitat")

#### Soil: Field ####

data_soil_field <- data_soil %>%
  filter_station_datalist(study_id == 990) %>%
  filter_wrapper(min_organisms = min_organisms)

data_soil_field %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_soil_field$Count_Data, "data/Count_Data_Filtered/Soil_Field/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_soil_field$Meta_Data, "data/Count_Data_Filtered/Soil_Field/Meta_Data/Prok/Meta_Data.tsv")

rm("data_soil_field")

#### Soil: Garden ####

data_soil_garden <- data_soil %>%
  filter_station_datalist(study_id == 1674) %>%
  filter_wrapper(min_organisms = min_organisms)

data_soil_garden %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_soil_garden$Count_Data, "data/Count_Data_Filtered/Soil_Garden/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_soil_garden$Meta_Data, "data/Count_Data_Filtered/Soil_Garden/Meta_Data/Prok/Meta_Data.tsv")

rm("data_soil_garden", "data_soil")

#### Indoor surface ####

data_surface_nonsaline_indoor <- import_data("../EMP/data/Count_Data/Surface_non_saline/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria") %>%
  filter_station_datalist(study_id == 2192) %>%
  filter_wrapper(min_organisms = min_organisms)

data_surface_nonsaline_indoor %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_surface_nonsaline_indoor$Count_Data, "data/Count_Data_Filtered/Surface_Nonsaline/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_surface_nonsaline_indoor$Meta_Data, "data/Count_Data_Filtered/Surface_Nonsaline/Meta_Data/Prok/Meta_Data.tsv")

rm("data_surface_nonsaline_indoor")

#### Freshwater: Timeseries ####

data_water_nonsaline <- import_data("../EMP/data/Count_Data/Water_non_saline/", kingdom = "Prok", min_counts = 2000) %>%
  filter_taxa_datalist(Class != "Chloroplast" & Family != "mitochondria")
    
data_water_nonsaline_timeseries <- data_water_nonsaline %>%
  filter_station_datalist(study_id == 1288) %>%
  filter_wrapper(min_organisms = min_organisms)

data_water_nonsaline_timeseries %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_water_nonsaline_timeseries$Count_Data, "data/Count_Data_Filtered/Water_Nonsaline_Timeseries/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_water_nonsaline_timeseries$Meta_Data, "data/Count_Data_Filtered/Water_Nonsaline_Timeseries/Meta_Data/Prok/Meta_Data.tsv")

rm("data_water_nonsaline_timeseries")

#### Freshwater: USA ####

data_water_nonsaline_usa <- data_water_nonsaline %>%
  filter_station_datalist(study_id == 1883 & env_feature == "freshwater habitat" & country == "GAZ:United States of America") %>%
  filter_wrapper(min_organisms = min_organisms)

data_water_nonsaline_usa %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_water_nonsaline_usa$Count_Data, "data/Count_Data_Filtered/Water_Nonsaline_USA/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_water_nonsaline_usa$Meta_Data, "data/Count_Data_Filtered/Water_Nonsaline_USA/Meta_Data/Prok/Meta_Data.tsv")

rm("data_water_nonsaline_usa")

#### Freshwater: Germany ####

data_water_nonsaline_germany <- data_water_nonsaline %>%
  filter_station_datalist(study_id == 945) %>%
  filter_wrapper(min_organisms = min_organisms)

data_water_nonsaline_germany %>% .$Meta_Data %>% .$doi %>% unique()

write_tsv(data_water_nonsaline_germany$Count_Data, "data/Count_Data_Filtered/Water_Nonsaline_Germany/Processed/Prok/Full_Prok_Count.tsv")
write_tsv(data_water_nonsaline_germany$Meta_Data, "data/Count_Data_Filtered/Water_Nonsaline_Germany/Meta_Data/Prok/Meta_Data.tsv")

rm("data_water_nonsaline_germany", "data_water_nonsaline")
