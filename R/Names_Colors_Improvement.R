## This script is only for improving plotting results by improving names of datasets and
# choosing nice colors for each dataset.
# 
# Author: Felix Milke
# Date: 26.02.2025

load(file = "output/clustering_datasets_04.RData")

new_dataset_names <- tibble(Type = purrr::map_chr(clustering_datasets_0.4, "Type"),
                            New_Type = c("Marine: CalCOFI", "Marine: Lat. transect", "Distalgut: Deer (EMP)", "Distalgut: Human (EMP)",
                                         "Distalgut: Kangaroo (EMP)", "Distalgut: Monkey (EMP)", "Distalgut: Rabbit (EMP)",
                                         "Secretion: Human saliva (EMP)", "Freshwater: Lake Mendota", "Secretion: Human nose (EMP)",
                                         "Plant: Rhizosphere (EMP)", "Plant: Roots (EMP)", "Plant: Surface (EMP)", 
                                         "Soil: Australia", "Soil: Cultivated (EMP)", "Soil: Field (EMP)",
                                         "Soil: Garden (EMP)", "Soil: Sandfilter (EMP)", "Marine: SPOT", "Skin: Human foot (EMP)",
                                         "Skin: Human hand (EMP)", "Indoor surface (EMP)", "Freshwater: Germany (EMP)", "Freshwater: Timeseries (EMP)",
                                         "Freshwater: USA (EMP)"))

biome_type <- tibble(
  Type = c("Distalgut_Deer","Distalgut_Human","Distalgut_Kangaroo","Distalgut_Monkey","Distalgut_Rabbit","Human_Saliva",
           "Lake_Mendota","Nose_Secretion","Plant_Rhizosphere","Plant_Roots","Plant_Surface",
           "Soil_Australia","Soil_Cultivated","Soil_Field","Soil_Garden","Soil_Sandfilter","Surface_Foot","Surface_Hand",
           "Surface_Nonsaline","Water_Nonsaline_Germany","Water_Nonsaline_Timeseries","Water_Nonsaline_USA", "SPOT", "Comparison", "CalCOFI"),
  Biome = c("Distalgut", "Distalgut", "Distalgut", "Distalgut", "Distalgut", "Human Secretion",
            "Freshwater", "Human Secretion", "Plant-associated", "Plant-associated", "Plant-associated",
            "Soil", "Soil", "Soil", "Soil", "Soil", "Human skin", "Human skin", 
            "Indoor surface", "Freshwater", "Freshwater", "Freshwater", "Marine", "Marine", "Marine"),
  New_Type = c("Distalgut: Deer (EMP)", "Distalgut: Human (EMP)",
               "Distalgut: Kangaroo (EMP)", "Distalgut: Monkey (EMP)", "Distalgut: Rabbit (EMP)", 
               "Secretion: Human saliva (EMP)", "Freshwater: Lake Mendota", "Secretion: Human nose (EMP)",
               "Plant: Rhizosphere (EMP)", "Plant: Roots (EMP)", "Plant: Surface (EMP)", 
               "Soil: Australia", "Soil: Cultivated (EMP)", "Soil: Field (EMP)",
               "Soil: Garden (EMP)", "Soil: Sandfilter (EMP)","Skin: Human foot (EMP)",
               "Skin: Human hand (EMP)", "Indoor surface (EMP)", "Freshwater: Germany (EMP)", "Freshwater: Timeseries (EMP)",
               "Freshwater: USA (EMP)", "Marine: SPOT", "Marine: Lat. transect", "Marine: CalCOFI"),
  Dataset = c("EMP", "EMP", "EMP", "EMP", "EMP", "EMP", "Lake_Mendota", "EMP", "EMP", "EMP", "EMP",
              "Soil_Australia", "EMP", "EMP", "EMP", "EMP", "EMP", "EMP", "EMP", "EMP", "EMP", "EMP", "SPOT", "Comparison", "CalCOFI"),
  Host = c("Deer", "Human", "Kangaroo", "Monkey", "Rabbit", "Human", "Lake", "Human", "Plant", "Plant", "Plant", 
           "Soil", "Soil", "Soil", "Soil", "Soil", "Human", "Human", "Indoor", "Lake", "Lake", "Lake", "Ocean", "Ocean", "Ocean")
)

# Colorblind-palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

distsalgut <- c("#3a30a8", "#9E56CF", "#56c4cf", "#cf5687", "#87cf56")
freshwater <- c("#009E73", "#F0E442", "#0072B2", "#D55E00")
human_secretion <- c("#ffc857", "#2d6a4f")
human_skin <- c("#023e7d", "#b42121")
indoor_surface <- "#8fce00"
marine <- c("#E69F00", "#bae600", "#e62c00", "#56B4E9")
plant_associated <- c("#6a329f", "#9f6a32", "#329f6a")
soil <- c("#0050ce", "#ce7e00",	"#b7ce00",	"#ce1700", "#7e00ce")

comb_palette <- c(distsalgut, freshwater, human_secretion, human_skin, indoor_surface, marine, plant_associated, soil)

Palette3 <- c("#f44336", "#8fce00", "#6a329f", "#ce7e00", "#2986cc")
