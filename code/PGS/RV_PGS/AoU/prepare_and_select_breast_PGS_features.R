# Setup #
library(pacman)
p_load(tidyverse)

# Functions #

get_results_from_STAAR = function(x){
  list = readRDS(x)
  data = NULL
  
  for(i in 1:length(list)){
    list_data = NULL
    list_data <-  rbind(list_data,unlist(list[[i]]$plof[,c(1:4,10)]))
    list_data <-  rbind(list_data,unlist(list[[i]]$plof_ds[,c(1:4,10)]))
    list_data <-  rbind(list_data,unlist(list[[i]]$missense[,c(1:4,10)]))
    list_data <-  rbind(list_data,unlist(list[[i]]$disruptive_missense[,c(1:4,10)]))
    list_data <-  rbind(list_data,unlist(list[[i]]$synonymous[,c(1:4,10)]))
    data = rbind.data.frame(data,list_data)
  }
  
  return(data)
}

get_files <- function(dir = ".") {
  pattern <- "^chr[0-9]+_min[0-9]+_max[0-9]+_STAAR_pipeline_results_full\\.rds$"
  list.files(dir, pattern = pattern, full.names = TRUE)
}

# Load STAAR pipeline results #
data = NULL
   dir = file.path("/lindstroem/austin_working/Dissertation/Aim_3/Rare_variant_PGS/AoU_trained/breast_results") 
   files = get_files(dir)
   for (j in files) {
     j_data = get_results_from_STAAR(j)
     data = rbind.data.frame(data,j_data)
   }

data$`STAAR-B` = as.numeric(data$`STAAR-B`)
data$`#SNV` = as.numeric(data$`#SNV`)

# Select burden features # 

## Known breast cancer predisposition genes ##

known_breast_genes = c(
  "ABRAXAS1", "APC", "ATM", "ATRIP", "BARD1", "BLM", "BRCA1", "BRCA2",
  "BRIP1", "CDH1", "CHEK2", "ERBB2", "FANCC", "FANCD2", "FANCM", "HOXB13",
  "MCPH1", "MEN1", "MLH1", "MRE11", "MSH2", "MSH6", "MUTYH", "NBN", "NF1",
  "PALB2", "PMS2", "POLG", "PPM1D", "PTEN", "RAD50", "RAD51B", "RAD51C",
  "RAD51D", "RBBP8", "RECQL", "RINT1", "SERPINA3", "STK11", "TEX15",
  "TP53", "XRCC2"
)


category_priority <- c(
  "plof",
  "plof_ds",
  "disruptive_missense",
  "missense",
  "synonymous"
)

selected_features <- data %>%
  filter(`STAAR-B` < 0.001 | (`Gene name` %in% known_breast_genes & `STAAR-B` < 0.1)) %>%
  mutate(
    Category_rank = match(Category, category_priority)
  ) %>%
  group_by(`Gene name`,Chr,`#SNV`) %>%
  arrange(Category_rank, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-Category_rank)

# Write RDS with selected burden features #

write_rds(x = selected_features, file = "/lindstroem/austin_working/Dissertation/Aim_3/Rare_variant_PGS/AoU_trained/selected_features/breast_selected_features_expanded.rds")
