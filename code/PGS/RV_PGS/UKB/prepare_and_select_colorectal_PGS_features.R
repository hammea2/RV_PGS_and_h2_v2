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
for (i in c("small","medium_small","medium_large","large")) {
  dir = file.path("/lindstroem/austin_working/Dissertation/Aim_3/Rare_variant_PGS/UKB_trained/colorectal_results",i) 
  files = get_files(dir)
  for (j in files) {
    j_data = get_results_from_STAAR(j)
    data = rbind.data.frame(data,j_data)
  }
}

data$`STAAR-B` = as.numeric(data$`STAAR-B`)
data$`#SNV` = as.numeric(data$`#SNV`)

# Select burden features # 

## Known colorectal cancer predisposition genes ##

known_colorectal_genes = c(
  "APC", "BLM", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "EGFR", "EPCAM",
  "FANCM", "KRAS", "MLH1", "MSH2", "MSH6", "MUTYH", "NTHL1", "PIK3CA",
  "PMS2", "POLD1", "POLE", "PTEN", "SMAD2", "SMAD4", "STK11",
  "TGFBR1", "TGFBR2", "TP53", "VHL"
)

category_priority <- c(
  "plof",
  "plof_ds",
  "disruptive_missense",
  "missense",
  "synonymous"
)

selected_features <- data %>%
  filter(`STAAR-B` < 0.001 | (`Gene name` %in% known_colorectal_genes & `STAAR-B` < 0.1)) %>%
  mutate(
    Category_rank = match(Category, category_priority)
  ) %>%
  group_by(`Gene name`,Chr,`#SNV`) %>%
  arrange(Category_rank, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-Category_rank)

# Write RDS with selected burden features #

write_rds(x = selected_features, file = "/lindstroem/austin_working/Dissertation/Aim_3/Rare_variant_PGS/UKB_trained/selected_features/colorectal_selected_features_expanded.rds")

