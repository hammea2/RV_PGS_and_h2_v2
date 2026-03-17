# Setup ##
rm(list = ls())
install.packages("tidyverse")
install.packages("data.table")
install.packages("R.utils")
library(tidyverse)
library(data.table)
library(R.utils)

## Load 1000G bim ##

# 1000G #
system("dx download file-GkqYB5jJfyQyzBFkgv28yp2p")
system("dx download file-GkqZ8P0JfyQxx78Q6VB68068")
system("dx download file-GkqZY90JfyQv68kF7xv1KZF0")
system("dx download file-GkqZgQ8JfyQxx78Q6VB6810V")
system("dx download file-GkqZqz0JfyQgyfq91K3BQkPF")
system("dx download file-GkqZzv8JfyQkG57GJbgFYy4Z")
system("dx download file-Gkqb41QJfyQzQFjZBjgK02p0")
system("dx download file-Gkqb7y0JfyQp44BV3q9Fgf1j")
system("dx download file-GkqbF6QJfyQpF6P9ZQZ2k4yB")
system("dx download file-GkqYKq0JfyQjKyxQzkXJZBf1")
system("dx download file-GkqYVvQJfyQXpP5pv7G0FJqb")
system("dx download file-GkqYZj0JfyQzkGg598G5K1K7")
system("dx download file-GkqYgZ0JfyQX64K1BJ40zV2y")
system("dx download file-GkqYq00JfyQjKyxQzkXJZGGX")
system("dx download file-GkqYx20JfyQgfF30bZQPXfyK")
system("dx download file-GkqYzZQJfyQX64K1BJ40zX3Q")
system("dx download file-GkqZ278JfyQgyfq91K3BQfqJ")
system("dx download file-GkqZ4GjJfyQvp1f0F053jfp8")
system("dx download file-GkqZ710JfyQxx78Q6VB6801G")
system("dx download file-GkqZJY8JfyQQ2jbKjvyP92yF")
system("dx download file-GkqZP8QJfyQQ2jbKjvyP934q")
system("dx download file-GkqZQFQJfyQp3QvBxVB80yz9")

# UKB #
system("dx download file-FxXZxjQJkF636yj326QqXPV5")
system("dx download file-FxXZxp8JkF67GV2z7v5ZqFzZ")
system("dx download file-FxXZxpQJkF69vjv312xj39xz")
system("dx download file-FxXZxq0JkF65bZBj13QpBZFZ")
system("dx download file-FxXZxq8JkF61JXY78Y2xp2Z6")
system("dx download file-FxXZxqjJkF6PpqqK0p2gZxvb")
system("dx download file-FxXZxv8JkF6457pxPkfgJ60Y")
system("dx download file-FxXZxvQJkF67Z48f0pjq53Q8")
system("dx download file-FxXZxx0JkF6BXzY51BBbz1f9")
system("dx download file-FxXZxZjJkF64ZZKv1B2V667F")
system("dx download file-FxXZxb0JkF69X2FQ19v6jXVQ")
system("dx download file-FxXZxbQJkF6G2F900kkp0YY4")
system("dx download file-FxXZxbjJkF6GbX1619P5YBj9")
system("dx download file-FxXZxf8JkF6BXzY51BBbz1bG")
system("dx download file-FxXZxfjJkF67Z48f0pjq53PG")
system("dx download file-FxXZxg0JkF68b43G0pf356Y0")
system("dx download file-FxXZxgQJkF66qV4X0p7KBgJG")
system("dx download file-FxXZxgjJkF6B6gJg1B6qXGjf")
system("dx download file-FxXZxj8JkF6B0kfB44f4PpVK")
system("dx download file-FxXZxk0JkF6GbKYb9VP89bYB")
system("dx download file-FxXZxk8JkF67GV2z7v5ZqFzQ")
system("dx download file-FxXZxkjJkF641JG196p7ggbF")

# High LD regions #
system("dx download file-Gkqg9y0JfyQvJqFzQfjPgzBx")
LD_regions  = fread("high_LD_regions_hg19_GRCh37.txt")

## Identify SNPs in common ##

for (chr in 1:22) {
  bim_UKB = fread(paste0("ukb22418_c",chr,"_b0_v2.bim"))
  
  bim_1000G = fread(paste0("chr",chr,"_1000G.bim"))
  
  common_UKB = bim_UKB %>% filter(V2 %in% bim_1000G$V2)
  
  ld_chr = LD_regions %>% filter(V1 == paste0("chr",chr))
  
  common_UKB$exclude = NA
  
  if(nrow(ld_chr)>0){
    for (i in 1:nrow(common_UKB)) {
      for (j in 1:nrow(ld_chr)) {
        common_UKB$exclude[i] = ifelse((common_UKB$V4[i] >= ld_chr$V2[j] & common_UKB$V4[i] <= ld_chr$V3[j])==TRUE,1,0)
        
      }
    }
    
    common_UKB = common_UKB %>% filter(exclude == 0)
    
  }
  
  common_UKB = common_UKB %>% select(V2)
  
  fwrite(common_UKB,paste0("chr",chr,"_UKB_common.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}


# Save the script, phenotype, covariate, and sample filtering files #
system("dx upload /home/rstudio-server/identify_variants_in_common.R")

for (i in 1:22) {
  system(paste0("dx upload /home/rstudio-server/chr",i,"_UKB_common.txt"))
}

