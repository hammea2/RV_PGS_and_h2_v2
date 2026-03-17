# Setup ##
rm(list = ls())
install.packages("tidyverse")
install.packages("data.table")
install.packages("R.utils")
library(tidyverse)
library(data.table)
library(R.utils)

## Load data ##

# 1000G .bim #
system("dx download file-Gkv7280JXvPV6fPQGx1zK0K9")

# UKB .bim #
system("dx download file-Gkqv38jJXqpvv7gv6KJjy9Yb")

## Identify SNPs that need to be flipped ##

bim_UKB = fread("UKB_1000G_harmonized_pruned.bim")

bim_1000G = fread("1000G_UKB_harmonized.bim")

SNPs_to_flip = NULL

for (i in 1:nrow(bim_1000G)) {
  if(bim_1000G[i,'V5']!=bim_UKB[i,'V5']){
    SNPs_to_flip = c(SNPs_to_flip,bim_1000G[i,'V2'])
  }
}

fwrite(as.data.frame(SNPs_to_flip),
       file = "1000G_SNPs_to_flip.txt", row.names = FALSE,
       col.names = FALSE, quote = FALSE)

# Save the script and other files #
system("dx upload /home/rstudio-server/1000G_variants_to_flip.R")
system("dx upload /home/rstudio-server/1000G_SNPs_to_flip.txt")
