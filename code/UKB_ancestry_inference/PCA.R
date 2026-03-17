# Setup ##
rm(list = ls())
install.packages("tidyverse")
install.packages("data.table")
install.packages("R.utils")
install.packages("devtools")
devtools::install_github("gabraham/flashpca/flashpcaR")
library(flashpcaR)
library(tidyverse)

## Load data ##

# 1000G #
system("dx download file-Gkv9p9QJpf03bYzyJ0bxFkgV")
system("dx download file-Gkv9p9QJpf0219G0fxzVBp02")
system("dx download file-Gkv9p9QJpf0BbYBZjVxXq1Q4")

# UKB #
system("dx download file-Gkv9p9QJpf03xgZq8V6b11ZB")
system("dx download file-Gkv9pBjJpf02f3Vz0jXxvGGg")
system("dx download file-Gkv9pB8Jpf0JyXb994zq4Kv0")

# PCA #

## Run FlashPCA on 1000G data ##
PCA_1000G = flashpca("1000G_final", ndim = 30,
                     do_loadings = TRUE,
                     return_scale = TRUE)

## Project UKB samples into 1000G PCA space ##
bim_1000G = read.table("1000G_final.bim")
bim_UKB = read.table("UKB_final.bim")

PCA_UKB = flashpcaR::project(X = "UKB_final",
                             loadings = PCA_1000G$loadings,
                             orig_mean = PCA_1000G$center,
                             orig_sd = PCA_1000G$scale,
                             ref_alleles=bim_1000G[,5],
                             verbose = TRUE)

# Prepare PC data sets to save ##

## 1000G ##
PCs_1000G = PCA_1000G$projection %>% as.data.frame()
PCs_1000G = PCs_1000G %>% rownames_to_column(var = "ID")
PCs_1000G = PCs_1000G %>% mutate(ID = sub(":.*", "", ID))
colnames(PCs_1000G) = c("Sample",paste0("PC_",1:30))
write.table(PCs_1000G,file = "PCs_1000G.txt", quote = FALSE,
            row.names = FALSE, col.names = TRUE,
            sep = "\t")
system("dx upload /home/rstudio-server/PCs_1000G.txt")

## UKB ##
PCs_UKB = PCA_UKB$projection %>% as.data.frame()
PCs_UKB = PCs_UKB %>% rownames_to_column(var = "ID")
PCs_UKB = PCs_UKB %>% mutate(ID = sub(":.*", "", ID))
PCs_UKB = PCs_UKB %>% filter(ID>0)
colnames(PCs_UKB) = c("eid",paste0("PC_",1:30))

write.table(PCs_UKB,file = "PCs_UKB.txt", quote = FALSE,
            row.names = FALSE, col.names = TRUE,
            sep = "\t")
system("dx upload /home/rstudio-server/PCs_UKB.txt")

# Save the script and other files #
system("dx upload /home/rstudio-server/PCA.R")
