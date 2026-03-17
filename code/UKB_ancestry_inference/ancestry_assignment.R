rm(list=ls())
install.packages('tidyverse')
library(tidyverse)

# Load data #

## 1000G sample information ##
system("dx download file-GkqbP40JfyQxQ6vpV9zKy6Q9")
samples_1000G = read_tsv("igsr_1000_genomes_samples.tsv")

## 1000G PCs ##
system("dx download file-GkzZXf8JfyQb8yQkJ063jfZj")
PCs_1000G = read_tsv("PCs_1000G.txt")

## UKB PCs ##
system("dx download file-GkzQPP8JfyQyfFQjvYp44zj1")
PCs_UKB = read_tsv("PCs_UKB.txt")

# Prepare and check 1000G data #

## Merge 1000G sample and PC data ##
PCs_1000G = merge(samples_1000G,PCs_1000G, by = "Sample")

## Check population clustering by top PCs ##

### PC1 vs. PC2 ###
PC1_2_1000G = ggplot(data = PCs_1000G, mapping = 
                       aes(x = PC_1, y= PC_2,
                           color= as.factor(Superpopulation_code)))+
  geom_point()+
  theme_classic()

ggsave("PC1_2_1000G.jpeg", plot = PC1_2_1000G,
       dpi = 800)

### PC3 vs. PC4 ###
PC3_4_1000G = ggplot(data = PCs_1000G, mapping = 
                       aes(x = PC_3, y= PC_4,
                           color= as.factor(Superpopulation_code)))+
  geom_point()+
  theme_classic()

ggsave("PC3_4_1000G.jpeg", plot = PC3_4_1000G,
       dpi = 800)

# Combine UKB and 1000G data
colnames(PCs_UKB) = c("Sample",paste0("PC_",1:30))
PCs_UKB = PCs_UKB %>% mutate(Sex = NA,
                             ID = NA,
                             Population_code = "UKB",
                             Population_name = "UK Biobank",
                             Superpopulation_code = NA,
                             Superpopulation_name = NA,
                             Population_elastic_ID = NA,
                             Data_collections = NA)
PCs_UKB = PCs_UKB %>% select(names(PCs_1000G))

PCA = rbind.data.frame(PCs_1000G,PCs_UKB)

# Assign 1000G populations without an unassigned group #
lst <- split(PCA[,10:39], PCA$Population_code)
lst <- lst[!names(lst) %in% "UKB"]

mds <- do.call(cbind, lapply(lst, function(x) {
  md <- mahalanobis(x, colMeans(x), cov = cov(x))
  md_sd <- (md - mean(md))/sd(md)
  x <- x[md_sd < 2,]
  return(mahalanobis(PCA[,10:39], colMeans(x), cov(x)))
}
))

PCA$Assigned_Population <- c(names(lst))[apply(mds, 1, which.min)]

mapper <- setNames(rep(c("EAS","EUR","AFR","AMR","SAS"), c(5,5,7,4,5)),
                   c("CHB", "JPT", "CHS", "CDX", "KHV","CEU", "TSI", "FIN", "GBR", "IBS",
                     "LWK", "YRI", "GWD", "MSL", "ESN", "ASW", "ACB","MXL", "PUR", "CLM", "PEL",
                     "GIH","PJL","BEB","STU","ITU"))

PCA$Assigned_Superpopulation <-mapper[PCA$Assigned_Population]

# Assign 1000G populations with an unassigned group #

# Function to calculate Mahalanobis distances and thresholds
calculate_distances_and_thresholds <- function(data, group_name) {
  # Calculate Mahalanobis distances for the group
  md <- mahalanobis(data, colMeans(data), cov = cov(data))
  md_sd <- (md - mean(md)) / sd(md)
  # Filter outliers based on standard deviation threshold
  filtered_data <- data[md_sd < 2,]
  
  # Calculate Mahalanobis distances for the PCA data
  distances <- mahalanobis(PCA[,10:39], colMeans(filtered_data), cov(filtered_data))
  
  # Calculate the 95th percentile threshold for the group's distances
  threshold <- quantile(distances, probs = 0.95, na.rm = TRUE)
  
  return(list(distances = distances, threshold = threshold))
}

# Calculate Mahalanobis distances and thresholds for each group
results <- lapply(lst, function(x) calculate_distances_and_thresholds(x, names(lst)[which(lst == x)]))
names(results) <- names(lst)

# Extract distances and thresholds for each group
distances_list <- lapply(results, `[[`, "distances")
thresholds_list <- sapply(results, `[[`, "threshold")

# Determine the minimum Mahalanobis distance and corresponding thresholds for each individual
min_distances <- apply(do.call(cbind, distances_list), 1, min, na.rm = TRUE)

# Assign groups or "Unassigned" based on thresholds
PCA$Assigned_Population_2 <- apply(do.call(cbind, distances_list), 1, function(row) {
  # Check which group has the minimum distance
  group_with_min_distance <- which.min(row)
  
  # Check if this distance is above the threshold for all groups
  if (row[group_with_min_distance] > thresholds_list[group_with_min_distance]) {
    return("Unassigned")
  } else {
    return(names(lst)[group_with_min_distance])
  }
})

mapper <- setNames(rep(c("EAS","EUR","AFR","AMR","SAS","Unassigned"), c(5,5,7,4,5,1)),
                   c("CHB", "JPT", "CHS", "CDX", "KHV","CEU", "TSI", "FIN", "GBR", "IBS",
                     "LWK", "YRI", "GWD", "MSL", "ESN", "ASW", "ACB","MXL", "PUR", "CLM", "PEL",
                     "GIH","PJL","BEB","STU","ITU","Unassigned"))

PCA$Assigned_Superpopulation_2 <-mapper[PCA$Assigned_Population_2]

# Prepare the final dataset #
`%!in%` = Negate(`%in%`)

final_assignments = PCA %>%
  filter(Sample %!in% samples_1000G$Sample) %>% 
  select(c(Sample,Assigned_Population, Assigned_Population_2,
           Assigned_Superpopulation, Assigned_Superpopulation_2))
write.table(final_assignments, file = "UKB_ancestry_assignments.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Save results #
system("dx upload /home/rstudio-server/ancestry_assignment.R")
system("dx upload /home/rstudio-server/UKB_ancestry_assignments.txt")
