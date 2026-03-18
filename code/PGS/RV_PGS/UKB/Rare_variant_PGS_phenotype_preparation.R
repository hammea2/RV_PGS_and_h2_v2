# Setup #
install.packages("pacman")
library(pacman)

p_load(rlang, tidyverse, data.table,
       R.utils)

system("dx download file-J4Jvbq8JfyQj6pz5KBQ6KBP3")
system("dx download file-J4jzYQ8Jpjk8Jp1599p562K6")
system("dx download file-J4jzYQ8Jpjk25xGB43yfxj8G")
system("dx download file-J4jzYQ8Jpjk5412kPffqfVyZ")

# Functions #
normalize_pgs_pca <- function(data, pgs_col, pca_cols, ref_filter = NULL) {
  # --- 1. Define reference population ---
  if (!is.null(ref_filter)) {
    ref_data <- data %>% filter(!!parse_expr(ref_filter))
  } else {
    ref_data <- data
  }
  
  # --- 2. Fit regression model for mean adjustment (Z_norm1) ---
  formula_mean <- as.formula(paste(pgs_col, "~", paste(pca_cols, collapse = " + ")))
  mean_model <- lm(formula_mean, data = ref_data)
  
  # Predicted mean PGS for all individuals
  data$pred_mean <- predict(mean_model, newdata = data)
  
  # Residuals
  data$resid_mean <- data[[pgs_col]] - data$pred_mean
  
  # SD of residuals in reference population
  if (!is.null(ref_filter)) {
    sd_ref <- sd(data$resid_mean[with(data, eval(parse(text = ref_filter)))], na.rm = TRUE)
  } else {
    sd_ref <- sd(data$resid_mean, na.rm = TRUE)
  }
  
  # Z_norm1: mean-adjusted
  data$Z_norm1 <- data$resid_mean / sd_ref
  
  # --- 3. Fit regression model for variance adjustment (Z_norm2) ---
  # Create ref_data with residuals for variance model
  ref_data <- if (!is.null(ref_filter)) {
    data %>% filter(!!parse_expr(ref_filter))
  } else {
    data
  }
  
  ref_data <- ref_data %>% mutate(resid_sq = resid_mean^2)
  formula_sd <- as.formula(paste("resid_sq ~", paste(pca_cols, collapse = " + ")))
  sd_model <- lm(formula_sd, data = ref_data)
  
  # Predicted SD for all individuals
  data$pred_var <- predict(sd_model, newdata = data)
  data$pred_sd <- sqrt(pmax(data$pred_var, 1e-8))
  
  # Z_norm2: mean + variance normalized
  data$Z_norm2 <- data$resid_mean / data$pred_sd
  
  return(data)
}


# Load phenotypes and prepare PCA-normalized common variant PGS scores
pheno = fread("aim_3_phenotypes.txt")
breast_pheno = pheno %>%
  filter(breast_cancer > -1)
colorectal_pheno = pheno %>%
  filter(colorectal_cancer > -1)
prostate_pheno = pheno %>%
  filter(prostate_cancer > -1)

breast_score = fread("breast_scores.sscore") %>%
  select(IID,SCORE1_AVG)
colnames(breast_score) = c("Participant_ID","breast_PGS")
colorectal_score = fread("colorectal_scores.sscore")%>%
  select(IID,SCORE1_AVG)
colnames(colorectal_score) = c("Participant_ID","colorectal_PGS")
prostate_score = fread("prostate_scores.sscore")%>%
  select(IID,SCORE1_AVG)
colnames(prostate_score) = c("Participant_ID","prostate_PGS")

breast_pheno = merge(breast_pheno,breast_score)
colorectal_pheno = merge(colorectal_pheno,colorectal_score)
prostate_pheno = merge(prostate_pheno,prostate_score)

breast_pheno <- normalize_pgs_pca(
  breast_pheno,
  pgs_col = "breast_PGS",
  pca_cols = paste0("UKB_PC", 1:20),
  ref_filter = "breast_cancer == 0"
)

colorectal_pheno <- normalize_pgs_pca(
  colorectal_pheno,
  pgs_col = "colorectal_PGS",
  pca_cols = paste0("UKB_PC", 1:20),
  ref_filter = "colorectal_cancer == 0"
)

prostate_pheno <- normalize_pgs_pca(
  prostate_pheno,
  pgs_col = "prostate_PGS",
  pca_cols = paste0("UKB_PC", 1:20),
  ref_filter = "prostate_cancer == 0"
)

# Prepare data for STAARpipeline test #

STAARpipeline_test_prostate_training_data = prostate_pheno %>%
  filter(Training_for_workflow_tests == "Training_for_workflow_tests") %>% 
  select(c(Participant_ID,prostate_cancer,Age_at_recruitment,Z_norm2,paste0("UKB_PC",1:20)))

colnames(STAARpipeline_test_prostate_training_data) = c("Participant_ID","Prostate_cancer","Age","Standardized_CV_PRS",paste0("UKB_PC",1:20))

write.csv(STAARpipeline_test_prostate_training_data,"STAAR_pipeline_test_prostate_cancer_phenotypes.csv", row.names = FALSE, quote = FALSE)

system("dx upload STAAR_pipeline_test_prostate_cancer_phenotypes.csv --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/phenotypes/")

# Prepare data for full STAARpipeline analyses #

## Breast ##

STAARpipeline_breast_training_data = breast_pheno %>%
  filter(internal_set == "Training") %>% 
  select(c(Participant_ID,breast_cancer,Age_at_recruitment,Z_norm2,paste0("UKB_PC",1:20)))

colnames(STAARpipeline_breast_training_data) = c("Participant_ID","breast_cancer","Age","Standardized_CV_PRS",paste0("UKB_PC",1:20))

write.csv(STAARpipeline_breast_training_data,"STAAR_pipeline_breast_cancer_phenotypes.csv", row.names = FALSE, quote = FALSE)

system("dx upload STAAR_pipeline_breast_cancer_phenotypes.csv --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/phenotypes/")

## Colorectal ##

STAARpipeline_colorectal_training_data = colorectal_pheno %>%
  filter(internal_set == "Training") %>% 
  select(c(Participant_ID,colorectal_cancer,Age_at_recruitment,Female,Z_norm2,paste0("UKB_PC",1:20)))

colnames(STAARpipeline_colorectal_training_data) = c("Participant_ID","colorectal_cancer","Age","Female","Standardized_CV_PRS",paste0("UKB_PC",1:20))

write.csv(STAARpipeline_colorectal_training_data,"STAAR_pipeline_colorectal_cancer_phenotypes.csv", row.names = FALSE, quote = FALSE)

system("dx upload STAAR_pipeline_colorectal_cancer_phenotypes.csv --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/phenotypes/")


## Prostate ##

STAARpipeline_prostate_training_data = prostate_pheno %>%
  filter(internal_set == "Training") %>% 
  select(c(Participant_ID,prostate_cancer,Age_at_recruitment,Z_norm2,paste0("UKB_PC",1:20)))

colnames(STAARpipeline_prostate_training_data) = c("Participant_ID","Prostate_cancer","Age","Standardized_CV_PRS",paste0("UKB_PC",1:20))

write.csv(STAARpipeline_prostate_training_data,"STAAR_pipeline_prostate_cancer_phenotypes.csv", row.names = FALSE, quote = FALSE)

system("dx upload STAAR_pipeline_prostate_cancer_phenotypes.csv --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/phenotypes/")


system("dx upload Rare_variant_PGS_phenotype_preparation.R --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/scripts/Rare_variant_PGS/phenotype_preparation/")
