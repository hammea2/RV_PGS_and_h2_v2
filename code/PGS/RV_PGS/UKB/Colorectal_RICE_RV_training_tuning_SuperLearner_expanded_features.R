# Setup #
install.packages("pacman")
library(pacman)

p_load(rlang, RISCA, gdsfmt, SeqArray,
       SeqVarTools, dplyr,
       Matrix, glmnet,
       caret, GenomicFeatures, 
       TxDb.Hsapiens.UCSC.hg38.knownGene,
       BiocManager, tidyverse, data.table,
       R.utils, ranger, SuperLearner)

BiocManager::install("GENESIS")
library(GENESIS)

p_load_gh("xihaoli/STAAR",
          "xihaoli/MultiSTAAR",
          "zilinli1988/SCANG",
          "xihaoli/STAARpipeline")

p_load_gh('jwilliams10/RICE')

# Download aGDS files #
system("dx download -r project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/aGDS/")

# Download phenotype file #
system("dx download file-J4Jvbq8JfyQj6pz5KBQ6KBP3")

# Download selected features file #
system("dx download file-J5PQ6f8JfyQXj0Zz04vzbpJk")

# Load phenotype data #
pheno = data.table::fread("aim_3_phenotypes.txt")

colorectal_pheno = pheno %>%
  dplyr::filter(colorectal_cancer > -1) %>%
  dplyr::select(c(Participant_ID,Age_at_recruitment,Female,paste0("UKB_PC",1:20),
                  Superpopulation,internal_set,colorectal_cancer))

rm(pheno)

# Load selected features #

selected_features <- readRDS("colorectal_selected_features_expanded.rds")
colnames(selected_features) <- c("Gene","Chr","Category","N_SNV","STAAR_B_p_value")

# Define paramaters to extract burden scores from GDS files #

## Define script variables ##
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category",
          "MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF",
          "aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info",
         "/genecode_comprehensive_exonic_category","/metasvm_pred",
         "/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed","/apc_epigenetics_transcription",
         "/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")
Annotation_name_catalog <- data.frame(name=name,dir=dir)
QC_label <- "annotation/filter"
variant_type <- "variant"
geno_missing_imputation <- "mean"
Annotation_dir <- "annotation/info/FunctionalAnnotation"
Use_annotation_weights <- TRUE
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

# Generate burden scores for selected starting features #
Sig_Burdens <- NULL
for(chr in c(1:22)){
  print(paste0("Processing chromosome ", chr))
  sig_results = selected_features %>%
    filter(Chr == chr)
  if(nrow(sig_results) == 0){
    next
  } 
  genofile <- seqOpen(paste0("./aGDS/chr",chr,".gds"))
  Chr_Sig_Burdens <- NULL
  for(i in 1:nrow(sig_results)){
    max = nrow(sig_results)
    chr <- sig_results$Chr[i]
    gene_name <- sig_results$Gene[i]
    category <- sig_results$Category[i]
    print(paste0("[",i,"/", max,"] Processing gene ", gene_name, " with annotation mask ", category))
    a <- Burden_Scores(region = "Coding",chr = chr, gene_name = gene_name,category = category,
                       genofile = genofile,IDs = colorectal_pheno$Participant_ID,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,silent = TRUE) 
    
    ### Combine them column-wise
    Chr_Sig_Burdens <- cbind(Chr_Sig_Burdens,a)
  }
  Sig_Burdens <- cbind(Sig_Burdens,Chr_Sig_Burdens)
  rm(Chr_Sig_Burdens)
  seqClose(genofile)
}

# Merge burden scores with phenotypes and align #
X <- data.frame(Participant_ID = colorectal_pheno$Participant_ID,Sig_Burdens)
Phenotype_Data <- inner_join(colorectal_pheno, X, by = "Participant_ID") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])

train_set = Phenotype_Data %>%
  filter(internal_set == "Training") %>%
  pull(Participant_ID)

validation_set = Phenotype_Data %>%
  filter(internal_set %in% c("Test","Tuning")) %>%
  pull(Participant_ID)

X_train <- X_Aligned[Phenotype_Data$Participant_ID %in% train_set,] %>% as.data.frame()

# Train SuperLearner model #
pheno_train = colorectal_pheno %>%
  filter(internal_set == "Training") %>% 
  as.data.frame()
trait = "colorectal_cancer"
Y <- as.numeric(pheno_train[, trait])

set.seed(1814879)

V <- 5 
folds <- createFolds(factor(Y), k = V, list = TRUE, returnTrain = FALSE)


SL.glmnet_0 <- function(...) {
  SL.glmnet(..., alpha = 0)
}

SL.glmnet_1 <- function(...) {
  SL.glmnet(..., alpha = 1)
}

SL.logit_glm <- function(Y, X, newX, family, obsWeights, ...) {
  
  dat <- data.frame(Y = Y, X)
  
  fit <- tryCatch(
    glm(
      formula = Y ~ .,
      family = binomial(link = "logit"),
      data = dat,
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    pred <- rep(NA_real_, nrow(newX))
  } else {
    pred <- tryCatch(
      as.numeric(predict(fit, newdata = newX, type = "response")),
      error = function(e) rep(NA_real_, nrow(newX))
    )
  }
  
  list(pred = pred, fit = fit)
}

SL.library <- c(
  "SL.glmnet_0",
  "SL.glmnet",
  "SL.glmnet_1",
  "SL.logit_glm"
)

sl_fit <- SuperLearner(
  Y = Y,
  X = X_train,
  family = binomial(),
  SL.library = SL.library,
  method = "method.AUC",     
  cvControl = list(
    V = V,
    validRows = folds         
  ),
  verbose = TRUE
)

# Check ensemble model contributions from each base learner #
sl_fit$coef
sl_weights <- sort(sl_fit$coef, decreasing = TRUE)
print(sl_weights)

# Estimate feature coefficients from ensemble model predictions using linear regression #
y_train_ridge <- sl_fit$SL.predict
X_train_mat <- as.matrix(X_train)           

df_train_eval <- data.frame(
  Y = pheno_train[, trait],
  SL_pred = y_train_ridge
)

lambda_grid = 10^seq(3, -3, length.out = 50) 

cv_ridge_pos <- cv.glmnet(
  X_train_mat,
  y_train_ridge,
  alpha = 0,
  lambda = lambda_grid,
  nfolds = 10,
  standardize = FALSE,
  lower.limits = 0
)

best_lambda <- cv_ridge_pos$lambda.min   
ridge_coef <- as.numeric(coef(cv_ridge_pos, s = best_lambda))[-1] 
ridge_intercept <- as.numeric(coef(cv_ridge_pos, s = best_lambda))[1]

names(ridge_coef) <- colnames(X_train_mat)

ridge_score_train <- as.numeric(X_train_mat %*% ridge_coef)

cor(ridge_score_train, y_train_ridge)

Final_coefficients = cbind.data.frame(selected_features, ridge_coef)
colnames(Final_coefficients) = c("Gene","Chr","Category","N_SNV","STAAR_B","BETA")

df_train_eval$RV_PGS <- ridge_score_train

mod_fit_pred = glm(Y~SL_pred, family = binomial, data = df_train_eval)
summary(mod_fit_pred)

mod_fit_approx_score = glm(Y~RV_PGS, family = binomial, data = df_train_eval)
summary(mod_fit_approx_score)

# Apply coefficients to the validation set #

X_valid_df <- X_Aligned[Phenotype_Data$Participant_ID %in% validation_set,] %>% as.data.frame()

pheno_valid = colorectal_pheno %>%
  filter(internal_set %in% c("Tuning","Test")) %>% 
  as.data.frame()

valid_features <- intersect(names(ridge_coef), colnames(X_valid_df))
feature_coefs <- ridge_coef[valid_features]

X_valid_df <- X_valid_df[, valid_features, drop = FALSE]

RV_PGS <- as.numeric(as.matrix(X_valid_df) %*% matrix(ridge_coef,ncol = 1))

pheno_valid$RV_PGS <- RV_PGS

colnames(pheno_valid) = c("Participant_ID", "Age", "Female", paste0("PC",1:20),
                          "Pop","Set","colorectal_cancer","raw_RV_PGS")

ggplot(data = pheno_valid, aes(x=raw_RV_PGS, fill = Pop)) +
  geom_density(alpha= 0.4) +
  theme_classic()


boxplot(raw_RV_PGS~colorectal_cancer, data = pheno_valid)

test_mod = glm(colorectal_cancer~raw_RV_PGS+Female+Age, family = binomial,
               data = pheno_valid)

summary(test_mod)


fwrite(Final_coefficients, "UKB_colorectal_RV_PGS_coefficients.txt", row.names = FALSE)

system("dx upload UKB_colorectal_RV_PGS_coefficients.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/")

fwrite(pheno_valid, "UKB_RV_colorectal_PGS_internal_validation_dataset.txt", row.names = FALSE)

system("dx upload UKB_RV_colorectal_PGS_internal_validation_dataset.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/")

# Upload script #
system("dx upload Colorectal_RICE_RV_training_tuning_SuperLearner_expanded_features.R --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/scripts/Rare_variant_PGS/RICE/")
