# Setup #
#install.packages("pacman")
library(pacman)

p_load(rlang, RISCA, gdsfmt, SeqArray,
       SeqVarTools, dplyr,
       Matrix, glmnet,
       caret, GenomicFeatures, 
       TxDb.Hsapiens.UCSC.hg38.knownGene,
       BiocManager, tidyverse, data.table,
       R.utils, SuperLearner, ranger)

#BiocManager::install("GENESIS")
library(GENESIS)

p_load_gh("xihaoli/STAAR",
          "xihaoli/MultiSTAAR",
          "zilinli1988/SCANG",
          "xihaoli/STAARpipeline")

p_load_gh('jwilliams10/RICE')

pheno = data.table::fread("./Aim_3/phenotypes/aim_3_phenotypes.txt")

breast_pheno = pheno %>%
  dplyr::filter(breast_cancer > -1) %>%
  dplyr::select(c(person_id,age_at_sample,paste0("PC",1:30),
                  Superpopulation,breast_cancer,internal_set))

colorectal_pheno = pheno %>%
  dplyr::filter(colorectal_cancer > -1) %>%
  dplyr::select(c(person_id,Female,age_at_sample,paste0("PC",1:30),
                  Superpopulation,colorectal_cancer,internal_set))

prostate_pheno = pheno %>%
  dplyr::filter(prostate_cancer > -1) %>%
  dplyr::select(c(person_id,age_at_sample,paste0("PC",1:30),
                  Superpopulation,prostate_cancer,internal_set))

rm(pheno)

breast_selected_features <- readRDS("./Aim_3/Rare_variant_PGS/RICE/breast_selected_features_expanded.rds")
colorectal_selected_features <- readRDS("./Aim_3/Rare_variant_PGS/RICE/colorectal_selected_features_expanded.rds")
prostate_selected_features <- readRDS("./Aim_3/Rare_variant_PGS/RICE/prostate_selected_features_expanded.rds")

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
  sig_results = breast_selected_features %>%
    filter(Chr == chr)
  if(nrow(sig_results) == 0){
    next
  } 
  genofile <- seqOpen(paste0("./Aim_3/Rare_variant_PGS/aGDS/chr",chr,"_annotated.gds"))
  Chr_Sig_Burdens <- NULL
  for(i in 1:nrow(sig_results)){
    chr <- sig_results$Chr[i]
    gene_name <- sig_results$`Gene name`[i]
    category <- sig_results$Category[i]
    print(paste0("[",i,"/",nrow(sig_results),"] ","Processing gene ", gene_name, " with annotation mask ", category))
    a <- Burden_Scores(region = "Coding",chr = chr, gene_name = gene_name,category = category,
                       genofile = genofile,IDs = breast_pheno$person_id,
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

X <- data.frame(person_id = breast_pheno$person_id,Sig_Burdens)
Phenotype_Data <- inner_join(breast_pheno, X, by = "person_id") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])

train_set = Phenotype_Data %>%
  filter(internal_set == "Training") %>%
  pull(person_id)

validation_set = Phenotype_Data %>%
  filter(internal_set %in% c("Tuning","Test")) %>%
  pull(person_id)

X_train <- X_Aligned[Phenotype_Data$person_id %in% train_set,] %>% as.data.frame()

pheno_train = breast_pheno %>%
  filter(internal_set == "Training") %>% 
  as.data.frame()
trait = "breast_cancer"
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

sl_fit$coef
sl_weights <- sort(sl_fit$coef, decreasing = TRUE)
print(sl_weights)

y_train_ridge <- sl_fit$SL.predict
X_train_mat <- as.matrix(X_train)           

df_train_eval <- data.frame(
  Y = pheno_train[, trait],
  SL_pred = y_train_ridge
)

lambda_grid <- 10^seq(3, -3, length.out = 50)

set.seed(298285)
cv_ridge <- cv.glmnet(
  X_train_mat,
  y_train_ridge,
  alpha = 0,                
  lambda = lambda_grid,     
  nfolds = 10,              
  standardize = FALSE,
  lower.limits = 0
)

best_lambda <- cv_ridge$lambda.min 

ridge_coef <- as.numeric(coef(cv_ridge, s = best_lambda))[-1] 
ridge_intercept <- as.numeric(coef(cv_ridge, s = best_lambda))[1]

names(ridge_coef) <- colnames(X_train_mat)

sort(abs(ridge_coef), decreasing = TRUE)

ridge_score_train <- as.numeric(X_train_mat %*% ridge_coef)

cor(ridge_score_train, y_train_ridge)

Final_coefficients = cbind.data.frame(breast_selected_features, ridge_coef)
colnames(Final_coefficients) = c("Gene","Chr","Category","N_SNV","STAAR_B","BETA")

df_train_eval$RV_PGS <- ridge_score_train

mod_fit_pred = glm(Y~SL_pred, family = binomial, data = df_train_eval)
summary(mod_fit_pred)

mod_fit_approx_score = glm(Y~RV_PGS, family = binomial, data = df_train_eval)
summary(mod_fit_approx_score)

X_valid_df <- X_Aligned[Phenotype_Data$person_id %in% validation_set,] %>% as.data.frame()

pheno_valid = breast_pheno %>%
  filter(internal_set %in% c("Tuning","Test")) %>% 
  as.data.frame()

RV_PGS <- as.numeric(as.matrix(X_valid_df) %*% matrix(ridge_coef,ncol = 1))


pheno_valid$RV_PGS <- RV_PGS

colnames(pheno_valid) = c("person_id", "Age", paste0("PC",1:30),
                     "Superpopulation","breast_cancer","Set","raw_RV_PGS")

ggplot(data = pheno_valid, aes(x=raw_RV_PGS, fill = Superpopulation)) +
  geom_density(alpha= 0.4) +
  theme_classic()

boxplot(raw_RV_PGS~breast_cancer, data = pheno_valid)

test_mod = glm(breast_cancer~raw_RV_PGS+Age, family = binomial,
               data = pheno_valid)

summary(test_mod)

fwrite(Final_coefficients, "./Aim_3/Rare_variant_PGS/RICE/AoU_breast_RV_PGS_coefficients.txt", row.names = FALSE)

system("gsutil cp ./Aim_3/Rare_variant_PGS/RICE/AoU_breast_RV_PGS_coefficients.txt gs://fc-secure-f2ecfa45-d7db-4232-929f-8cf417be1eb0/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/" )

fwrite(pheno_valid, "./Aim_3/Rare_variant_PGS/RICE/AoU_RV_breast_PGS_internal_testing_dataset.txt", row.names = FALSE)

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
  sig_results = colorectal_selected_features %>%
    filter(Chr == chr)
  if(nrow(sig_results) == 0){
    next
  } 
  genofile <- seqOpen(paste0("./Aim_3/Rare_variant_PGS/aGDS/chr",chr,"_annotated.gds"))
  Chr_Sig_Burdens <- NULL
  for(i in 1:nrow(sig_results)){
    chr <- sig_results$Chr[i]
    gene_name <- sig_results$`Gene name`[i]
    category <- sig_results$Category[i]
    print(paste0("[",i,"/",nrow(sig_results),"] ","Processing gene ", gene_name, " with annotation mask ", category))
    a <- Burden_Scores(region = "Coding",chr = chr, gene_name = gene_name,category = category,
                       genofile = genofile,IDs = colorectal_pheno$person_id,
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

X <- data.frame(person_id = colorectal_pheno$person_id,Sig_Burdens)
Phenotype_Data <- inner_join(colorectal_pheno, X, by = "person_id") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])

train_set = Phenotype_Data %>%
  filter(internal_set == "Training") %>%
  pull(person_id)

validation_set = Phenotype_Data %>%
  filter(internal_set %in% c("Tuning","Test")) %>%
  pull(person_id)

X_train <- X_Aligned[Phenotype_Data$person_id %in% train_set,] %>% as.data.frame()

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

sl_fit$coef
sl_weights <- sort(sl_fit$coef, decreasing = TRUE)
print(sl_weights)

y_train_ridge <- sl_fit$SL.predict
X_train_mat <- as.matrix(X_train)           

df_train_eval <- data.frame(
  Y = pheno_train[, trait],
  SL_pred = y_train_ridge
)

lambda_grid <- 10^seq(3, -3, length.out = 50)

set.seed(298285)
cv_ridge <- cv.glmnet(
  X_train_mat,
  y_train_ridge,
  alpha = 0,                
  lambda = lambda_grid,     
  nfolds = 10,              
  standardize = FALSE,
  lower.limits = 0
)

best_lambda <- cv_ridge$lambda.min 

ridge_coef <- as.numeric(coef(cv_ridge, s = best_lambda))[-1] 
ridge_intercept <- as.numeric(coef(cv_ridge, s = best_lambda))[1]

names(ridge_coef) <- colnames(X_train_mat)

ridge_score_train <- as.numeric(X_train_mat %*% ridge_coef)

cor(ridge_score_train, y_train_ridge)

Final_coefficients = cbind.data.frame(colorectal_selected_features, ridge_coef)
colnames(Final_coefficients) = c("Gene","Chr","Category","N_SNV","STAAR_B","BETA")

df_train_eval$RV_PGS <- ridge_score_train

mod_fit_pred = glm(Y~SL_pred, family = binomial, data = df_train_eval)
summary(mod_fit_pred)

mod_fit_approx_score = glm(Y~RV_PGS, family = binomial, data = df_train_eval)
summary(mod_fit_approx_score)

X_valid_df <- X_Aligned[Phenotype_Data$person_id %in% validation_set,] %>% as.data.frame()

pheno_valid = colorectal_pheno %>%
  filter(internal_set %in% c("Tuning","Test")) %>% 
  as.data.frame()

RV_PGS <- as.numeric(as.matrix(X_valid_df) %*% matrix(ridge_coef,ncol = 1))


pheno_valid$RV_PGS <- RV_PGS

colnames(pheno_valid) = c("person_id", "Age", "Female", paste0("PC",1:30),
                     "Superpopulation","colorectal_cancer","Set","raw_RV_PGS")

ggplot(data = pheno_valid, aes(x=raw_RV_PGS, fill = Superpopulation)) +
  geom_density(alpha= 0.4) +
  theme_classic()

boxplot(raw_RV_PGS~colorectal_cancer, data = pheno_valid)

test_mod = glm(colorectal_cancer~raw_RV_PGS+Age+Female, family = binomial,
               data = pheno_valid)

summary(test_mod)

fwrite(Final_coefficients, "./Aim_3/Rare_variant_PGS/RICE/AoU_colorectal_RV_PGS_coefficients.txt", row.names = FALSE)

system("gsutil cp ./Aim_3/Rare_variant_PGS/RICE/AoU_colorectal_RV_PGS_coefficients.txt gs://fc-secure-f2ecfa45-d7db-4232-929f-8cf417be1eb0/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/" )

fwrite(pheno_valid, "./Aim_3/Rare_variant_PGS/RICE/AoU_RV_colorectal_PGS_internal_testing_dataset.txt", row.names = FALSE)

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
  sig_results = prostate_selected_features %>%
    filter(Chr == chr)
  if(nrow(sig_results) == 0){
    next
  } 
  genofile <- seqOpen(paste0("./Aim_3/Rare_variant_PGS/aGDS/chr",chr,"_annotated.gds"))
  Chr_Sig_Burdens <- NULL
  for(i in 1:nrow(sig_results)){
    chr <- sig_results$Chr[i]
    gene_name <- sig_results$`Gene name`[i]
    category <- sig_results$Category[i]
    print(paste0("[",i,"/",nrow(sig_results),"] ","Processing gene ", gene_name, " with annotation mask ", category))
    a <- Burden_Scores(region = "Coding",chr = chr, gene_name = gene_name,category = category,
                       genofile = genofile,IDs = prostate_pheno$person_id,
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

X <- data.frame(person_id = prostate_pheno$person_id,Sig_Burdens)
Phenotype_Data <- inner_join(prostate_pheno, X, by = "person_id") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])

train_set = Phenotype_Data %>%
  filter(internal_set == "Training") %>%
  pull(person_id)

validation_set = Phenotype_Data %>%
  filter(internal_set %in% c("Tuning","Test")) %>%
  pull(person_id)

X_train <- X_Aligned[Phenotype_Data$person_id %in% train_set,] %>% as.data.frame()

pheno_train = prostate_pheno %>%
  filter(internal_set == "Training") %>% 
  as.data.frame()
trait = "prostate_cancer"
Y <- as.numeric(pheno_train[, trait])

set.seed(1814879)

V <- 10 
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

sl_fit$coef
sl_weights <- sort(sl_fit$coef, decreasing = TRUE)
print(sl_weights)

y_train_ridge <- sl_fit$SL.predict
X_train_mat <- as.matrix(X_train)           

df_train_eval <- data.frame(
  Y = pheno_train[, trait],
  SL_pred = y_train_ridge
)

lambda_grid <- 10^seq(3, -3, length.out = 50)

set.seed(298285)
cv_ridge <- cv.glmnet(
  X_train_mat,
  y_train_ridge,
  alpha = 0,                
  lambda = lambda_grid,     
  nfolds = 10,              
  standardize = FALSE,
  lower.limits = 0
)

best_lambda <- cv_ridge$lambda.min 

ridge_coef <- as.numeric(coef(cv_ridge, s = best_lambda))[-1] 
ridge_intercept <- as.numeric(coef(cv_ridge, s = best_lambda))[1]

names(ridge_coef) <- colnames(X_train_mat)

ridge_score_train <- as.numeric(X_train_mat %*% ridge_coef)

cor(ridge_score_train, y_train_ridge)

Final_coefficients = cbind.data.frame(prostate_selected_features, ridge_coef)
colnames(Final_coefficients) = c("Gene","Chr","Category","N_SNV","STAAR_B","BETA")

df_train_eval$RV_PGS <- ridge_score_train

mod_fit_pred = glm(Y~SL_pred, family = binomial, data = df_train_eval)
summary(mod_fit_pred)

mod_fit_approx_score = glm(Y~RV_PGS, family = binomial, data = df_train_eval)
summary(mod_fit_approx_score)

X_valid_df <- X_Aligned[Phenotype_Data$person_id %in% validation_set,] %>% as.data.frame()

pheno_valid = prostate_pheno %>%
  filter(internal_set %in% c("Tuning","Test")) %>% 
  as.data.frame()

RV_PGS <- as.numeric(as.matrix(X_valid_df) %*% matrix(ridge_coef,ncol = 1))


pheno_valid$RV_PGS <- RV_PGS

colnames(pheno_valid) = c("person_id", "Age", paste0("PC",1:30),
                     "Superpopulation","prostate_cancer","Set","raw_RV_PGS")

ggplot(data = pheno_valid, aes(x=raw_RV_PGS, fill = Superpopulation)) +
  geom_density(alpha= 0.4) +
  theme_classic()

boxplot(raw_RV_PGS~prostate_cancer, data = pheno_valid)

test_mod = glm(prostate_cancer~raw_RV_PGS+Age, family = binomial,
               data = pheno_valid)

summary(test_mod)

fwrite(Final_coefficients, "./Aim_3/Rare_variant_PGS/RICE/AoU_prostate_RV_PGS_coefficients.txt", row.names = FALSE)

system("gsutil cp ./Aim_3/Rare_variant_PGS/RICE/AoU_prostate_RV_PGS_coefficients.txt gs://fc-secure-f2ecfa45-d7db-4232-929f-8cf417be1eb0/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/" )

fwrite(pheno_valid, "./Aim_3/Rare_variant_PGS/RICE/AoU_RV_prostate_PGS_internal_testing_dataset.txt", row.names = FALSE)
