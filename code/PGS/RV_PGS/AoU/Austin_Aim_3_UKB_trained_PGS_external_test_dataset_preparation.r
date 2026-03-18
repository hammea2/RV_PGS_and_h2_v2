# Setup #
#install.packages("pacman")
library(pacman)

p_load(rlang, RISCA, gdsfmt, SeqArray,
       SeqVarTools, dplyr,
       Matrix, glmnet,
       caret, GenomicFeatures, 
       TxDb.Hsapiens.UCSC.hg38.knownGene,
       BiocManager, tidyverse, data.table,
       R.utils)

#BiocManager::install("GENESIS")
library(GENESIS)

p_load_gh("xihaoli/STAAR",
          "xihaoli/MultiSTAAR",
          "zilinli1988/SCANG",
          "xihaoli/STAARpipeline")

p_load_gh('jwilliams10/RICE')

pheno = data.table::fread("./Aim_3/phenotypes/aim_3_phenotypes.txt")

pheno = data.table::fread("./Aim_3/phenotypes/aim_3_phenotypes.txt")

breast_pheno = pheno %>%
  dplyr::filter(breast_cancer > -1) %>%
  dplyr::select(c(person_id,age_at_sample,paste0("PC",1:30),
                  Superpopulation,breast_cancer))

colorectal_pheno = pheno %>%
  dplyr::filter(colorectal_cancer > -1) %>%
  dplyr::select(c(person_id,Female,age_at_sample,paste0("PC",1:30),
                  Superpopulation,colorectal_cancer))

prostate_pheno = pheno %>%
  dplyr::filter(prostate_cancer > -1) %>%
  dplyr::select(c(person_id,age_at_sample,paste0("PC",1:30),
                  Superpopulation,prostate_cancer))

rm(pheno)

UKB_breast_RV_PGS_coefficients = data.table::fread("./Aim_3/Rare_variant_PGS/RICE/UKB_breast_RV_PGS_coefficients.txt") %>%
filter(is.na(BETA) == FALSE)

UKB_colorectal_RV_PGS_coefficients = data.table::fread("./Aim_3/Rare_variant_PGS/RICE/UKB_colorectal_RV_PGS_coefficients.txt") %>%
filter(is.na(BETA) == FALSE)

UKB_prostate_RV_PGS_coefficients = data.table::fread("./Aim_3/Rare_variant_PGS/RICE/UKB_prostate_RV_PGS_coefficients.txt") %>%
filter(is.na(BETA) == FALSE)

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
  sig_results = UKB_breast_RV_PGS_coefficients %>%
    filter(Chr == chr)
  if(nrow(sig_results) == 0){
    next
  } 
  genofile <- seqOpen(paste0("./Aim_3/Rare_variant_PGS/aGDS/chr",chr,"_annotated.gds"))
  Chr_Sig_Burdens <- NULL
  for(i in 1:nrow(sig_results)){
    chr <- sig_results$Chr[i]
    gene_name <- sig_results$Gene[i]
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

# Merge burden scores with phenotypes and align #
X <- data.frame(person_id = breast_pheno$person_id,Sig_Burdens)
Phenotype_Data <- inner_join(breast_pheno, X, by = "person_id") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])


X_valid <- X_Aligned[Phenotype_Data$person_id %in% breast_pheno$person_id,]

# Calculate PRS in the training dataset

# Step 0 — create unique predictor names for training coefficients
UKB_breast_RV_PGS_coefficients <- UKB_breast_RV_PGS_coefficients %>%
  mutate(Predictor = paste(Gene, Category, sep = "_")) %>%
  arrange(Chr, Gene, Category)

# Step 1 — rename X_valid columns using Gene_Category
colnames(X_valid) <- UKB_breast_RV_PGS_coefficients$Predictor

# Step 2 — remove NA weights (predictors with NA beta)
Filtered_UKB_breast_RV_PGS_coefficients <- UKB_breast_RV_PGS_coefficients %>%
  filter(!is.na(BETA))

# List of valid predictor names
valid_predictors <- Filtered_UKB_breast_RV_PGS_coefficients$Predictor

# Step 3 — filter X_valid to predictors with weights
X_valid_filt <- X_valid[, valid_predictors, drop = FALSE]

# Step 4 — reorder coefficients to match column order of X_valid_filt
Filtered_UKB_breast_RV_PGS_coefficients <- Filtered_UKB_breast_RV_PGS_coefficients %>%
  filter(Predictor %in% valid_predictors) %>%
  arrange(match(Predictor, colnames(X_valid_filt)))

# Step 5 — extract aligned weight vector
weights_vec <- Filtered_UKB_breast_RV_PGS_coefficients$BETA

# Diagnostic checks
cat("Columns in X_valid_filt: ", ncol(X_valid_filt), "\n")
cat("Length of weights_vec: ", length(weights_vec), "\n")

# Compute PRS = X_valid_filt * weights
PRS_Validation <- as.matrix(X_valid_filt) %*% matrix(weights_vec, ncol = 1)

# Attach IDs
PRS_Validation <- data.frame(
  person_id = breast_pheno$person_id,
  PRS = PRS_Validation[,1]
)

# Prep final rare PGS for the testing dataset #
RV_PRS <- inner_join(breast_pheno,PRS_Validation)
colnames(RV_PRS) = c("person_id", "Age", paste0("PC",1:30),
                     "Superpopulation","breast_cancer","raw_RV_PGS")

ggplot(data = RV_PRS, aes(x=raw_RV_PGS, fill = Superpopulation)) +
  geom_density(alpha= 0.4) +
  theme_classic()

boxplot(raw_RV_PGS~breast_cancer, data = RV_PRS)

test_mod = glm(breast_cancer~raw_RV_PGS+Age, family = binomial,
               data = RV_PRS)
summary(test_mod)

fwrite(RV_PRS, "./Aim_3/Rare_variant_PGS/RICE/UKB_RV_breast_PGS_external_testing_dataset.txt", row.names = FALSE)

system("gsutil cp ./Aim_3/Rare_variant_PGS/RICE/UKB_RV_breast_PGS_external_testing_dataset.txt gs://fc-secure-f2ecfa45-d7db-4232-929f-8cf417be1eb0/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/" )

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
  sig_results = UKB_colorectal_RV_PGS_coefficients %>%
    filter(Chr == chr)
  if(nrow(sig_results) == 0){
    next
  } 
  genofile <- seqOpen(paste0("./Aim_3/Rare_variant_PGS/aGDS/chr",chr,"_annotated.gds"))
  Chr_Sig_Burdens <- NULL
  for(i in 1:nrow(sig_results)){
    chr <- sig_results$Chr[i]
    gene_name <- sig_results$Gene[i]
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

# Merge burden scores with phenotypes and align #
X <- data.frame(person_id = colorectal_pheno$person_id,Sig_Burdens)
Phenotype_Data <- inner_join(colorectal_pheno, X, by = "person_id") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])


X_valid <- X_Aligned[Phenotype_Data$person_id %in% colorectal_pheno$person_id,]

# Calculate PRS in the training dataset

# Step 0 — create unique predictor names for training coefficients
UKB_colorectal_RV_PGS_coefficients <- UKB_colorectal_RV_PGS_coefficients %>%
  mutate(Predictor = paste(Gene, Category, sep = "_")) %>%
  arrange(Chr, Gene, Category)

# Step 1 — rename X_valid columns using Gene_Category
colnames(X_valid) <- UKB_colorectal_RV_PGS_coefficients$Predictor

# Step 2 — remove NA weights (predictors with NA beta)
Filtered_UKB_colorectal_RV_PGS_coefficients <- UKB_colorectal_RV_PGS_coefficients %>%
  filter(!is.na(BETA))

# List of valid predictor names
valid_predictors <- Filtered_UKB_colorectal_RV_PGS_coefficients$Predictor

# Step 3 — filter X_valid to predictors with weights
X_valid_filt <- X_valid[, valid_predictors, drop = FALSE]

# Step 4 — reorder coefficients to match column order of X_valid_filt
Filtered_UKB_colorectal_RV_PGS_coefficients <- Filtered_UKB_colorectal_RV_PGS_coefficients %>%
  filter(Predictor %in% valid_predictors) %>%
  arrange(match(Predictor, colnames(X_valid_filt)))

# Step 5 — extract aligned weight vector
weights_vec <- Filtered_UKB_colorectal_RV_PGS_coefficients$BETA

# Diagnostic checks
cat("Columns in X_valid_filt: ", ncol(X_valid_filt), "\n")
cat("Length of weights_vec: ", length(weights_vec), "\n")

# Compute PRS = X_valid_filt * weights
PRS_Validation <- as.matrix(X_valid_filt) %*% matrix(weights_vec, ncol = 1)

# Attach IDs
PRS_Validation <- data.frame(
  person_id = colorectal_pheno$person_id,
  PRS = PRS_Validation[,1]
)

# Prep final rare PGS for the testing dataset #
RV_PRS <- inner_join(colorectal_pheno,PRS_Validation)
colnames(RV_PRS) = c("person_id", "Female", "Age", paste0("PC",1:30),
                     "Superpopulation","colorectal_cancer","raw_RV_PGS")

ggplot(data = RV_PRS, aes(x=raw_RV_PGS, fill = Superpopulation)) +
  geom_density(alpha= 0.4) +
  theme_classic()

boxplot(raw_RV_PGS~colorectal_cancer, data = RV_PRS)

test_mod = glm(colorectal_cancer~raw_RV_PGS+Age+Female, family = binomial,
               data = RV_PRS)
summary(test_mod)

fwrite(RV_PRS, "./Aim_3/Rare_variant_PGS/RICE/UKB_RV_colorectal_PGS_external_testing_dataset.txt", row.names = FALSE)

system("gsutil cp ./Aim_3/Rare_variant_PGS/RICE/UKB_RV_colorectal_PGS_external_testing_dataset.txt gs://fc-secure-f2ecfa45-d7db-4232-929f-8cf417be1eb0/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/" )

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
  sig_results = UKB_prostate_RV_PGS_coefficients %>%
    filter(Chr == chr)
  if(nrow(sig_results) == 0){
    next
  } 
  genofile <- seqOpen(paste0("./Aim_3/Rare_variant_PGS/aGDS/chr",chr,"_annotated.gds"))
  Chr_Sig_Burdens <- NULL
  for(i in 1:nrow(sig_results)){
    chr <- sig_results$Chr[i]
    gene_name <- sig_results$Gene[i]
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

# Merge burden scores with phenotypes and align #
X <- data.frame(person_id = prostate_pheno$person_id,Sig_Burdens)
Phenotype_Data <- inner_join(prostate_pheno, X, by = "person_id") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])


X_valid <- X_Aligned[Phenotype_Data$person_id %in% prostate_pheno$person_id,]

# Calculate PRS in the training dataset

# Step 0 — create unique predictor names for training coefficients
UKB_prostate_RV_PGS_coefficients <- UKB_prostate_RV_PGS_coefficients %>%
  mutate(Predictor = paste(Gene, Category, sep = "_")) %>%
  arrange(Chr, Gene, Category)

# Step 1 — rename X_valid columns using Gene_Category
colnames(X_valid) <- UKB_prostate_RV_PGS_coefficients$Predictor

# Step 2 — remove NA weights (predictors with NA beta)
Filtered_UKB_prostate_RV_PGS_coefficients <- UKB_prostate_RV_PGS_coefficients %>%
  filter(!is.na(BETA))

# List of valid predictor names
valid_predictors <- Filtered_UKB_prostate_RV_PGS_coefficients$Predictor

# Step 3 — filter X_valid to predictors with weights
X_valid_filt <- X_valid[, valid_predictors, drop = FALSE]

# Step 4 — reorder coefficients to match column order of X_valid_filt
Filtered_UKB_prostate_RV_PGS_coefficients <- Filtered_UKB_prostate_RV_PGS_coefficients %>%
  filter(Predictor %in% valid_predictors) %>%
  arrange(match(Predictor, colnames(X_valid_filt)))

# Step 5 — extract aligned weight vector
weights_vec <- Filtered_UKB_prostate_RV_PGS_coefficients$BETA

# Diagnostic checks
cat("Columns in X_valid_filt: ", ncol(X_valid_filt), "\n")
cat("Length of weights_vec: ", length(weights_vec), "\n")

# Compute PRS = X_valid_filt * weights
PRS_Validation <- as.matrix(X_valid_filt) %*% matrix(weights_vec, ncol = 1)

# Attach IDs
PRS_Validation <- data.frame(
  person_id = prostate_pheno$person_id,
  PRS = PRS_Validation[,1]
)

# Prep final rare PGS for the testing dataset #
RV_PRS <- inner_join(prostate_pheno,PRS_Validation)
colnames(RV_PRS) = c("person_id", "Age", paste0("PC",1:30),
                     "Superpopulation","prostate_cancer","raw_RV_PGS")

ggplot(data = RV_PRS, aes(x=raw_RV_PGS, fill = Superpopulation)) +
  geom_density(alpha= 0.4) +
  theme_classic()

boxplot(raw_RV_PGS~prostate_cancer, data = RV_PRS)

test_mod = glm(prostate_cancer ~ raw_RV_PGS + Age, family = binomial, data = RV_PRS)
summary(test_mod)

fwrite(RV_PRS, "./Aim_3/Rare_variant_PGS/RICE/UKB_RV_prostate_PGS_external_testing_dataset.txt", row.names = FALSE)

system("gsutil cp ./Aim_3/Rare_variant_PGS/RICE/UKB_RV_prostate_PGS_external_testing_dataset.txt gs://fc-secure-f2ecfa45-d7db-4232-929f-8cf417be1eb0/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/" )


