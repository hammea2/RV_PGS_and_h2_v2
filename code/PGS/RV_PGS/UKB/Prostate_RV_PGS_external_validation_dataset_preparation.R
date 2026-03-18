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

# Download AoU coefficients file #
system("dx download file-J5Q0YG0JfyQzFBJGjX7YbYBX")

# Load phenotype data #
pheno = data.table::fread("aim_3_phenotypes.txt")

prostate_pheno = pheno %>%
  dplyr::filter(prostate_cancer > -1) %>%
  dplyr::select(c(Participant_ID,Age_at_recruitment,paste0("UKB_PC",1:20),
                  Superpopulation,prostate_cancer))

rm(pheno)

# Load coefficients #

AoU_prostate_RV_PGS_coefficients = data.table::fread("./AoU_prostate_RV_PGS_coefficients.txt") %>%
  filter(is.na(BETA) == FALSE)

# Extract burden scores #

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

for(chr in 1:22){
  message("Processing chromosome ", chr)
  
  sig_results <- AoU_prostate_RV_PGS_coefficients %>%
    filter(Chr == chr)
  
  if(nrow(sig_results) == 0){
    next
  }
  
  genofile <- seqOpen(paste0("./aGDS/chr", chr, ".gds"))
  Chr_Sig_Burdens <- NULL
  
  for(i in seq_len(nrow(sig_results))){
    
    gene_name <- sig_results$Gene[i]
    category  <- sig_results$Category[i]
    
    message("[", i, "/", nrow(sig_results), 
            "] Processing gene ", gene_name,
            " with annotation mask ", category)
    
    a <- tryCatch(
      Burden_Scores(
        region = "Coding",
        chr = chr,
        gene_name = gene_name,
        category = category,
        genofile = genofile,
        IDs = prostate_pheno$Participant_ID,
        rare_maf_cutoff = 0.01,
        rv_num_cutoff = 2,
        QC_label = QC_label,
        variant_type = variant_type,
        geno_missing_imputation = geno_missing_imputation,
        Annotation_dir = Annotation_dir,
        Annotation_name_catalog = Annotation_name_catalog,
        silent = TRUE
      ),
      error = function(e) NULL
    )
    
    ## If no variants selected OR error occurred → add zero column
    if(is.null(a) || ncol(as.matrix(a)) == 0){
      a <- matrix(
        0,
        nrow = length(prostate_pheno$Participant_ID),
        ncol = 1
      )
    }
    
    Chr_Sig_Burdens <- cbind(Chr_Sig_Burdens, a)
  }
  
  Sig_Burdens <- cbind(Sig_Burdens, Chr_Sig_Burdens)
  rm(Chr_Sig_Burdens)
  seqClose(genofile)
}


# Prepare RV PGS dataset #

## Merge burden scores with phenotypes and align ##
X <- data.frame(Participant_ID = prostate_pheno$Participant_ID,Sig_Burdens)
Phenotype_Data <- inner_join(prostate_pheno, X, by = "Participant_ID") %>%
  as.data.frame()
X_Aligned <- as.matrix(Phenotype_Data[,paste0("X",1:ncol(Sig_Burdens))])

X_valid <- X_Aligned[Phenotype_Data$Participant_ID %in% prostate_pheno$Participant_ID,]

## Calculate PRS in the training dataset ##

# Step 0 — create unique predictor names for training coefficients
AoU_prostate_RV_PGS_coefficients <- AoU_prostate_RV_PGS_coefficients %>%
  mutate(Predictor = paste(Gene, Category, sep = "_")) %>%
  arrange(Chr, Gene, Category)

# Step 1 — rename X_valid columns using Gene_Category
colnames(X_valid) <- AoU_prostate_RV_PGS_coefficients$Predictor

# Step 2 — remove NA weights (predictors with NA beta)
Filtered_AoU_prostate_RV_PGS_coefficients <- AoU_prostate_RV_PGS_coefficients %>%
  filter(!is.na(BETA))

# List of valid predictor names
valid_predictors <- Filtered_AoU_prostate_RV_PGS_coefficients$Predictor

# Step 3 — filter X_valid to predictors with weights
X_valid_filt <- X_valid[, valid_predictors, drop = FALSE]

# Step 4 — reorder coefficients to match column order of X_valid_filt
Filtered_AoU_prostate_RV_PGS_coefficients <- Filtered_AoU_prostate_RV_PGS_coefficients %>%
  filter(Predictor %in% valid_predictors) %>%
  arrange(match(Predictor, colnames(X_valid_filt)))

# Step 5 — extract aligned weight vector
weights_vec <- Filtered_AoU_prostate_RV_PGS_coefficients$BETA

# Diagnostic checks
cat("Columns in X_valid_filt: ", ncol(X_valid_filt), "\n")
cat("Length of weights_vec: ", length(weights_vec), "\n")

# Compute PRS = X_valid_filt * weights
PRS_Validation <- as.matrix(X_valid_filt) %*% matrix(weights_vec, ncol = 1)

# Attach IDs
PRS_Validation <- data.frame(
  Participant_ID = prostate_pheno$Participant_ID,
  PRS = PRS_Validation[,1]
)

RV_PRS <- inner_join(prostate_pheno,PRS_Validation)
colnames(RV_PRS) = c("Participant_ID", "Age", paste0("PC",1:20),
                     "Superpopulation","prostate_cancer","raw_RV_PGS")

# Check PGS #

ggplot(data = RV_PRS, aes(x=raw_RV_PGS, fill = Superpopulation)) +
  geom_density(alpha= 0.4) +
  theme_classic()

boxplot(raw_RV_PGS~prostate_cancer, data = RV_PRS)

test_mod = glm(prostate_cancer~raw_RV_PGS+Age+Superpopulation, family = binomial,
               data = RV_PRS)
summary(test_mod)

# Upload external validation dataset #

fwrite(RV_PRS, "AoU_RV_prostate_PGS_external_validation_dataset.txt", row.names = FALSE)

system("dx upload AoU_RV_prostate_PGS_external_validation_dataset.txt --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/Rare_variant_PGS/RICE/")

# Upload script #
system("dx upload Prostate_RV_PGS_external_validation_dataset_preparation.R --path project-GZqFkYQJfyQzB1G6YBg89ZG3:/austin_working/dissertation_project/Aim_3/scripts/Rare_variant_PGS/RICE/")
