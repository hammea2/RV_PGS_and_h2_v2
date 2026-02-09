# Setup #
library(pacman)

p_load(seqinr,fastDummies,
       rlang, tidyverse, data.table,
       R.utils, devtools)

## Install or load BHR ##

#devtools::install_github("ajaynadig/bhr", force = TRUE)
library(bhr)

# Load variant annotations #
UKB_annotations = NULL
for (i in 1:22) {
  chr_variants = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/variant_sets/final_variant_sets/chr",i,"_final_variants.txt")) %>%
    unlist() %>%
    as.vector()
  chr_annotations = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/variant_sets/filtered_gene_set_variants_to_QC/filtered_chr",i,".tsv.gz")) %>%
    select(c(Uploaded_variation,Gene,SYMBOL,Consequence)) %>%
    filter(Uploaded_variation %in% chr_variants)
  UKB_annotations = rbind.data.frame(UKB_annotations,chr_annotations)
}
rm(chr_annotations)

# Load SAIGE SKAT-O results to define significant genes #
EUR_breast = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_breast_SAIGE_results.txt"))
EUR_colorectal = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_colorectal_SAIGE_results.txt"))
EUR_prostate = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_prostate_SAIGE_results.txt"))


pLoF_sig_genes_breast = EUR_breast %>%
  filter(Group  == "pLoF", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)
pLoF_sig_genes_colorectal = EUR_colorectal %>%
  filter(Group  == "pLoF", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)
pLoF_sig_genes_prostate = EUR_prostate %>%
  filter(Group  == "pLoF", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

missense_sig_genes_breast = EUR_breast %>%
  filter(Group  == "missense", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)
missense_sig_genes_colorectal = EUR_colorectal %>%
  filter(Group  == "missense", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)
missense_sig_genes_prostate = EUR_prostate %>%
  filter(Group  == "missense", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

synonymous_sig_genes_breast = EUR_breast %>%
  filter(Group  == "synonymous", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)
synonymous_sig_genes_colorectal = EUR_colorectal %>%
  filter(Group  == "synonymous", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)
synonymous_sig_genes_prostate = EUR_prostate %>%
  filter(Group  == "synonymous", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

# Load EUR MAF reports #
process_MAF_report <- function(pop, cancer){
  MAF_report = NULL
  for(i in 1:22){
    report = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/PLINK_MAF_reports/EUR_breast_colorectal_prostate/chr",i,"_",cancer,"_",pop,".frq.cc.gz")) %>%
      filter(MAF_U + MAF_A > 0)
    report = report %>%
      mutate(prevalence = NCHROBS_A / (NCHROBS_A + NCHROBS_U),
             AF = (MAF_A*prevalence) + (MAF_U*(1-prevalence)))
    MAF_report = rbind.data.frame(MAF_report,report)
  }
  return(MAF_report)
}

EUR_breast_single = process_MAF_report("EUR","breast")
EUR_colorectal_single = process_MAF_report("EUR","colorectal")
EUR_prostate_single = process_MAF_report("EUR","prostate")


# Merge with annotations #
EUR_breast_single = merge(EUR_breast_single,
                          UKB_annotations,
                          by.x = "SNP",
                          by.y = "Uploaded_variation",
                          all.x = TRUE)

EUR_breast_single$SYMBOL[EUR_breast_single$SYMBOL == "" | is.na(EUR_breast_single$SYMBOL)] <- EUR_breast_single$Gene[EUR_breast_single$SYMBOL == "" | is.na(EUR_breast_single$SYMBOL)]


EUR_colorectal_single = merge(EUR_colorectal_single,
                              UKB_annotations,
                              by.x = "SNP",
                              by.y = "Uploaded_variation",
                              all.x = TRUE)

EUR_colorectal_single$SYMBOL[EUR_colorectal_single$SYMBOL == "" | is.na(EUR_colorectal_single$SYMBOL)] <- EUR_colorectal_single$Gene[EUR_colorectal_single$SYMBOL == "" | is.na(EUR_colorectal_single$SYMBOL)]


EUR_prostate_single = merge(EUR_prostate_single,
                            UKB_annotations,
                            by.x = "SNP",
                            by.y = "Uploaded_variation",
                            all.x = TRUE)

EUR_prostate_single$SYMBOL[EUR_prostate_single$SYMBOL == "" | is.na(EUR_prostate_single$SYMBOL)] <- EUR_prostate_single$Gene[EUR_prostate_single$SYMBOL == "" | is.na(EUR_prostate_single$SYMBOL)]


# Get GRCh38 coding region lengths for each gene from Ensembl CDS FASTA #

fasta_file <- "/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/BHR_baseline_model/Homo_sapiens.GRCh38.cds.all.fa.gz"
cds_seqs <- read.fasta(file = fasta_file, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)

cds_info <- lapply(cds_seqs, function(seq_obj) {
  ## Extract full annotation string##
  annot <- attr(seq_obj, "Annot")
  
  ## Extract gene identifier using regex ##
  gene_match <- regmatches(annot, regexpr("gene:ENSG[0-9]+\\.[0-9]+", annot))
  gene_id <- sub("gene:", "", gene_match)
  
  ## Sequence length ##
  seq_length <- nchar(seq_obj[[1]]) - 1
  
  data.frame(Gene_ID = gene_id, CDS_Length = seq_length, stringsAsFactors = FALSE)
})

## Combine into a single data frame ##
cds_df <- do.call(rbind, cds_info) %>%
  mutate(gene = sub("\\.[0-9]+$", "", Gene_ID)) %>%
  group_by(gene) %>%
  summarise(Median_CDS_Length = median(CDS_Length), .groups = "drop")

## Filter to tested genes ##
cds_data = cds_df %>% 
  filter(gene %in% UKB_annotations$Gene)

# Get BHR baseline model file and expected pLoF quintile #

baseline = fread("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/BHR_baseline_model/ms_baseline_oe5.txt")
baseline = baseline %>%
  mutate(pLoF_quin = case_when(baseline_oe1_total5 == 1 ~ 1,
                               baseline_oe2_total5 == 1 ~ 2,
                               baseline_oe3_total5 == 1 ~ 3,
                               baseline_oe4_total5 == 1 ~ 4,
                               (baseline_oe1_total5 & baseline_oe2_total5 & baseline_oe3_total5 & baseline_oe4_total5) == 0 ~ 5)) %>%
  select(c(gene,pLoF_quin))

# Merge baseline BHR with CDS lengths #
baseline = merge(baseline, cds_data,
                 by.x = "gene",
                 by.y = "gene",
                 all.y = TRUE)

# Generate deciles of median CDS length to the power of pLoF quintile #
baseline = baseline %>%
  mutate(product = case_when(is.na(pLoF_quin) == FALSE ~ Median_CDS_Length^pLoF_quin,
                             is.na(pLoF_quin) == TRUE ~ Median_CDS_Length^2.5),
         decile = ntile(product,10))

baseline = baseline %>%
  select(c(gene,decile))

baseline = fastDummies::dummy_cols(.data = baseline, select_columns = "decile",
                          remove_first_dummy = TRUE) %>%
  select(-c(decile))

# Merge baseline with gene symbols

SYMBOL_and_GENE <- UKB_annotations %>%
  select(Gene, SYMBOL) %>%
  distinct()

colnames(SYMBOL_and_GENE) = c("gene","SYMBOL")

baseline <- baseline %>%
  mutate(Gene_original = gene)

baseline <- merge(SYMBOL_and_GENE, baseline, by.x = "gene", by.y = "gene", all = FALSE) %>%
  mutate(
    gene = ifelse(is.na(SYMBOL) | SYMBOL == "", Gene_original, SYMBOL)
  ) %>%
  select(-SYMBOL, -Gene_original) %>%
  as.data.frame()


# Prepare summary statistics for BHR #

## Define function to prepare summary statistics ##
prep_BHR_sumstats <- function(data, pheno, pop) {
  
  # Calculate per-sd betas
  data <- data %>%
    mutate(
      lin_beta = (2 * (MAF_A - AF) * prevalence) /
        sqrt(2 * AF * (1 - AF) * prevalence * (1 - prevalence)),
      twopq = 2 * AF * (1 - AF),
      beta_perallele = lin_beta / sqrt(twopq),
      gene_position_raw = suppressWarnings(parse_number(str_extract(SNP, "(?<=:)[0-9]+")))
    )
  
  # Summarize to gene-level
  sumstats <- data %>%
    group_by(SYMBOL) %>%
    mutate(gene_position = median(gene_position_raw, na.rm = TRUE)) %>%
    ungroup() %>%
    transmute(
      gene = SYMBOL,
      AF = AF,
      beta = beta_perallele,
      gene_position,
      chromosome = as.numeric(CHR),
      N = (NCHROBS_A + NCHROBS_U)/2,
      phenotype_key = pheno
    )
  
  return(sumstats)
}

# Define annotation categories
pLoF_variants <- c(
  "transcript_ablation", "splice_acceptor_variant",
  "splice_donor_variant", "stop_gained", 
  "frameshift_variant", "stop_lost",
  "start_lost", "transcript_elongation",
  "feature_elongation", "feature_truncation"
)

# Define annotation helper
annotate_variants <- function(df) {
  df %>%
    mutate(annotation = case_when(
      Consequence %in% pLoF_variants ~ "pLoF",
      Consequence == "missense_variant" ~ "missense",
      Consequence == "synonymous_variant" ~ "synonymous",
      TRUE ~ NA_character_
    ))
}

# Define AF bins
bin_AF <- function(df) {
  df %>%
    mutate(
      MAC = AF * (NCHROBS_A + NCHROBS_U)  # calculate MAC
    ) %>%
    mutate(
      AF_bin = case_when(
        AF < 0.00005 ~ 1,
        AF >= 0.00005 & AF < 0.0001 ~ 2,
        AF >= 0.0001 ~ 3,
        TRUE    ~ NA_real_
      )
    )
  
}

# Main processing function that extracts metadata from dataset name
process_single_dataset <- function(data, dataset_name) {
  # Expect dataset_name in format: POP_PHENO_single
  tokens <- unlist(strsplit(dataset_name, "_"))
  if (length(tokens) < 3) stop("Dataset name format should be <POP>_<PHENO>_single")
  
  pop <- tokens[1]
  pheno <- tokens[2]
  
  process_by_annotation_and_AF(data, pheno, pop)
}

# General processing function
process_by_annotation_and_AF <- function(df, pheno, pop) {
  df <- df %>%
    annotate_variants() %>%
    filter(AF <= 0.001) %>%
    bin_AF()
  
  # Split by annotation and AF_ctrl bin
  results <- df %>%
    filter(!is.na(annotation), !is.na(AF_bin)) %>%
    group_by(annotation, AF_bin) %>%
    group_split()
  
  if (length(results) == 0) {
    message("No valid groups for ", pheno, " / ", pop)
    return(NULL)
  }
  
  # Apply prep_BHR_sumstats to each subset
  lapply(results, function(subset_df) {
    annotation <- unique(subset_df$annotation)
    af_bin <- unique(subset_df$AF_bin)
    name <- paste0("res_", pop, "_", pheno, "_", annotation, "_AF", af_bin)
    message("Processing: ", name)
    assign(name, prep_BHR_sumstats(subset_df, pheno = pheno, pop = pop), envir = .GlobalEnv)
  })
}

# Run for breast, colorectal, and prostate
process_single_dataset(EUR_breast_single,"EUR_breast_single")
process_single_dataset(EUR_colorectal_single, "EUR_colorectal_single")
process_single_dataset(EUR_prostate_single,"EUR_prostate_single")

# Run BHR #

## Helper function to parse output from BHR

parse_BHR_output <- function(bhr_result, phenotype, annotation, AF_bin) {
  add_metadata <- function(df, component_name) {
    df %>%
      mutate(
        component = component_name,
        phenotype = phenotype,
        annotation = annotation,
        AF_bin = AF_bin,
        .before = 1
      )
  }
  
  parsed <- list()
  
  # --- Extract mixed_model components ---
  mm <- bhr_result$mixed_model
  
  if (!is.null(mm$heritabilities)) {
    h2_df <- as.data.frame(mm$heritabilities)
    h2_df$metric <- rownames(mm$heritabilities)
    parsed[["heritabilities"]] <- add_metadata(h2_df, "heritabilities") %>%
      pivot_longer(-c(component, phenotype, annotation, AF_bin, metric),
                   names_to = "annotation_feature", values_to = "value")
  }
  
  if (!is.null(mm$fractions)) {
    fr_df <- as.data.frame(mm$fractions)
    fr_df$metric <- rownames(mm$fractions)
    parsed[["fractions"]] <- add_metadata(fr_df, "fractions") %>%
      pivot_longer(-c(component, phenotype, annotation, AF_bin, metric),
                   names_to = "annotation_feature", values_to = "value")
  }
  
  if (!is.null(mm$enrichments_final)) {
    enrich_df <- as.data.frame(mm$enrichments_final)
    enrich_df$metric <- rownames(mm$enrichments_final)
    parsed[["enrichments"]] <- add_metadata(enrich_df, "enrichments") %>%
      pivot_longer(-c(component, phenotype, annotation, AF_bin, metric),
                   names_to = "annotation_feature", values_to = "value")
  }
  
  if (!is.null(mm$intercept)) {
    intercept_df <- data.frame(
      metric = c("intercept", "intercept_se"),
      value = c(mm$intercept, mm$intercept_se)
    )
    parsed[["intercept"]] <- add_metadata(intercept_df, "intercept")
  }
  
  if (!is.null(mm$genomewide_h2)) {
    parsed[["genomewide_h2"]] <- add_metadata(
      data.frame(metric = "genomewide_h2", value = mm$genomewide_h2),
      "genomewide_h2"
    )
  }
  
  # --- Extract significant_genes ---
  sg <- bhr_result$significant_genes
  if (!is.null(sg$number_significant_genes)) {
    sig_genes_df <- data.frame(
      metric = c("number_significant_genes", "fraction_burdenh2_significant", "fraction_burdenh2_significant_se"),
      value = c(sg$number_significant_genes, sg$fraction_burdenh2_significant, sg$fraction_burdenh2_significant_se)
    )
    parsed[["significant_genes"]] <- add_metadata(sig_genes_df, "significant_genes")
  }
  
  # --- Extract qc metrics ---
  qc <- bhr_result$qc
  if (!is.null(qc)) {
    qc_df <- as.data.frame(qc)
    qc_df_long <- pivot_longer(qc_df, everything(), names_to = "metric", values_to = "value")
    parsed[["qc"]] <- add_metadata(qc_df_long, "qc")
  }
  
  # --- Combine all parsed outputs ---
  bind_rows(parsed)
}

## Breast ##

### Get significant genes ###

pLoF_sig = pLoF_sig_genes_breast %>%
  unlist() %>% as.vector()

missense_sig = missense_sig_genes_breast %>%
  unlist() %>% as.vector()

synonymous_sig = synonymous_sig_genes_breast %>%
  unlist() %>% as.vector()

### BHR ###

#### MAF stratified ####
all_stratified <- ls(pattern = "^res_")
results_list <- list()
  
for (pheno in c("breast")) {
  objects <- all_stratified[grepl(paste0("^res_EUR_", pheno, "_"), all_stratified)]
    
  for (obj_name in objects) {
    if (!exists(obj_name, envir = .GlobalEnv)) next
      
    df <- get(obj_name, envir = .GlobalEnv)
      
    parts <- strsplit(obj_name, "_")[[1]]
    annotation <- parts[4]
    AF_bin <- if ("AF" %in% names(df)) {
      af_min <- min(df$AF, na.rm = TRUE)
      af_max <- max(df$AF, na.rm = TRUE)
      paste0(format(af_min, digits = 1, scientific = FALSE), " - ", 
             format(af_max, digits = 1, scientific = FALSE))
    } else {
      "all"
    }    
    
    fixed_genes <- switch(annotation,
                          pLoF = pLoF_sig,
                          missense = missense_sig,
                          synonymous = synonymous_sig)
      
    result <- tryCatch({
      bhr::BHR(
        mode = "univariate",
        trait1_sumstats = df,
        annotations = list(baseline),
        fixed_genes = fixed_genes
      )
    }, error = function(e) NULL)
      
    if (!is.null(result)) {
      parsed <- parse_BHR_output(result, phenotype = pheno, annotation = annotation, AF_bin = AF_bin)
      results_list[[obj_name]] <- parsed
        
      assign(sub("^res_", "BHR_", obj_name), result, envir = .GlobalEnv)
    }
  }
}

breast_stratified_results = dplyr::bind_rows(results_list)

### Aggregate ###

aggregate_breast_pLoF_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("breast"),
  ss_list = list(res_EUR_breast_pLoF_AF1, 
                 res_EUR_breast_pLoF_AF2, 
                 res_EUR_breast_pLoF_AF3),
  annotations = list(baseline),
  fixed_genes = pLoF_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "breast",
         Group = "pLoF")


aggregate_breast_missense_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("breast"),
  ss_list = list(res_EUR_breast_missense_AF1, 
                 res_EUR_breast_missense_AF2, 
                 res_EUR_breast_missense_AF3),
  annotations = list(baseline),
  fixed_genes = missense_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "breast",
         Group = "missense")

aggregate_breast_synonymous_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("breast"),
  ss_list = list(res_EUR_breast_synonymous_AF1, 
                 res_EUR_breast_synonymous_AF2, 
                 res_EUR_breast_synonymous_AF3),
  annotations = list(baseline),
  fixed_genes = synonymous_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "breast",
         Group = "synonymous")

aggregate_breast_all_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("breast"),
  ss_list = list(res_EUR_breast_pLoF_AF1, 
                 res_EUR_breast_pLoF_AF2, 
                 res_EUR_breast_pLoF_AF3,
                 res_EUR_breast_missense_AF1, 
                 res_EUR_breast_missense_AF2, 
                 res_EUR_breast_missense_AF3,
                 res_EUR_breast_synonymous_AF1, 
                 res_EUR_breast_synonymous_AF2,
                 res_EUR_breast_synonymous_AF3),
  annotations = list(baseline),
  fixed_genes = c(pLoF_sig, missense_sig, synonymous_sig)
) %>%
  as.data.frame() %>%
  mutate(Cancer = "breast",
         Group = "All")

breast_aggregate_results = rbind.data.frame(aggregate_breast_all_results,
                                            aggregate_breast_missense_results,
                                            aggregate_breast_pLoF_results,
                                            aggregate_breast_synonymous_results)
  
## Colorectal ##

### Get significant genes ###

pLoF_sig = pLoF_sig_genes_colorectal %>%
  unlist() %>% as.vector()

missense_sig = missense_sig_genes_colorectal %>%
  unlist() %>% as.vector()

synonymous_sig = synonymous_sig_genes_colorectal %>%
  unlist() %>% as.vector()

### BHR ###

#### MAF Stratified ####

all_stratified <- ls(pattern = "^res_")
results_list <- list()

for (pheno in c("colorectal")) {
  objects <- all_stratified[grepl(paste0("^res_EUR_", pheno, "_"), all_stratified)]
  
  for (obj_name in objects) {
    if (!exists(obj_name, envir = .GlobalEnv)) next
    
    df <- get(obj_name, envir = .GlobalEnv)
    
    parts <- strsplit(obj_name, "_")[[1]]
    annotation <- parts[4]
    AF_bin <- if ("AF" %in% names(df)) {
      af_min <- min(df$AF, na.rm = TRUE)
      af_max <- max(df$AF, na.rm = TRUE)
      paste0(format(af_min, digits = 1, scientific = FALSE), " - ", 
             format(af_max, digits = 1, scientific = FALSE))
    } else {
      "all"
    }    
    
    fixed_genes <- switch(annotation,
                          pLoF = pLoF_sig,
                          missense = missense_sig,
                          synonymous = synonymous_sig)
    
    result <- tryCatch({
      bhr::BHR(
        mode = "univariate",
        trait1_sumstats = df,
        annotations = list(baseline),
        fixed_genes = fixed_genes
      )
    }, error = function(e) NULL)
    
    if (!is.null(result)) {
      parsed <- parse_BHR_output(result, phenotype = pheno, annotation = annotation, AF_bin = AF_bin)
      results_list[[obj_name]] <- parsed
      
      assign(sub("^res_", "BHR_", obj_name), result, envir = .GlobalEnv)
    }
  }
}

colorectal_stratified_results = dplyr::bind_rows(results_list)

#### Aggregate ####

aggregate_colorectal_pLoF_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("colorectal"),
  ss_list = list(res_EUR_colorectal_pLoF_AF1, 
                 res_EUR_colorectal_pLoF_AF2, 
                 res_EUR_colorectal_pLoF_AF3),
  annotations = list(baseline),
  fixed_genes = pLoF_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "colorectal",
         Group = "pLoF")


aggregate_colorectal_missense_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("colorectal"),
  ss_list = list(res_EUR_colorectal_missense_AF1, 
                 res_EUR_colorectal_missense_AF2, 
                 res_EUR_colorectal_missense_AF3),
  annotations = list(baseline),
  fixed_genes = missense_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "colorectal",
         Group = "missense")

aggregate_colorectal_synonymous_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("colorectal"),
  ss_list = list(res_EUR_colorectal_synonymous_AF1, 
                 res_EUR_colorectal_synonymous_AF2, 
                 res_EUR_colorectal_synonymous_AF3),
  annotations = list(baseline),
  fixed_genes = synonymous_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "colorectal",
         Group = "synonymous")

aggregate_colorectal_all_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("colorectal"),
  ss_list = list(res_EUR_colorectal_pLoF_AF1, 
                 res_EUR_colorectal_pLoF_AF2, 
                 res_EUR_colorectal_pLoF_AF3,
                 res_EUR_colorectal_missense_AF1, 
                 res_EUR_colorectal_missense_AF2, 
                 res_EUR_colorectal_missense_AF3,
                 res_EUR_colorectal_synonymous_AF1, 
                 res_EUR_colorectal_synonymous_AF2, 
                 res_EUR_colorectal_synonymous_AF3),
  annotations = list(baseline),
  fixed_genes = c(pLoF_sig, missense_sig, synonymous_sig)
) %>%
  as.data.frame() %>%
  mutate(Cancer = "colorectal",
         Group = "All")

colorectal_aggregate_results = rbind.data.frame(aggregate_colorectal_all_results,
                                            aggregate_colorectal_missense_results,
                                            aggregate_colorectal_pLoF_results,
                                            aggregate_colorectal_synonymous_results)

## Prostate ##

### Get significant genes ###

pLoF_sig = pLoF_sig_genes_prostate %>%
  unlist() %>% as.vector()

missense_sig = missense_sig_genes_prostate %>%
  unlist() %>% as.vector()

synonymous_sig = synonymous_sig_genes_prostate %>%
  unlist() %>% as.vector()

### BHR ###

#### MAF Stratified ####

all_stratified <- ls(pattern = "^res_")
results_list <- list()

for (pheno in c("prostate")) {
  objects <- all_stratified[grepl(paste0("^res_EUR_", pheno, "_"), all_stratified)]
  
  for (obj_name in objects) {
    if (!exists(obj_name, envir = .GlobalEnv)) next
    
    df <- get(obj_name, envir = .GlobalEnv)
    
    parts <- strsplit(obj_name, "_")[[1]]
    annotation <- parts[4]
    AF_bin <- if ("AF" %in% names(df)) {
      af_min <- min(df$AF, na.rm = TRUE)
      af_max <- max(df$AF, na.rm = TRUE)
      paste0(format(af_min, digits = 1, scientific = FALSE), " - ", 
             format(af_max, digits = 1, scientific = FALSE))
    } else {
      "all"
    }    
    
    fixed_genes <- switch(annotation,
                          pLoF = pLoF_sig,
                          missense = missense_sig,
                          synonymous = synonymous_sig)
    
    result <- tryCatch({
      bhr::BHR(
        mode = "univariate",
        trait1_sumstats = df,
        annotations = list(baseline),
        fixed_genes = fixed_genes
      )
    }, error = function(e) NULL)
    
    if (!is.null(result)) {
      parsed <- parse_BHR_output(result, phenotype = pheno, annotation = annotation, AF_bin = AF_bin)
      results_list[[obj_name]] <- parsed
      
      assign(sub("^res_", "BHR_", obj_name), result, envir = .GlobalEnv)
    }
  }
}

prostate_stratified_results = dplyr::bind_rows(results_list)

#### Aggregate ####

aggregate_prostate_pLoF_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("prostate"),
  ss_list = list(res_EUR_prostate_pLoF_AF1, 
                 res_EUR_prostate_pLoF_AF2, 
                 res_EUR_prostate_pLoF_AF3),
  annotations = list(baseline),
  fixed_genes = pLoF_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "prostate",
         Group = "pLoF")


aggregate_prostate_missense_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("prostate"),
  ss_list = list(res_EUR_prostate_missense_AF1, 
                 res_EUR_prostate_missense_AF2, 
                 res_EUR_prostate_missense_AF3),
  annotations = list(baseline),
  fixed_genes = missense_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "prostate",
         Group = "missense")

aggregate_prostate_synonymous_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("prostate"),
  ss_list = list(res_EUR_prostate_synonymous_AF1, 
                 res_EUR_prostate_synonymous_AF2, 
                 res_EUR_prostate_synonymous_AF3),
  annotations = list(baseline),
  fixed_genes = synonymous_sig
) %>%
  as.data.frame() %>%
  mutate(Cancer = "prostate",
         Group = "synonymous")

aggregate_prostate_all_results <- bhr::BHR(
  mode = "aggregate",
  trait_list = c("prostate"),
  ss_list = list(res_EUR_prostate_pLoF_AF1, 
                 res_EUR_prostate_pLoF_AF2, 
                 res_EUR_prostate_pLoF_AF3,
                 res_EUR_prostate_missense_AF1, 
                 res_EUR_prostate_missense_AF2, 
                 res_EUR_prostate_missense_AF3,
                 res_EUR_prostate_synonymous_AF1, 
                 res_EUR_prostate_synonymous_AF2, 
                 res_EUR_prostate_synonymous_AF3),
  annotations = list(baseline),
  fixed_genes = c(pLoF_sig, missense_sig, synonymous_sig)
) %>%
  as.data.frame() %>%
  mutate(Cancer = "prostate",
         Group = "All")

prostate_aggregate_results = rbind.data.frame(aggregate_prostate_all_results,
                                                aggregate_prostate_missense_results,
                                                aggregate_prostate_pLoF_results,
                                                aggregate_prostate_synonymous_results)

# Save MAF stratified results
saveRDS(breast_stratified_results, file = "/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/EUR_breast_stratified_h2_results.rds")
saveRDS(breast_aggregate_results, file = "/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/EUR_breast_aggregate_h2_results.rds")
saveRDS(colorectal_stratified_results, file = "/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/EUR_colorectal_stratified_h2_results.rds")
saveRDS(colorectal_aggregate_results, file = "/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/EUR_colorectal_aggregate_h2_results.rds")
saveRDS(prostate_stratified_results, file = "/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/EUR_prostate_stratified_h2_results.rds")
saveRDS(prostate_aggregate_results, file = "/lindstroem/austin_working/Dissertation/UKBB/Aim_2/BHR/EUR_prostate_aggregate_h2_results.rds")

