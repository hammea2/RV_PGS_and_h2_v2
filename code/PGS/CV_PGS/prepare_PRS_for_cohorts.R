# Setup #
library(pacman)
p_load(data.table,
       tidyverse,
       R.utils)

# Functions #
get_score_data <- function(file) {
  # Ensure file exists
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  # Read all lines from gzipped file
  lines <- readLines(gzfile(file))
  
  # Split into metadata and data sections
  meta_lines <- lines[grepl("^#", lines)]
  data_lines <- lines[!grepl("^#", lines)]
  
  # Process metadata
  # Remove leading '#' and optional space
  metadata <- data.frame(
    key = sub("^#\\s*", "", meta_lines),
    stringsAsFactors = FALSE
  )
  
  # Split metadata by " - " or "="
  if (any(grepl(" - |=|=", metadata$key))) {
    metadata <- do.call(rbind, lapply(metadata$key, function(line) {
      parts <- strsplit(line, " - |=|=", perl = TRUE)[[1]]
      parts <- trimws(parts)
      if (length(parts) == 2) {
        data.frame(Meta = parts[1], Value = parts[2], stringsAsFactors = FALSE)
      } else {
        data.frame(Meta = parts[1], Value = NA, stringsAsFactors = FALSE)
      }
    }))
  } else {
    colnames(metadata) <- "Line"
  }
  
  # Read the non-metadata lines as a data.frame
  if (length(data_lines) > 0) {
    # Read from textConnection so we can use read.table()
    data <- read.table(textConnection(data_lines), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    data <- data.frame()
  }
  
  # Return both as a list
  list(metadata = metadata, data = data)
}

read_VEP <- function(x, cols = NULL) {
  # Read once
  lines <- readLines(x)
  
  # Remove '##' comment lines and single '#' header marker
  lines <- sub("^#", "", lines[!startsWith(lines, "##")])
  
  # fread directly from in-memory text
  data.table::fread(
    text = lines, sep = "\t", na.strings = "-",
    data.table = TRUE, select = cols, showProgress = FALSE
  )
}


annotate <- function(data, chrs) {
  temp <- NULL
  for (chr in chrs) {
    chr_temp <- NULL
    data_chr <- data %>%
      filter(CHR == chr)
    
    annot = fread(paste0("/lindstroem/UKBB_AoU_WGS_annotations/Meta_analyses/UKB_AoUv8/chr",chr,"_variants_in_common.tsv.gz"),
                  nThread = (parallel::detectCores()-1))
    annot = annot %>% 
      filter(POS %in% data_chr$POS) %>%
      select("Gene","SYMBOL", "Consequence", "AF", "gnomADg_AF",
             "CADD_PHRED", "CADD_RAW", "CHR", "POS",
             "allele_pair_id", "UKB_ID", "UKB_REF", "UKB_ALT",
             "AoU_ID", "AoU_REF", "AoU_ALT")
    colnames(annot) = c(
      "Gene", "SYMBOL", "Consequence", "AF_1000G", "AF_gnomADg", "CADD_PHRED", "CADD_RAW",
      "CHR", "POS", "allele_pair_id", "UKB_ID", "UKB_A2", "UKB_A1",
      "AoU_ID", "AoU_A2", "AoU_A1"
    )
    annot = annot %>%
      select(-c("Gene","SYMBOL","Consequence","CADD_PHRED","CADD_RAW"))
    chr_temp = merge(data_chr, annot, by = c("CHR", "POS", "allele_pair_id"), all.x = TRUE)
    temp <- rbind.data.frame(temp, chr_temp)
  }
  results <- temp %>%
    mutate(
      flip_UKB = case_when(
        A1 == UKB_A1 & A2 == UKB_A2 ~ FALSE,
        A1 == UKB_A2 & A2 == UKB_A1 ~ TRUE,
        TRUE ~ NA
      ),
      flip_AoU = case_when(
        A1 == AoU_A1 & A2 == AoU_A2 ~ FALSE,
        A1 == AoU_A2 & A2 == AoU_A1 ~ TRUE,
        TRUE ~ NA
      )
    )
  
  return(results)
}

# Process breast cancer PGS #
breast_data = get_score_data("/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/scoring_file_processing/PGS000004_hmPOS_GRCh38.txt.gz")
breast_meta = breast_data$metadata %>%
  as.data.frame()
breast_scoring = breast_data$data %>%
  as.data.frame() %>%
  select(c(hm_chr,hm_pos,effect_allele,other_allele,effect_weight,allelefrequency_effect))
colnames(breast_scoring) = c("CHR","POS","A1","A2","A1_weight","A1_freq")
breast_scoring = breast_scoring %>%
  mutate(cancer = "breast")
breast_scoring = breast_scoring %>%
  mutate(MAF = ifelse(A1_freq >= 0.5, (1-A1_freq), A1_freq))
breast_MAF_drop = breast_scoring %>%
  filter(MAF < 0.01)
fwrite(breast_MAF_drop,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/breast_score_drop.txt",
       sep = " ")


# Process colorectal cancer PGS #
colorectal_data = get_score_data("/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/scoring_file_processing/PGS000055_hmPOS_GRCh38.txt.gz")
colorectal_meta = colorectal_data$metadata %>%
  as.data.frame()
colorectal_scoring = colorectal_data$data %>%
  as.data.frame() %>%
  select(c(hm_chr,hm_pos,effect_allele,other_allele,effect_weight))
colnames(colorectal_scoring) = c("CHR","POS","A1","A2","A1_weight")
colorectal_scoring = colorectal_scoring %>%
  mutate(cancer = "colorectal")

# Process prostate cancer PGS #
prostate_data = get_score_data("/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/scoring_file_processing/PGS000662_hmPOS_GRCh38.txt.gz")
prostate_meta = prostate_data$metadata %>%
  as.data.frame()
prostate_scoring = prostate_data$data %>%
  as.data.frame() %>%
  select(c(hm_chr,hm_pos,effect_allele,other_allele,effect_weight,
           allelefrequency_effect_European, allelefrequency_effect_African,
           allelefrequency_effect_Asian,allelefrequency_effect_Hispanic))
colnames(prostate_scoring) = c("CHR","POS","A1","A2","A1_weight",
                               "AF_European", "AF_African", "AF_Asian",
                               "AF_Hispanic")
prostate_scoring = prostate_scoring %>%
  mutate(cancer = "prostate")
prostate_scoring = prostate_scoring %>%
  mutate(
    European_MAF  = ifelse(AF_European  >= 0.5, 1 - AF_European,  AF_European),
    African_MAF   = ifelse(AF_African   >= 0.5, 1 - AF_African,   AF_African),
    Asian_MAF     = ifelse(AF_Asian     >= 0.5, 1 - AF_Asian,     AF_Asian),
    Hispanic_MAF  = ifelse(AF_Hispanic  >= 0.5, 1 - AF_Hispanic,  AF_Hispanic),
    Max_MAF = pmax(European_MAF, African_MAF, Asian_MAF, Hispanic_MAF, na.rm = TRUE)
  )
prostate_MAF_drop = prostate_scoring %>%
  filter(Max_MAF < 0.01)
fwrite(prostate_MAF_drop,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/prostate_score_drop.txt",
       sep = " ")


# Merge and allele pairs #

breast_scoring = breast_scoring %>%
  select(c("CHR","POS","A1","A2","A1_weight","cancer"))

prostate_scoring = prostate_scoring %>%
  select(c("CHR","POS","A1","A2","A1_weight","cancer"))

scoring = rbind.data.frame(breast_scoring,colorectal_scoring,prostate_scoring)

scoring <- scoring %>%
  mutate(
    allele_pair_id = mapply(
      function(ref, alt) paste(sort(c(ref, alt)), collapse = "_"),
      A2, A1
    )
  )

# Annotate with UKB and AoU variant IDs
scoring_annotated = annotate(scoring, chrs = unique(scoring$CHR))

final_scoring = scoring_annotated %>%
  filter(is.na(flip_AoU) != TRUE, is.na(flip_UKB) != TRUE) %>%
  group_by(POS) %>%
  slice_min(
    order_by = (flip_UKB + flip_AoU),  
    with_ties = FALSE
  ) %>%
  ungroup()

# Prepare score files for each cancer in each cohort#

breast_UKB_score = final_scoring %>%
  filter(cancer == "breast") %>%
  select("UKB_ID","A1","A1_weight")

fwrite(breast_UKB_score,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/breast_UKB_score.txt",
       sep = " ")

breast_AoU_score = final_scoring %>%
  filter(cancer == "breast") %>%
  select("AoU_ID","A1","A1_weight")

fwrite(breast_AoU_score,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/breast_AoU_score.txt",
       sep = " ")

colorectal_UKB_score = final_scoring %>%
  filter(cancer == "colorectal") %>%
  select("UKB_ID", "A1","A1_weight")

fwrite(colorectal_UKB_score,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/colorectal_UKB_score.txt",
       sep = " ")

colorectal_AoU_score = final_scoring %>%
  filter(cancer == "colorectal") %>%
  select("AoU_ID","A1","A1_weight")

fwrite(colorectal_AoU_score,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/colorectal_AoU_score.txt",
       sep = " ")

prostate_UKB_score = final_scoring %>%
  filter(cancer == "prostate") %>%
  select("UKB_ID","A1","A1_weight")

fwrite(prostate_UKB_score,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/prostate_UKB_score.txt",
       sep = " ")

prostate_AoU_score = final_scoring %>%
  filter(cancer == "prostate") %>%
  select("AoU_ID","A1","A1_weight")

fwrite(prostate_AoU_score,"/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/prostate_AoU_score.txt",
       sep = " ")

# Prepare chromosome-stratified variant ID sets for each cohort #
for(i in c(1:22,"X")){
  chr_variants = final_scoring %>%
    filter(CHR == i)
  chr_UKB = chr_variants %>% select(UKB_ID)
  chr_AoU = chr_variants %>% select(AoU_ID)
  fwrite(chr_UKB,paste0("/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/UKB_PGS_chr",i,"_variants.txt"), col.names = FALSE)
  fwrite(chr_AoU,paste0("/lindstroem/austin_working/Dissertation/Aim_3/Common_variant_PGS/final_scoring_files/AoU_PGS_chr",i,"_variants.txt"), col.names = FALSE)
  
}
