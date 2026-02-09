suppressMessages({
  library(gdsfmt)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(STAAR)
  library(MultiSTAAR)
  library(SCANG)
  library(STAARpipeline)
})

############################################################
# Parse command-line arguments using base R
############################################################
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  pos <- which(args == flag)
  if (length(pos) == 0) stop(paste("Missing argument:", flag))
  return(args[pos + 1])
}

chr              <- as.integer(get_arg("--chr"))
null_model_file  <- get_arg("--null_model")
agds_path        <- get_arg("--agds_file")
vep_path         <- get_arg("--vep_file")
min_rows         <- as.integer(get_arg("--min_rows"))
max_rows         <- as.integer(get_arg("--max_rows"))

message("Running STAAR pipeline for chr ", chr)
message("Row thresholds: min = ", min_rows, ", max = ", max_rows)

############################################################
# Load null model + fix ID structure
############################################################
message("Loading null model...")
Null_Model <- readRDS(null_model_file)

fix_id_int <- function(x) as.integer(sub("_.*$", "", x))

# Sigma_i
dn <- Null_Model$Sigma_i@Dimnames
if (!is.null(dn[[1]])) Null_Model$Sigma_i@Dimnames[[1]] <- fix_id_int(dn[[1]])
if (!is.null(dn[[2]])) Null_Model$Sigma_i@Dimnames[[2]] <- fix_id_int(dn[[2]])

# Sigma_iX
dn <- Null_Model$Sigma_iX@Dimnames
if (!is.null(dn[[1]])) Null_Model$Sigma_iX@Dimnames[[1]] <- fix_id_int(dn[[1]])

# XSigma_i
dn <- Null_Model$XSigma_i@Dimnames
if (!is.null(dn[[2]])) Null_Model$XSigma_i@Dimnames[[2]] <- fix_id_int(dn[[2]])

# id_include
Null_Model$id_include <- fix_id_int(Null_Model$id_include)

############################################################
# Load VEP annotations
############################################################
message("Loading VEP file: ", vep_path)
VEP_annotations <- read.table(
  gzfile(vep_path),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE,       
  comment.char = ""  
)
# Ensure SYMBOL column exists
if (!"SYMBOL" %in% colnames(VEP_annotations)) VEP_annotations$SYMBOL <- VEP_annotations$Gene
VEP_annotations$SYMBOL[is.na(VEP_annotations$SYMBOL) | VEP_annotations$SYMBOL == ""] <- VEP_annotations$Gene[is.na(VEP_annotations$SYMBOL) | VEP_annotations$SYMBOL == ""]

# Count rows per gene
gene_sizes <- aggregate(VEP_annotations$SYMBOL, by = list(SYMBOL = VEP_annotations$SYMBOL), FUN = length)
colnames(gene_sizes)[2] <- "N"

############################################################
# Filter genes based on thresholds
############################################################
gene_keep <- gene_sizes[gene_sizes$N >= min_rows & gene_sizes$N <= max_rows, ]
unique_genes <- unique(gene_keep$SYMBOL)

message("Genes passing thresholds: ", length(unique_genes))
rm(VEP_annotations); gc()

############################################################
# Open GDS file
############################################################
message("Opening AGDS: ", agds_path)
geno <- seqOpen(agds_path, readonly = FALSE)

# Ensure sample IDs stored as character
IDs <- seqGetData(geno, "sample.id")
seqAddValue(geno, "sample.id", as.character(IDs), replace = TRUE, compress = "LZMA_RA")

############################################################
# STAAR annotation directories
############################################################
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category",
          "MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF",
          "aPC.EpigeneticActive","aPC.EpigeneticRepressed",
          "aPC.EpigeneticTranscription","aPC.Conservation",
          "aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info",
         "/genecode_comprehensive_exonic_category","/metasvm_pred",
         "/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed",
         "/apc_epigenetics_transcription","/apc_conservation",
         "/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")

Annotation_name_catalog <- data.frame(name = name, dir = dir)
QC_label <- "annotation/filter"
variant_type <- "variant"
geno_missing_imputation <- "mean"
Annotation_dir <- "annotation/info/FunctionalAnnotation"
Use_annotation_weights <- TRUE
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive",
                     "aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity",
                     "aPC.Mappability","aPC.TF","aPC.Protein")

############################################################
# Run STAAR pipeline
############################################################
results_coding <- list()
count <- 1
message("Beginning gene-by-gene STAAR analysis...")

for (gene in unique_genes) {
  message(sprintf("%s (%d / %d)", gene, count, length(unique_genes)))
  count <- count + 1

  res <- Gene_Centric_Coding(
    chr = chr,
    gene_name = gene,
    genofile = geno,
    obj_nullmodel = Null_Model,
    rare_maf_cutoff = 0.01,
    rv_num_cutoff = 2,
    QC_label = QC_label,
    variant_type = variant_type,
    geno_missing_imputation = geno_missing_imputation,
    Annotation_dir = Annotation_dir,
    Annotation_name_catalog = Annotation_name_catalog,
    Use_annotation_weights = Use_annotation_weights,
    Annotation_name = Annotation_name,
    silent = TRUE
  )

  results_coding <- append(results_coding, list(res))
}

seqClose(geno)

############################################################
# Save full results
############################################################
full_out <- sprintf(
  "chr%d_min%d_max%d_STAAR_pipeline_results_full.rds",
  chr, min_rows, max_rows
)
saveRDS(results_coding, full_out)
message("Saved full results: ", full_out)

message("Finished chromosome ", chr)
