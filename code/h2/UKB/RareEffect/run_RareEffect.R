source("RareEffect.R")

# Read source and libraries
library(SAIGE)
library(data.table)
library(Matrix)
library(optparse)

option_list <- list(
    make_option(c("--rdaFile"), type="character", default=""),
    make_option(c("--chrom"), type="character", default=""),
    make_option(c("--geneListFile"), type="character", default=""),
    make_option(c("--groupFile"), type="character", default=""),
    make_option(c("--traitType"), type="character", default=""),
    make_option(c("--bedFile"), type="character", default=""),
    make_option(c("--bimFile"), type="character", default=""),
    make_option(c("--famFile"), type="character", default=""),
    make_option(c("--macThreshold"), type="integer", default=10),
    make_option(c("--collapseLoF"), type="logical", default=FALSE),
    make_option(c("--collapsemis"), type="logical", default=FALSE),
    make_option(c("--collapsesyn"), type="logical", default=FALSE),
    make_option(c("--apply_AR"), type="logical", default=FALSE),
    make_option(c("--outputPrefix"), type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = FALSE)
print(args)

gene_list <- readLines(args$geneListFile)

combined_effect_file <- paste0(args$outputPrefix, "_all_genes_effect.txt")
combined_h2_file <- paste0(args$outputPrefix, "_all_genes_h2.txt")
error_log <- paste0(args$outputPrefix, "_error_log.txt")

for (gene in gene_list) {
    cat("Running RareEffect for gene:", gene, "\n")

    tryCatch({
        gene_output_prefix <- paste0(args$outputPrefix, "_", gene)

        run_RareEffect(
            rdaFile = args$rdaFile,
            chrom = args$chrom,
            geneName = gene,
            groupFile = args$groupFile,
            traitType = args$traitType,
            bedFile = args$bedFile,
            bimFile = args$bimFile,
            famFile = args$famFile,
            macThreshold = args$macThreshold,
            collapseLoF = args$collapseLoF,
            collapsemis = args$collapsemis,
            collapsesyn = args$collapsesyn,
            apply_AR = args$apply_AR,
            outputPrefix = gene_output_prefix
        )

        effect_out <- paste0(gene_output_prefix, "_effect.txt")
        h2_out <- paste0(gene_output_prefix, "_h2.txt")

        if (file.exists(effect_out)) {
            effect_data <- fread(effect_out)
            effect_data[, Gene := gene]
            fwrite(
                effect_data,
                combined_effect_file,
                append = file.exists(combined_effect_file),
                col.names = !file.exists(combined_effect_file)
            )
            file.remove(effect_out)
        } else {
            warning_msg <- paste("Missing effect file for gene:", gene)
            message(warning_msg)
            write(warning_msg, file = error_log, append = TRUE)
        }

        if (file.exists(h2_out)) {
            h2_data <- fread(h2_out)
            h2_data[, Gene := gene]
            fwrite(
                h2_data,
                combined_h2_file,
                append = file.exists(combined_h2_file),
                col.names = !file.exists(combined_h2_file)
            )
            file.remove(h2_out)
        } else {
            warning_msg <- paste("Missing h2 file for gene:", gene)
            message(warning_msg)
            write(warning_msg, file = error_log, append = TRUE)
        }

        rm(effect_data, h2_data, gene_output_prefix, effect_out, h2_out)
        gc(verbose = FALSE)

    }, error = function(e) {
        err_msg <- paste("Error processing gene:", gene, ":", conditionMessage(e))
        message(err_msg)
        write(err_msg, file = error_log, append = TRUE)
    })
}
