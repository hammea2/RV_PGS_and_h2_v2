#!/bin/bash

# Exit on error, unset variable, or failed pipe
set -euo pipefail

# Enable command tracing if DEBUG is set
[[ ${DEBUG:-0} == 1 ]] && set -x

# Usage check
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <chr> <cancer_type|all> <population>"
  exit 1
fi

chr=$1
cancer_arg=$2
population=$3

# Determine cancer types to run
if [ "$cancer_arg" == "all" ]; then
  cancer_types=("breast" "colorectal" "prostate")
else
  cancer_types=("$cancer_arg")
fi

# Prepare gene list: preserve order and remove duplicates (first occurrence kept)
tmp_gene_list=$(mktemp)
awk '!seen[$1]++ {print $1}' "chr${chr}_group_file.txt" > "$tmp_gene_list"
total_genes=$(wc -l < "$tmp_gene_list")

# Print summary
echo "Total unique genes: $total_genes"

# Split gene list into chunks of 10 genes each
split_prefix=$(mktemp -u)
split -l 10 "$tmp_gene_list" "${split_prefix}_chunk_"
rm "$tmp_gene_list"

# Run R script for each cancer type and each chunk
for cancer in "${cancer_types[@]}"; do
  for chunk_file in ${split_prefix}_chunk_*; do
    echo "Processing $(wc -l < "$chunk_file") genes from $chunk_file for ${cancer}..."

    Rscript run_RareEffect.R \
      --rdaFile="${cancer}_${population}_null_model.rda" \
      --chrom="${chr}" \
      --geneListFile="${chunk_file}" \
      --groupFile="chr${chr}_group_file.txt" \
      --traitType=binary \
      --bedFile="c${chr}.bed" \
      --bimFile="c${chr}.bim" \
      --famFile="c${chr}.fam" \
      --macThreshold=10 \
      --collapseLoF=FALSE \
      --outputPrefix="${cancer}_${population}_chr${chr}"

  done
done

# Clean up chunk files
rm ${split_prefix}_chunk_*
