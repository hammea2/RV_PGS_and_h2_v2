#!/bin/bash

# Exit immediately if any command fails
set -e

# Usage message
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <chr>"
  exit 1
fi

chr=$1

echo "Running plink2 for chromosome ${chr}..."

plink2 \
  --bfile "ukb23158_c${chr}_b0_v1" \
  --no-psam-pheno \
  --extract "chr${chr}_final_variants.txt" \
  --set-all-var-ids @:#:\$r:\$a \
  --new-id-max-allele-len 10000 \
  --make-bed \
  --out "c${chr}"

echo "plink2 completed for chromosome ${chr}"
