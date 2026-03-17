#!/bin/bash
for i in {1..22}
    do
        CHR=${i}  
        plink2 --bfile ukb22418_c${CHR}_b0_v2 --extract chr${CHR}_UKB_1000G_harmonized.prune.in --make-bed --out chr${CHR}_UKB_1000G_harmonized_pruned
    done


