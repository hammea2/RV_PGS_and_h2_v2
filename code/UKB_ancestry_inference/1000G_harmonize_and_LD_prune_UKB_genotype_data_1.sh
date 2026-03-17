#!/bin/bash
for i in {1..22}
    do
        CHR=${i}
        plink2 --bfile ukb22418_c${CHR}_b0_v2 --extract chr${CHR}_UKB_common.txt --indep-pairwise 500 10 0.2 --out chr${CHR}_UKB_1000G_harmonized   
    done


