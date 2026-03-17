#!/bin/bash
for i in {1..22}
    do
        CHR=${i}  
        plink --bfile chr${CHR}_1000G --extract chr${CHR}_UKB_1000G_harmonized.prune.in --make-bed --out chr${CHR}
    done