#!/bin/bash
plink --bfile chr1 --merge-list 1000G_files.txt --make-bed --out 1000G_UKB_harmonized
