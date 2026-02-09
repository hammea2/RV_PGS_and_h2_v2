# Setup #
rm(list = ls())
library(pacman)
p_load(rlang, tidyverse, data.table, R.utils,
       ggplot2, ggpubr)


# Consolidate chromosome-specific files into a single file per cancer #

## Prepare loop variables ##
cancers = c("breast","colorectal","prostate")
chrom = 1:22
pop = "EUR"

## Loop to consolidate results for each cancer analysis ##
for (i in cancers) {
  gene_results = NULL
  single_results = NULL
  for (j in chrom) {
      chr_gene_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/Aim_2/EUR_breast_colorectal_prostate/chr",j,"_",i,"_",pop,"_SAIGE_results.txt"))
      chr_single_results = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/Aim_2/EUR_breast_colorectal_prostate/chr",j,"_",i,"_",pop,"_SAIGE_results.txt.singleAssoc.txt"))
      chr_gene_results = chr_gene_results %>% mutate(CHR = j)
      gene_results = rbind.data.frame(gene_results,chr_gene_results)
      single_results = rbind.data.frame(single_results, chr_single_results)

    }

  fwrite(gene_results,paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/Aim_2/",pop,"_",i,"_SAIGE_results.txt"),
         row.names = FALSE, quote = FALSE, sep = "\t")
  fwrite(single_results,paste0("/lindstroem/austin_working/Dissertation/UKBB/gene_based_results/Aim_2/",pop,"_",i,"_SAIGE_results.txt.singleAssoc.txt"),
         row.names = FALSE, quote = FALSE, sep = "\t")
}

# Custom QQ plot function
gg_qq <- function(pvector, lambda, data) {
  observed <- -log10(sort(pvector))
  expected <- -log10(ppoints(length(pvector)))
  
  df <- data.frame(Expected = expected, Observed = observed)
  
  ggplot(df, aes(x = Expected, y = Observed)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = paste0(" λ = ", lambda, "; N genes tested =", length(unique(data$Region))),
         x = "Expected -log10(p)",
         y = "Observed -log10(p)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
}


# Breast #

## Load data ##
EUR_breast_gene = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_breast_SAIGE_results.txt")) 
EUR_breast_single = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_breast_SAIGE_results.txt.singleAssoc.txt")) 

## QQ plots ##
qq_plots <- EUR_breast_gene %>%
  group_by(Group) %>%
  group_split() %>%
  imap(~{
    chisq_vals <- qchisq(1 - .x$Pvalue, df = 1)  
    lambda_gc <- round(median(chisq_vals, na.rm = TRUE) / qchisq(0.5, 1), 3)
    gg_qq(pvector = .x$Pvalue, lambda = lambda_gc, data = .x)
  })

plot = ggarrange(plotlist = qq_plots,
                 ncol = 2, nrow = 2,
                 labels = "AUTO")

ggsave("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/scripts/figures/QQ/UKB_EUR_breast_gene_QQ.jpeg", plot = plot, width = 6, height = 6, dpi = 400)



## Get significant genes for each annotation type ##
pLoF_sig_genes = EUR_breast_gene %>%
  filter(Group  == "pLoF", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)
  
miss_sig_genes = EUR_breast_gene %>%
  filter(Group  == "missense", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

syn_sig_genes = EUR_breast_gene %>%
  filter(Group  == "synonymous", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)


fwrite(pLoF_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_breast_pLoF_sig_genes.txt",
       row.names = FALSE)

fwrite(miss_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_breast_missense_sig_genes.txt",
       row.names = FALSE)

fwrite(syn_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_breast_synonymous_sig_genes.txt",
       row.names = FALSE)

## Report P-values for singificant genes ##

sig_tresh = (0.05/(length(unique(EUR_breast_gene$Region))))
EUR_gene_sig = EUR_breast_gene %>% filter(Pvalue < sig_tresh, Group != "Cauchy")
EUR_gene_sig = EUR_gene_sig %>% select(c(Region, Group, Pvalue))
colnames(EUR_gene_sig) = c("Gene","Group","EUR_P")

fwrite(EUR_gene_sig,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_breast_sig.csv")


# Colorectal #

## Load data ##
EUR_colorectal_gene = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_colorectal_SAIGE_results.txt")) 
EUR_colorectal_single = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_colorectal_SAIGE_results.txt.singleAssoc.txt")) 

## QQ plots ##
qq_plots <- EUR_colorectal_gene %>%
  group_by(Group) %>%
  group_split() %>%
  imap(~{
    chisq_vals <- qchisq(1 - .x$Pvalue, df = 1)  
    lambda_gc <- round(median(chisq_vals, na.rm = TRUE) / qchisq(0.5, 1), 3)
    gg_qq(pvector = .x$Pvalue, lambda = lambda_gc, data = .x)
  })

plot = ggarrange(plotlist = qq_plots,
                 ncol = 2, nrow = 2,
                 labels = "AUTO")

ggsave("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/scripts/figures/QQ/UKB_EUR_colorectal_gene_QQ.jpeg", plot = plot, width = 6, height = 6, dpi = 400)



## Get significant genes for each annotation type ##
pLoF_sig_genes = EUR_colorectal_gene %>%
  filter(Group  == "pLoF", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

miss_sig_genes = EUR_colorectal_gene %>%
  filter(Group  == "missense", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

syn_sig_genes = EUR_colorectal_gene %>%
  filter(Group  == "synonymous", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)


fwrite(pLoF_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_colorectal_pLoF_sig_genes.txt",
       row.names = FALSE)

fwrite(miss_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_colorectal_missense_sig_genes.txt",
       row.names = FALSE)

fwrite(syn_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_colorectal_synonymous_sig_genes.txt",
       row.names = FALSE)

## Report P-values for singificant genes ##

sig_tresh = (0.05/(length(unique(EUR_colorectal_gene$Region))))
EUR_gene_sig = EUR_colorectal_gene %>% filter(Pvalue < sig_tresh, Group != "Cauchy")
EUR_gene_sig = EUR_gene_sig %>% select(c(Region, Group, Pvalue))
colnames(EUR_gene_sig) = c("Gene","Group","EUR_P")

fwrite(EUR_gene_sig,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_colorectal_sig.csv")


# Prostate #

## Load data ##
EUR_prostate_gene = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_prostate_SAIGE_results.txt")) 
EUR_prostate_single = fread(paste0("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/gene_based_test_results/EUR_prostate_SAIGE_results.txt.singleAssoc.txt")) 

## QQ plots ##
qq_plots <- EUR_prostate_gene %>%
  group_by(Group) %>%
  group_split() %>%
  imap(~{
    chisq_vals <- qchisq(1 - .x$Pvalue, df = 1)  
    lambda_gc <- round(median(chisq_vals, na.rm = TRUE) / qchisq(0.5, 1), 3)
    gg_qq(pvector = .x$Pvalue, lambda = lambda_gc, data = .x)
  })

plot = ggarrange(plotlist = qq_plots,
                 ncol = 2, nrow = 2,
                 labels = "AUTO")

ggsave("/lindstroem/austin_working/Dissertation/UKBB/Aim_2/scripts/figures/QQ/UKB_EUR_prostate_gene_QQ.jpeg", plot = plot, width = 6, height = 6, dpi = 400)



## Get significant genes for each annotation type ##
pLoF_sig_genes = EUR_prostate_gene %>%
  filter(Group  == "pLoF", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

miss_sig_genes = EUR_prostate_gene %>%
  filter(Group  == "missense", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)

syn_sig_genes = EUR_prostate_gene %>%
  filter(Group  == "synonymous", Pvalue < (0.05 / length(unique(Region)))) %>% 
  select(Region)


fwrite(pLoF_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_prostate_pLoF_sig_genes.txt",
       row.names = FALSE)

fwrite(miss_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_prostate_missense_sig_genes.txt",
       row.names = FALSE)

fwrite(syn_sig_genes,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_EUR_prostate_synonymous_sig_genes.txt",
       row.names = FALSE)

## Report P-values for singificant genes ##

sig_tresh = (0.05/(length(unique(EUR_prostate_gene$Region))))
EUR_gene_sig = EUR_prostate_gene %>% filter(Pvalue < sig_tresh, Group != "Cauchy")
EUR_gene_sig = EUR_gene_sig %>% select(c(Region, Group, Pvalue))
colnames(EUR_gene_sig) = c("Gene","Group","EUR_P")

fwrite(EUR_gene_sig,"/lindstroem/austin_working/Dissertation/UKBB/Aim_2/sig_genes/UKB_prostate_sig.csv")

