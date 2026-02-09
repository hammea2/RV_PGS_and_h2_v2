version 1.0

workflow SAIGE_GENE_sparse_EUR_breast_colorectal_prostate {
  input {
    Array[File] bed_files
    Array[File] bim_files
    Array[File] fam_files
    Array[File] group_files
    Array[File] variant_lists
    File pheno_file
    Array[File] model_files
    Array[File] VR_files
    File sparseGRM_file
    File sparseGRM_samples_file
    Array[String] phenotypes = ["breast_EUR", "colorectal_EUR", "prostate_EUR"]
    Array[Int] chromosome_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]  
  }

  scatter (chr_index in chromosome_indices) {
    call plink_extract {
      input:
        fam = fam_files[chr_index],
        bim = bim_files[chr_index],
        bed = bed_files[chr_index],
        variants = variant_lists[chr_index],
        chr = chr_index + 1  
    }

    scatter (pheno_index in range(length(phenotypes))) {
      call SAIGE_2 {
        input:
          bed = plink_extract.out_bed,
          bim = plink_extract.out_bim,
          fam = plink_extract.out_fam,
          pheno_file = pheno_file,
          phenotype = phenotypes[pheno_index],
          sparseGRM = sparseGRM_file,
          samples = sparseGRM_samples_file,
          model = model_files[pheno_index],
          VR = VR_files[pheno_index],
          group = group_files[chr_index], 
          chr = chr_index + 1  
      }
    }
  }

  output {
    Array[Array[File]] output_txt = SAIGE_2.results_txt
    Array[Array[File]] output_txt_markers = SAIGE_2.results_txt_markers
    Array[Array[File]] output_txt_single = SAIGE_2.results_single_txt
    Array[Array[File]] log_files = SAIGE_2.log_file  
  }
}

task plink_extract {
  input {
    File fam
    File bim
    File bed
    File variants
    Int chr
  }

  command <<<
    plink2 \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --no-psam-pheno \
      --extract ~{variants} \
      --set-all-var-ids @:#:\$r:\$a \
      --new-id-max-allele-len 10000 \
      --make-bed \
      --out c~{chr}
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd2_v2_x8"
    dx_timeout: "4h"
  }

  output {
    File out_bed = "c~{chr}.bed"
    File out_bim = "c~{chr}.bim"
    File out_fam = "c~{chr}.fam"
  }
}

task SAIGE_2 {
  input {
    File fam
    File bim
    File bed
    File pheno_file
    String phenotype
    File sparseGRM
    File samples
    File model
    File VR
    File group
    Int chr
  }

  command <<<
    set -euo pipefail
    /usr/bin/time -o chr~{chr}_~{phenotype}_SAIGE_2.runinfo.txt -v step2_SPAtests.R \
    --bedFile=~{bed} \
    --bimFile=~{bim} \
    --famFile=~{fam} \
    --AlleleOrder=alt-first \
    --SAIGEOutputFile="chr~{chr}_~{phenotype}_SAIGE_results.txt" \
    --chrom=~{chr} \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=~{model} \
    --varianceRatioFile=~{VR} \
    --sparseGRMFile=~{sparseGRM} \                          
    --sparseGRMSampleIDFile=~{samples} \        
    --is_Firth_beta=TRUE \
    --groupFile=~{group} \
    --annotation_in_groupTest="pLoF,missense,synonymous" \
    --maxMAF_in_groupTest=0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --is_single_in_groupTest=TRUE \
    --is_output_moreDetails=TRUE \
    --LOCO=FALSE \
    --is_fastTest=TRUE
  >>>

  runtime {
    docker: "wzhou88/saige:1.3.6"
    dx_instance_type: "mem2_ssd2_v2_x4"
  }

  output {
    File results_txt = "chr~{chr}_~{phenotype}_SAIGE_results.txt"
    File results_txt_markers = "chr~{chr}_~{phenotype}_SAIGE_results.txt.markerList.txt"
    File results_single_txt = "chr~{chr}_~{phenotype}_SAIGE_results.txt.singleAssoc.txt"
    File log_file = "chr~{chr}_~{phenotype}_SAIGE_2.runinfo.txt"
  }
}