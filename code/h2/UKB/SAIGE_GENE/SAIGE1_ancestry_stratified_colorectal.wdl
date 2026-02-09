version 1.0

workflow SAIGE_sparse_ancestry_stratified_colorectal {
  input {
    File sparseGRM_file
    File sparseGRM_samples_file
    File bed_file
    File bim_file
    File fam_file
    File WES_bed_file
    File WES_bim_file
    File WES_fam_file
    File pheno_file
    Array[String] phenotypes = ["colorectal_AFR", "colorectal_AMR", "colorectal_EAS", "colorectal_EUR", "colorectal_SAS"]
    Array[String] instance_types = ["mem2_ssd2_v2_x8", "mem2_ssd2_v2_x4", "mem2_ssd2_v2_x4", "mem2_ssd2_v2_x16", "mem2_ssd2_v2_x8"]  
  }

  call plink_merge {
    input:
      fam = fam_file,
      bim = bim_file,
      bed = bed_file,
      WES_bed = WES_bed_file,
      WES_bim = WES_bim_file,
      WES_fam = WES_fam_file
  }

  scatter (i in range(length(phenotypes))) {
    call SAIGE_1_sparse_GRM {
      input:
        sparseGRM = sparseGRM_file,
        samples = sparseGRM_samples_file,
        fam = plink_merge.out_fam_file,
        bim = plink_merge.out_bim_file,
        bed = plink_merge.out_bed_file,
        pheno_file = pheno_file,
        pheno = phenotypes[i],
        instance_type = instance_types[i]
    }
  }

  output {
    Array[File] model_files = flatten(SAIGE_1_sparse_GRM.model_file)  
    Array[File] variance_ratio_files = flatten(SAIGE_1_sparse_GRM.variance_ratio_file)  
    Array[File] log_files = flatten(SAIGE_1_sparse_GRM.log_file)  
  }
}


task plink_merge {
  input {
    File fam
    File bim
    File bed
    File WES_fam
    File WES_bim
    File WES_bed
  }

  command <<<
    plink \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --bmerge ~{WES_bed} ~{WES_bim} ~{WES_fam} \
      --make-bed \
      --out final
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd1_v2_x4"
    dx_timeout: "10"
  }

  output {
    File out_fam_file = "final.fam"
    File out_bim_file = "final.bim"
    File out_bed_file = "final.bed" 
  }
}


task SAIGE_1_sparse_GRM {
  input {
    File sparseGRM
    File samples
    File fam
    File bim
    File bed
    File pheno_file
    String pheno
    String instance_type
  }

  command <<<
    set -euo pipefail

    /usr/bin/time -o SAIGE_step1.runinfo.txt -v step1_fitNULLGLMM.R \
        --sparseGRMFile=~{sparseGRM} \                          
        --sparseGRMSampleIDFile=~{samples} \ 
        --bedFile=~{bed} \
        --bimFile=~{bim} \
        --famFile=~{fam} \                              
        --phenoFile=~{pheno_file} \                                   
        --phenoCol=~{pheno} \                                     
        --covarColList=Age_at_recruitment,Female,UKB_PC1,UKB_PC2,UKB_PC3,UKB_PC4,UKB_PC5 \                                
        --sampleIDColinphenoFile=Participant_ID \                             
        --traitType=binary \                                          
        --outputPrefix=~{pheno}_null_model \                                
        --useSparseGRMtoFitNULL=TRUE  \
        --isCateVarianceRatio=TRUE \
        --useSparseGRMforVarRatio=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: instance_type
  }

  output {
    Array[File] model_file = glob("*.rda")  
    Array[File] variance_ratio_file = glob("*varianceRatio.txt")  
    Array[File] log_file = glob("SAIGE_step1_~{pheno}.runinfo.txt")  
  }
}