version 1.0

workflow PLINK_MAF_breast_colorectal_prostate {
  input { 
    Array[File] bed_files       
    Array[File] bim_files 
    Array[File] fam_files 
    Array[File] pheno_files 
    Array[File] variant_lists 
    Array[String] phenotypes = [
      "breast_EUR",
      "colorectal_EUR",
      "prostate_EUR"
    ] 
    Array[Int] chr_index = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21] 
  }
    
  scatter (chr in chr_index) {
    call plink_extract_chr {
      input:
        fam      = fam_files[chr-1],
        bim      = bim_files[chr-1],
        bed      = bed_files[chr-1],
        variants = variant_lists[chr-1],
        chr      = chr
    }

    scatter (k in range(length(pheno_files))) {
      call plink_report_maf {
        input:
          chr_bed    = plink_extract_chr.out_bed,
          chr_bim    = plink_extract_chr.out_bim,
          chr_fam    = plink_extract_chr.out_fam,    
          pheno_fam  = pheno_files[k],              
          pheno      = phenotypes[k],
          chr        = chr
      }
    }
  }

  output {
    Array[Array[File]] output_txt = flatten(plink_report_maf.MAF_report)
  }
}

task plink_extract_chr {
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

task plink_report_maf {
  input {
    File chr_bed
    File chr_bim
    File chr_fam     
    File pheno_fam   
    String pheno
    Int chr
  }

  command <<< 
    plink \
      --bed ~{chr_bed} \
      --bim ~{chr_bim} \
      --fam ~{chr_fam} \
      --pheno ~{pheno_fam} \
      --mpheno 4 \
      --make-bed \
      --allow-no-sex \
      --out chr~{chr}_withpheno

    plink \
      --bfile chr~{chr}_withpheno \
      --freq case-control gz \
      --allow-no-sex \
      --out chr~{chr}_~{pheno}
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd2_v2_x4"
  }

  output {
    Array[File] MAF_report = glob("chr~{chr}_~{pheno}*.frq.cc.gz")
  }
}
