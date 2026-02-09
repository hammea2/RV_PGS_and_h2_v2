version 1.0

workflow STAAR_pipeline_medium_large {

  input {
    File script
    Array[Int] chromosomes
    File null_model_rds
    Array[File] agds_files         
    Array[File] vep_files 
    Int min_rows
    Int max_rows
  }

  scatter (chr in chromosomes) {

    call run_staar_pipeline {
      input:
        script = script,
        chr            = chr,
        null_model_rds = null_model_rds,
        agds_file      = agds_files[chr - 1],  
        vep_file       = vep_files[chr - 1],    
        min_rows       = min_rows,
        max_rows       = max_rows
      }
  }

    output {
        Array[File] result_rds = run_staar_pipeline.full_results
    }
}


task run_staar_pipeline {
    input {
        File script
        Int chr
        File null_model_rds
        File agds_file
        File vep_file
        Int min_rows
        Int max_rows
    }

    command <<<
        Rscript ~{script} \
            --chr ~{chr} \
            --null_model ~{null_model_rds} \
            --agds_file ~{agds_file} \
            --vep_file ~{vep_file} \
            --min_rows ~{min_rows} \
            --max_rows ~{max_rows}
    >>>

    runtime {
        docker: "zilinli/staarpipeline:0.9.9"
        dx_instance_type: "mem2_ssd2_v2_x4"
        dx_timeout: "10h"
  }

  output {
        File full_results   = "chr~{chr}_min~{min_rows}_max~{max_rows}_STAAR_pipeline_results_full.rds"
  }
}