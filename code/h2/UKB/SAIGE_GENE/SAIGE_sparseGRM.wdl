version 1.0

workflow SAIGE_sparseGRM {
  input {
    Array[File] fam_files
    Array[File] bim_files
    Array[File] bed_files
    File sample_file
    File high_ld_regions_hg19
    Array[String] chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                  "11", "12", "13", "14", "15", "16", "17", "18",
                                  "19", "20", "21", "22"]
  }

  scatter (chr_index in range(22)) {
    call plink_filter_and_LD_prune {
      input:
        fam = fam_files[chr_index],
        bim = bim_files[chr_index],
        bed = bed_files[chr_index],
        samples = sample_file,
        high_ld = high_ld_regions_hg19,
        chr = chromosomes[chr_index]
    }
  }

  call plink_merge {
    input:
      input_beds = plink_filter_and_LD_prune.output_bed,
      input_bims = plink_filter_and_LD_prune.output_bim,
      input_fams = plink_filter_and_LD_prune.output_fam
  }

  call SAIGE_step0 {
    input:
      merged_bed = plink_merge.merged_bed,
      merged_bim = plink_merge.merged_bim,
      merged_fam = plink_merge.merged_fam
  }

  output {
    File merge_log_file = plink_merge.merge_log
    File merged_bed_file = plink_merge.merged_bed
    File merged_bim_file = plink_merge.merged_bim
    File merged_fam_file = plink_merge.merged_fam
    File matrix_file = SAIGE_step0.matrix_file
    File matrix_sample_file = SAIGE_step0.matrix_sample_file
    File SAIGE_log_file = SAIGE_step0.log_file
  }
}

task plink_filter_and_LD_prune {
  input {
    File fam
    File bim
    File bed
    File samples
    File high_ld
    String chr
  }

  command <<<
    plink \
      --bed ~{bed} \
      --bim ~{bim} \
      --fam ~{fam} \
      --keep ~{samples} \
      --geno 0.01 \
      --maf 0.01 \
      --make-bed \
      --out c~{chr}

    plink \
      --bfile c~{chr} \
      --exclude range ~{high_ld} \
      --write-snplist \
      --out include

    plink2 \
      --bfile c~{chr} \
      --extract include.snplist \
      --indep-pairwise 500 50 0.2 \
      --write-snplist \
      --out c~{chr}

    plink \
      --bfile c~{chr} \
      --extract c~{chr}.prune.in \
      --make-bed \
      --out c~{chr}_ld_pruned
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd1_v2_x4"
    dx_timeout: "2h"
  }

  output {
    File output_bed = "c~{chr}_ld_pruned.bed"
    File output_bim = "c~{chr}_ld_pruned.bim"
    File output_fam = "c~{chr}_ld_pruned.fam"
    File plink_log = "c~{chr}_ld_pruned.log"
  }
}

task plink_merge {
  input {
    Array[File] input_beds
    Array[File] input_bims
    Array[File] input_fams
  }

  command <<<
    cat ~{write_lines(input_beds)} | sed -e 's/.bed//g' > merge_list.txt
    cat merge_list.txt

    plink \
      --merge-list merge_list.txt \
      --make-bed \
      --out merged
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd2_v2_x8"
    dx_timeout: "2h"
  }

  output {
    File merged_bed = "merged.bed"
    File merged_bim = "merged.bim"
    File merged_fam = "merged.fam"
    File merge_log = "merged.log"
  }
}

task SAIGE_step0 {
  input {
    File merged_bed
    File merged_bim
    File merged_fam
  }

  command <<<
    set -euo pipefail

    /usr/bin/time -o SAIGE_step0.runinfo.txt -v createSparseGRM.R \
        --bedFile=~{merged_bed} \
        --bimFile=~{merged_bim} \
        --famFile=~{merged_fam} \
        --outputPrefix=test \
        --nThreads=32 \
        --numRandomMarkerforSparseKin=5000 \
        --relatednessCutoff=0.05
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd2_v2_x32"
    dx_timeout: "24h"
  }

  output {
    File matrix_file = glob("test_relatednessCutoff_*_randomMarkersUsed.sparseGRM.mtx")[0]
    File matrix_sample_file = glob("test_relatednessCutoff_*_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")[0]
    File log_file = "SAIGE_step0.runinfo.txt"
  }
}
