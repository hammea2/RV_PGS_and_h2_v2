version 1.0

workflow PLINK_to_GDS {
  input {
    Array[File] bed_files
    Array[File] bim_files
    Array[File] fam_files
    Array[File] variant_lists
    Array[Int] subset_chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
  }

  scatter (chr in subset_chromosomes) {
    call plink_extract {
      input:
        fam = fam_files[chr - 1],
        bim = bim_files[chr - 1],
        bed = bed_files[chr - 1],
        variants = variant_lists[chr - 1],
        chr = chr
    }

    call plink_to_gds {
      input:
        bed = plink_extract.out_bed,
        bim = plink_extract.out_bim,
        fam = plink_extract.out_fam,
        chr = chr
    }
  }

  output {
    # Only GDS files are output, PLINK files are not included
    Array[File] gds_files = plink_to_gds.gds_out
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
      --no-pheno \
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

task plink_to_gds {
  input {
    File bed
    File bim
    File fam
    Int chr
  }

  command <<<
    set -euo pipefail
    Rscript -e '
      library(SeqArray)
      bed <- "~{bed}"
      bim <- "~{bim}"
      fam <- "~{fam}"
      gds_out <- sprintf("chr%s.gds", ~{chr})
      seqBED2GDS(bed.fn=bed, bim.fn=bim, fam.fn=fam,
                    out.gdsfn=gds_out, compress.annotation="LZMA_RA",
                    compress.geno="LZMA_RA", verbose=TRUE)
    '
  >>>

  runtime {
    docker: "bugsam/vcf2gds"
    dx_instance_type: "mem1_ssd2_v2_x2"
  }

  output {
    File gds_out = "chr~{chr}.gds"
  }
}
