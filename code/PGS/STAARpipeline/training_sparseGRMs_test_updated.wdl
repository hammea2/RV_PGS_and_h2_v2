version 1.0

workflow FastSparseGRM_Batch {
  input {
    File plink_bed
    File plink_bim
    File plink_fam
    Array[File] sample_lists
    Array[String] sample_names
  }

  scatter(i in range(length(sample_lists))) {
    File current_sample_list   = sample_lists[i]
    String current_sample_name = sample_names[i]

    ###############################################
    # Subset PLINK
    ###############################################
    call SubsetPLINK {
      input:
        plink_bed_file       = plink_bed,
        plink_bim_file       = plink_bim,
        plink_fam_file       = plink_fam,
        sample_keep_list     = current_sample_list,
        subset_output_prefix = current_sample_name + "_subset"
    }

    ###############################################
    # KING IBD
    ###############################################
    call RunKING {
      input:
        subset_bed_file = SubsetPLINK.subset_bed_file,
        subset_bim_file = SubsetPLINK.subset_bim_file,
        subset_fam_file = SubsetPLINK.subset_fam_file,
        ibd_degree      = 4
    }

    ###############################################
    # Get Ancestry Divergence Estimates
    ###############################################
    call GetDivergence {
      input:
        grm_bed_file      = SubsetPLINK.subset_bed_file,
        grm_bim_file      = SubsetPLINK.subset_bim_file,
        grm_fam_file      = SubsetPLINK.subset_fam_file,
        king_seg_file     = RunKING.king_seg_file,
        degree            = 4,
        divThresh         = -0.02209709,
        nRandomSNPs       = 0,
        prefix_out        = current_sample_name + "_divergence"
    }

    ###############################################
    # Extract unrelated samples
    ###############################################
    call ExtractUnrelated {
      input:
        grm_bed_file      = SubsetPLINK.subset_bed_file,
        grm_bim_file      = SubsetPLINK.subset_bim_file,
        grm_fam_file      = SubsetPLINK.subset_fam_file,
        king_seg_file     = RunKING.king_seg_file,
        unrelated_degree  = 4,
        divergence_file   = GetDivergence.divergence_file,
        unrelated_prefix  = current_sample_name + "_unrelated"
    }

    ###############################################
    # Calculate PCA
    ###############################################
    call CalcPCA {
      input:
        grm_bed_file      = SubsetPLINK.subset_bed_file,
        grm_bim_file      = SubsetPLINK.subset_bim_file,
        grm_fam_file      = SubsetPLINK.subset_fam_file,
        unrelated_list_file   = ExtractUnrelated.unrelated_list_file,
        prefix_out            = current_sample_name + "_PCs",
        no_iter               = 10
    }

    ###############################################
    # FastSparseGRM
    ###############################################
    call CalcSparseGRM {
      input:
        grm_bed_file      = SubsetPLINK.subset_bed_file,
        grm_bim_file      = SubsetPLINK.subset_bim_file,
        grm_fam_file      = SubsetPLINK.subset_fam_file,
        grm_king_seg_file = RunKING.king_seg_file,
        grm_train_file    = ExtractUnrelated.unrelated_list_file,
        grm_score_file    = CalcPCA.pc_score_file,
        grm_degree        = 4,
        grm_output_prefix = current_sample_name + "_SparseGRM",
        grm_bed_prefix    = current_sample_name + "_subset"
    }
  }

  output {
    Array[File] sparse_rdata_files     = CalcSparseGRM.sparse_rdata_file
    Array[File] king_seg_files         = RunKING.king_seg_file
    Array[File] divergence_files       = GetDivergence.divergence_file
    Array[File] unrelated_list_files   = ExtractUnrelated.unrelated_list_file
    Array[File] pc_score_files         = CalcPCA.pc_score_file
  }
}

###############################################
# Task: Subset PLINK
###############################################
task SubsetPLINK {
  input {
    File plink_bed_file
    File plink_bim_file
    File plink_fam_file
    File sample_keep_list
    String subset_output_prefix
  }

  command <<<
    set -euo pipefail
    plink2 \
      --bed ~{plink_bed_file} \
      --bim ~{plink_bim_file} \
      --fam ~{plink_fam_file} \
      --keep ~{sample_keep_list} \
      --make-bed \
      --out ~{subset_output_prefix}
  >>>

  runtime {
    docker: "austinhamm/saige_v1.4.3_plus_utilities:latest"
    dx_instance_type: "mem2_ssd2_v2_x8"
    dx_timeout: "4h"
  }

  output {
    File subset_bed_file = "~{subset_output_prefix}.bed"
    File subset_bim_file = "~{subset_output_prefix}.bim"
    File subset_fam_file = "~{subset_output_prefix}.fam"
  }
}

###############################################
# Task: KING IBD
###############################################
task RunKING {
  input {
    File subset_bed_file
    File subset_bim_file
    File subset_fam_file
    Int ibd_degree
  }

  command <<<
    set -euo pipefail
    chmod +x /FastSparseGRM/extdata/king
    export PATH="/FastSparseGRM/extdata:$PATH"
    king -b ~{subset_bed_file} --ibdseg --degree ~{ibd_degree} --cpus 16 --prefix king_out
  >>>

  runtime {
    docker: "mjohnsonngi/fastsparsegrm"
    dx_instance_type: "mem2_ssd2_v2_x16"
    dx_timeout: "4h"
  }

  output {
    File king_seg_file = "king_out.seg"
  }
}

###############################################
# Task: Get Ancestry Divergence Estimates
###############################################
task GetDivergence {
  input {
    File grm_bed_file
    File grm_bim_file
    File grm_fam_file
    File king_seg_file
    Int degree
    Float divThresh
    Int nRandomSNPs
    String prefix_out
  }

  command <<<
    set -euo pipefail
    plink_prefix="$(echo "~{grm_bed_file}" | sed 's/\.bed$//')"
    Rscript --vanilla /FastSparseGRM/extdata/getDivergence_wrapper.R \
      --prefix.in "${plink_prefix}" \
      --file.seg "~{king_seg_file}" \
      --num_threads 4 \
      --degree "~{degree}" \
      --divThresh "~{divThresh}" \
      --nRandomSNPs "~{nRandomSNPs}" \
      --prefix.out "~{prefix_out}"
  >>>

  runtime {
    docker: "mjohnsonngi/fastsparsegrm"
    dx_instance_type: "mem2_ssd2_v2_x4"
    dx_timeout: "4h"
  }

  output {
    File divergence_file = "~{prefix_out}.div"
  }
}

###############################################
# Task: Extract Unrelated Samples
###############################################
task ExtractUnrelated {
  input {
    File grm_bed_file
    File grm_bim_file
    File grm_fam_file
    File king_seg_file
    Int unrelated_degree
    File divergence_file
    String unrelated_prefix
  }

  command <<<
    set -euo pipefail
    plink_prefix="$(echo "~{grm_bed_file}" | sed 's/\.bed$//')"
    Rscript --vanilla /FastSparseGRM/extdata/extractUnrelated_wrapper.R \
      --prefix.in ${plink_prefix} \
      --file.seg ~{king_seg_file} \
      --degree ~{unrelated_degree} \
      --file.div ~{divergence_file} \
      --file.include "" \
      --prefix.out ~{unrelated_prefix}
  >>>

  runtime {
    docker: "mjohnsonngi/fastsparsegrm"
    dx_instance_type: "mem2_ssd2_v2_x4"
    dx_timeout: "4h"
  }

  output {
    File unrelated_list_file = "~{unrelated_prefix}.unrels"
  }
}

###############################################
# Task: Calculate PCA
###############################################
task CalcPCA {
  input {
    File grm_bed_file
    File grm_bim_file
    File grm_fam_file
    File unrelated_list_file
    String prefix_out
    Int no_iter
  }

  command <<<
    set -euo pipefail
    plink_prefix="$(echo "~{grm_bed_file}" | sed 's/\.bed$//')"
    Rscript --vanilla /FastSparseGRM/extdata/runPCA_wrapper.R \
      --prefix.in ${plink_prefix} \
      --file.unrels ~{unrelated_list_file} \
      --no_pcs 20 \
      --num_threads 8 \
      --prefix.out ~{prefix_out} \
      --no_iter ~{no_iter}
  >>>

  runtime {
    docker: "mjohnsonngi/fastsparsegrm"
    dx_instance_type: "mem2_ssd2_v2_x8"
    dx_timeout: "4h"
  }

  output {
    File pc_score_file = "~{prefix_out}.score"
  }
}

###############################################
# Task: FastSparseGRM
###############################################
task CalcSparseGRM {
  input {
    File grm_bed_file
    File grm_bim_file
    File grm_fam_file
    File grm_king_seg_file
    File grm_train_file
    File grm_score_file
    Int grm_degree
    String grm_output_prefix
    String grm_bed_prefix
  }

  command <<<
    set -euo pipefail
    plink_prefix="$(echo "~{grm_bed_file}" | sed 's/\.bed$//')"
    Rscript --vanilla /FastSparseGRM/extdata/calcSparseGRM_wrapper.R \
      --prefix.in ${plink_prefix} \
      --prefix.out ~{grm_output_prefix} \
      --file.train ~{grm_train_file} \
      --file.score ~{grm_score_file} \
      --file.seg ~{grm_king_seg_file} \
      --num_threads 16 \
      --no_pcs 20 \
      --block.size 5000 \
      --max.related.block 5000 \
      --KINGformat.out FALSE \
      --degree ~{grm_degree}
  >>>

  runtime {
    docker: "mjohnsonngi/fastsparsegrm"
    dx_instance_type: "mem2_ssd2_v2_x16"
    dx_timeout: "24h"
  }

  output {
    File sparse_rdata_file = "~{grm_output_prefix}.sGRM.RData"
  }
}