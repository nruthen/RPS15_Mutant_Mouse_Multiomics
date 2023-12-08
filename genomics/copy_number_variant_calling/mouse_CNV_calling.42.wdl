version 1.0

workflow mouse_CNV_calling {
	input {
        File target_regions
        File tumor_bam
        File normal_bam 
        Array[File] normal_bams
        File ref_fasta
        File ref_dict
        File ref_index
        File ref_annotation
        String sample_name
        String gatk_docker = "broadinstitute/gatk:4.1.8.1"
        String cnvkit_docker = "etal/cnvkit:0.9.8"
        Boolean use_gatk
        Boolean use_cnvkit
	}
    if (use_gatk) {
      call CollectReadCounts as NormalCounts {
          input:
              docker_path=gatk_docker,
              target_regions=target_regions,
              tumor_bam=normal_bam
      }
      call CreateReadCountPanelOfNormals {
    	  input:
        	  docker_path=gatk_docker,
              normal_hdf5=NormalCounts.tumor_hdf5
      }
      call CollectReadCounts as TumorCounts {
          input:
              docker_path=gatk_docker,
              target_regions=target_regions,
              tumor_bam=tumor_bam
      }
      call DenoiseReadCounts {
          input:
              docker_path=gatk_docker,
              pon_hdf5=CreateReadCountPanelOfNormals.pon_hdf5,
              sample_name=sample_name,
              tumor_hdf5=TumorCounts.tumor_hdf5
      }
      call PlotDenoisedCopyRatios {
          input:
              docker_path=gatk_docker,
              sample_name=sample_name,
              standardized_cr=DenoiseReadCounts.standardized_cr,
              denoised_cr=DenoiseReadCounts.denoised_cr

      }
      call ModelAndCallSegments {
          input:
              docker_path=gatk_docker,
              tumor_denoised_cr=DenoiseReadCounts.denoised_cr
      }
	}
    if (use_cnvkit) {
      call SamtoolsIndex as IndexTumor {
          input:
          	  docker_path=gatk_docker,
              bam_file=tumor_bam
      }
      scatter (normal_bam in normal_bams) {
      	  call SamtoolsIndex as IndexNormal {
              input:
          	  	  docker_path=gatk_docker,
                  bam_file=normal_bam
          }
      }
      call CNVKitWES {
          input:
              docker_path=cnvkit_docker,
              ref_fasta=ref_fasta,
              ref_dict=ref_dict,
              ref_annotation=ref_annotation,
              target_baits=target_regions,
              normal_bams = IndexNormal.output_bam,
              normal_bais = IndexNormal.index_file,
              tumor_bam = IndexTumor.output_bam,
              tumor_bai = IndexTumor.index_file
	  }
  	}
    output {
    	Array[File]? plots_and_metrics = PlotDenoisedCopyRatios.plots_and_metrics
        File? standardized_copy_ratios = DenoiseReadCounts.standardized_cr
        File? denoised_copy_ratios = DenoiseReadCounts.denoised_cr
        File? copy_ratio_segments = ModelAndCallSegments.copy_ratio_segments
        File? called_segments = ModelAndCallSegments.called_segments
      	Array[File]? all_other_files = ModelAndCallSegments.all_other_files
        Array[File]? cnvkit_wes_results = CNVKitWES.results
    }
}

task CollectReadCounts {
    input {
        File target_regions
        File tumor_bam
        String tumor_name = basename(tumor_bam, ".bam")
        String docker_path
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
    }
    command <<<
		samtools index ~{tumor_bam} &&\
		gatk CollectReadCounts \
			-I ~{tumor_bam} \
			-L ~{target_regions} \
			--interval-merging-rule OVERLAPPING_ONLY \
			-O "/cromwell_root/~{tumor_name}.counts.hdf5"
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        File tumor_hdf5 = "/cromwell_root/~{tumor_name}.counts.hdf5"
    }
}

task DenoiseReadCounts {
    input {
        File pon_hdf5
        File tumor_hdf5
        String docker_path
        String sample_name
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
    }
    command <<<
        gatk --java-options "-Xmx12g" DenoiseReadCounts \
            -I ~{tumor_hdf5} \
            --count-panel-of-normals ~{pon_hdf5} \
            --standardized-copy-ratios "/cromwell_root/~{sample_name}_clean.standardizedCR.tsv" \
            --denoised-copy-ratios "/cromwell_root/~{sample_name}T_clean.denoisedCR.tsv"
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        File standardized_cr = "/cromwell_root/~{sample_name}_clean.standardizedCR.tsv"
        File denoised_cr = "/cromwell_root/~{sample_name}T_clean.denoisedCR.tsv"
    }
}

task PlotDenoisedCopyRatios {
    input {
        File standardized_cr
        File denoised_cr
        File ref_dict
        File ref_subset_dict
        String docker_path
        String sample_name
        Int? contig_length = 1976
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
    }
    command <<<
		if [ -f "~{ref_subset_dict}" ]; then
			gatk PlotDenoisedCopyRatios \
				--standardized-copy-ratios ~{standardized_cr} \
				--denoised-copy-ratios ~{denoised_cr} \
				--sequence-dictionary ~{ref_subset_dict} \
				--minimum-contig-length ~{contig_length} \
				--output "/cromwell_root/plots" \
				--output-prefix ~{sample_name}
		else
			gatk PlotDenoisedCopyRatios \
				--standardized-copy-ratios ~{standardized_cr} \
				--denoised-copy-ratios ~{denoised_cr} \
				--sequence-dictionary ~{ref_dict} \
				--minimum-contig-length ~{contig_length} \
				--output "/cromwell_root/plots" \
				--output-prefix ~{sample_name}
		fi
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        Array[File] plots_and_metrics = glob("/cromwell_root/plots/*")
    }
}

task SamtoolsIndex {
    input {
        File bam_file
        String bam_name = basename(bam_file)
        String docker_path
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 8
        Int? max_retries = 3
    }
    command <<<
		samtools index ~{bam_file} &&\
		mv ~{bam_file} "/cromwell_root/~{bam_name}" &&\
		mv "$(dirname ~{bam_file})/~{bam_name}.bai" "/cromwell_root/~{bam_name}.bai"
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        File output_bam = "/cromwell_root/~{bam_name}"
        File index_file = "/cromwell_root/~{bam_name}.bai" 
    }
}

task ModelAndCallSegments {
    input {
        File tumor_denoised_cr
        String tumor_name = basename(tumor_denoised_cr, "_clean.denoisedCR.tsv")
        String docker_path
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
    }
    command <<<
		mkdir output_dir && \
		gatk ModelSegments \
			--denoised-copy-ratios ~{tumor_denoised_cr} \
			--output-prefix "~{tumor_name}" \
			-O "/cromwell_root/output_dir"
		gatk CallCopyRatioSegments \
			-I "/cromwell_root/output_dir/~{tumor_name}.cr.seg" \
			-O "/cromwell_root/output_dir/~{tumor_name}.called.seg"
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        File copy_ratio_segments = "/cromwell_root/output_dir/~{tumor_name}.cr.seg"
        File called_segments = "/cromwell_root/output_dir/~{tumor_name}.called.seg"
        Array[File] all_other_files = glob("/cromwell_root/output_dir/*")
    }
}

task CNVKitWES {
 input {
        File tumor_bam
        File tumor_bai
        Array[File] normal_bams
        Array[File] normal_bais
        File ref_fasta
        File ref_dict
        File ref_index
        File ref_annotation
        File? ref_cnn
        File target_baits 
        Int? preemptible = 3
        Int? max_retries = 2
        Int? memoryGB = 16
        Int? disk_size = 100
        Int? threads = 8
        String docker_path
        String fasta_basename = basename(ref_fasta)
        String fai_basename = basename(ref_index)
        String dict_basename = basename(ref_dict)
        String? output_dir = "/cromwell_root/results/"
        String? ref_name = "GRCm38_ensembl_102.cnn"
        Boolean need_ref
    }
    command <<<
		cp ~{ref_fasta} "/cromwell_root/~{fasta_basename}" && \
		cp ~{ref_index} "/cromwell_root/~{fai_basename}" && \
		cp ~{ref_dict} "/cromwell_root/~{dict_basename}" && \
		mkdir "/cromwell_root/normal" && \
		mkdir "/cromwell_root/tumor" && \
		cp ~{tumor_bam} "/cromwell_root/tumor" && \
		cp ~{tumor_bai} "/cromwell_root/tumor" && \
		cp ~{sep=" " normal_bais} "/cromwell_root/normal" && \
		cp ~{sep=" " normal_bams} "/cromwell_root/normal" && \
		cnvkit.py access "/cromwell_root/~{fasta_basename}" -s 10000 -o "/cromwell_root/access-10kb.grcm38.bed" && \
		cnvkit.py batch /cromwell_root/tumor/*.bam --normal /cromwell_root/normal/*.bam \
			--targets ~{target_baits} --annotate ~{ref_annotation} \
			--fasta "/cromwell_root/~{fasta_basename}" --access "/cromwell_root/access-10kb.grcm38.bed" \
			--output-reference ~{ref_name} --output-dir ~{output_dir} \
			--diagram --scatter -p ~{threads}
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        cpus: "~{threads}"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
   output {
        Array[File] results = glob("/cromwell_root/results/*")
    }
}
task CreateReadCountPanelOfNormals {
    input {
        File normal_hdf5
        Int principal_components
        String docker_path
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
    }
    command <<<
        gatk --java-options "-Xmx6500m" CreateReadCountPanelOfNormals -I ~{normal_hdf5} --minimum-interval-median-percentile 5.0 --number-of-eigensamples ~{principal_components} -O "/cromwell_root/cnvponC.pon.hdf5"
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        File pon_hdf5 = "/cromwell_root/cnvponC.pon.hdf5"
    }
}