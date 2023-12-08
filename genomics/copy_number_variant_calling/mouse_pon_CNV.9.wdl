version 1.0

workflow mouse_pon_CNV {
	input {
        File ref_fasta
        File ref_index
        File ref_dict
		Array[File] normal_bams
        String sample_name
        String gatk_docker = "broadinstitute/gatk:4.1.8.1"
	}
    call PreprocessIntervals {
    	input:
			docker_path=gatk_docker,
            ref_fasta=ref_fasta,
            ref_index=ref_index,
            ref_dict=ref_dict
    }
	scatter (normal_bam in normal_bams) {
    	call CollectReadCounts as NormalCounts {
    		input:
				docker_path=gatk_docker,
        		normal_bam=normal_bam,
            	target_regions=PreprocessIntervals.preprocess_interval_list
		}
		call CreateCommand1 {
        	input:
            	file_path=NormalCounts.normal_hdf5
		}
    }
    call CreateCommand2 {
    	input:
        	arguments=CreateCommand1.argument
    }
    call CreateReadCountPanelOfNormals {
    	input:
        	docker_path=gatk_docker,
            arguments=CreateCommand2.arguments,
            normal_hdf5s=NormalCounts.normal_hdf5
    }
    output {
    	File pon_hdf5 = CreateReadCountPanelOfNormals.pon_hdf5
        File preprocessed_intervals = PreprocessIntervals.preprocess_interval_list
    }
}

task PreprocessIntervals {
    input {
        File ref_fasta
        File ref_index
        File ref_dict
        String docker_path
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
    }
    command <<<
		#with index and dict in same directory
		#Default for Whole Genome Analysis
		#Default for Whole Genome Analysis
		gatk PreprocessIntervals \
			-R ~{ref_fasta} \
			--bin-length 1000 \
			--padding 0 \
			--interval-merging-rule OVERLAPPING_ONLY \
			-O "/cromwell_root/targets_C.preprocessed.interval_list"
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        File preprocess_interval_list = "/cromwell_root/targets_C.preprocessed.interval_list"
    }
}

task CollectReadCounts {
	input {
    	#targets_C.preprocessed.interval_list
		File target_regions
		File normal_bam
        String normal_name = basename(normal_bam, ".bam")
		String docker_path
		Int disk_size
		# Optional runtime variables (>= 16 GB RAM recommended)
		Int? preemptible = 3
		Int? memoryGB = 16
		Int? max_retries = 3
	}
	command <<<
		samtools index ~{normal_bam} && \
		gatk CollectReadCounts \
			-I ~{normal_bam} \
			-L ~{target_regions} \
			--interval-merging-rule OVERLAPPING_ONLY \
			-O "/cromwell_root/~{normal_name}.counts.hdf5"
	>>>
	runtime {
		docker : "~{docker_path}"
		preemptible: "${preemptible}"
		disks: "local-disk " + disk_size + " HDD"
		memory: "${memoryGB}GB"
		maxRetries: "${max_retries}"
	}
	output {
		File normal_hdf5 = "/cromwell_root/~{normal_name}.counts.hdf5"
        File normal_bai = "~{normal_bam}.bai"
	}
}

task CreateCommand1 {
	input {
		# Array[File] normal_hdf5
		# Changed from String file_path so that the file is localized in cromwell_root first.
		File file_path
	}
	command <<<
		set -euxo pipefail
		echo "-I ~{file_path}"
	>>>
	runtime {
		docker : "ubuntu:latest"
		memory: "1GB"
		maxRetries: "1"
	}
	output {
		String argument = read_string(stdout())
	}
}

task CreateCommand2 {
    input {
        #Array[File] normal_bams
        Array[String] arguments
    }
	
	command <<<
		python <<OEF
		args = ""
		len(['~{sep="', '" arguments}'])
		for x in ['~{sep="', '" arguments}']:
			args += x
			args += " "
		args = args[:-1]
		print(args)
		OEF
	>>>

    runtime {
    	docker : "python:3.4"
    	memory: "1GB"
    	maxRetries: "1"
    }
    output {
		String arguments = read_string(stdout())
    }
}

task CreateReadCountPanelOfNormals {
    input {
        String arguments
        Array[File] normal_hdf5s
        Int principal_components
        String docker_path
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
    }
    command <<<
        gatk --java-options "-Xmx6500m" CreateReadCountPanelOfNormals ~{arguments} --minimum-interval-median-percentile 5.0 --number-of-eigensamples ~{principal_components} -O "/cromwell_root/cnvponC.pon.hdf5"
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