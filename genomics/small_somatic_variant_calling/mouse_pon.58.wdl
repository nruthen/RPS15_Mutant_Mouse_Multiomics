## About:
##	This WDL workflow creates a panel of normals for the MWEGS Pipeline VariantCalling  
##	workflow. 

version 1.0 

workflow CreatePanelOfNormals {
	input {
        File ref_fasta
        File ref_index
        File ref_dict 
        Array[File] normal_bams
        #File inputNormalsFile 
        #Array[Array[File]] inputNormals = read_tsv(inputNormalsFile)
        String gatk_docker = "broadinstitute/gatk:4.1.8.1"
        String old_gatk_docker = "broadinstitute/gatk:4.0.1.2"
        String gatk_path = "/gatk/gatk"
        String samtools_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        String samtools_path = "/usr/gitc/samtools-1.3.1/"
	}
	
    scatter (normal in normal_bams) {
    	String sample_basename = basename(normal, ".bam")
    	call SamtoolsIndex {
        	input:
            	docker_path = samtools_docker,
                samtools_path = samtools_path,
                sample_name = sample_basename,
                normal_bam = normal
        }
        call Mutect2 {
            input:
                docker_path = gatk_docker,
                gatk_path = gatk_path,
                sample_name = sample_basename,
                ref_fasta = ref_fasta,
                ref_index = ref_index,
                ref_dict = ref_dict,
                normal_bam = SamtoolsIndex.output_bam,
                normal_bai = SamtoolsIndex.output_bai
        }
    }

    scatter (output_vcf in Mutect2.output_vcf) { # OR scatter(normal_bam in normal_bams) if we use an Array[File]
        call CreateCommand1 {
            input:
                file_path = output_vcf
        }
    }
    
	call CreateCommand2 {
		input:
			arguments = CreateCommand1.argument 
	}

    call CreateSomaticPanelOfNormals {
    	input: 
        	docker_path = old_gatk_docker,
        	gatk_path = gatk_path,
       		arguments = CreateCommand2.arguments,
            normal_vcfs = Mutect2.output_vcf
    }
}

task SamtoolsIndex {
	input {
		String docker_path
        String sample_name
        String samtools_path
		File normal_bam
        Int? preemptible = 1
		Int? max_retries = 1
		Int? memoryGB = 16
		#Int? diskGB_buffer = 10
		#Int diskGB = ceil(size(fastq1, 'G')*2 + size(fastq2, 'G') + diskGB_buffer)
        Int? diskGB = 100
    }
	command <<<
		set -euxo pipefail
		samtools index ~{normal_bam} ~{sample_name}.bam.bai && \
		mv ~{normal_bam} ~{sample_name}.bam
	>>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File output_bam = "~{sample_name}.bam"
		File output_bai = "~{sample_name}.bam.bai"
	}
}

task Mutect2 {
	input {
        String docker_path
        String gatk_path
        String sample_name
        File normal_bam
        File normal_bai
        File ref_fasta
        File ref_index
        File ref_dict 
		Int? preemptible = 1
		Int? max_retries = 1
		Int? memoryGB = 16
		#Int? diskGB_buffer = 10
		#Int diskGB = ceil(size(fastq1, 'G')*2 + size(fastq2, 'G') + diskGB_buffer)
        Int? diskGB = 100
	}

	command <<<
		# ~{gatk_path} samtools index ~{normal_bam} && \
		set -euxo pipefail
		~{gatk_path} Mutect2 \
		-R ~{ref_fasta} \
		-I ~{normal_bam}  \
		-tumor ~{sample_name} \
		-O "~{sample_name}.vcf"
		# --germline-resource \ # File germline vcf
	>>>

	runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}

	output {
		File output_vcf = "~{sample_name}.vcf"
	}
}

task CreateCommand1 {
    input {
        # Array[File] normal_bams
        # Changed from String file_path so that the file is localized in cromwell_root first.
        File file_path
    }
	command <<<
		set -euxo pipefail
		echo "--vcfs ~{file_path}"
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
        # Array[File] normal_bams
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



task CreateSomaticPanelOfNormals {
	input {
        String docker_path
        String gatk_path
        Array[File] normal_vcfs # localizes files for task
        String arguments
		Int? preemptible = 1
		Int? max_retries = 1
		Int? memoryGB = 16
		#Int? diskGB_buffer = 10
		#Int diskGB = ceil(size(fastq1, 'G')*2 + size(fastq2, 'G') + diskGB_buffer)
        Int? diskGB = 100
	}

	command <<<
		set -euxo pipefail
		~{gatk_path} CreateSomaticPanelOfNormals --duplicate-sample-strategy ALLOW_ALL --vcfs ~{sep=' --vcfs ' normal_vcfs} -O pon.vcf.gz
	>>>

	runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}

	output {
		File pon_vcf = "pon.vcf.gz"
	}
}