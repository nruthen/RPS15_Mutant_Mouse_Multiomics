## About:
##	The following workflow is the initial subworkflow of the Mouse WES/WGS Analysis Pipeline. 
##	The subworkflow performs the standard preprocessing steps for analyzing WES/WGS data according 
## 	to GATK best practices: alignment via BWA MEM, marking duplicates with Picard, and base 
## 	recalibration with GATK. 

version 1.0 


workflow mouse_preprocess {
	input {
    	# When specifying inputs via Terra tables sampleName is "this.sample_id"
        String sampleName
        # mouse_preprocess requires paired-end reads 
        File read1
        File read2
        # mouse_preprocess requires an Ensembl reference fasta from the Ensembl FTP Download server (see bwaMem task)
        File ref_ftp_link = "ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
        # mouse_preprocess requires variant calling files (and their index files) containing known germiline variants
        File known_vcf
        File known_vcf_idx
        File known_vcf_tbi
        # Task runtime variables for bwa-mem
        Int cores 
        Int threads = cores*2
        # Task disk size should exceed (2x or 3x) the size of your largest sample file 
        Int disk_size
        # Task command components (paths to directories within docker containers, complex commands)
        String bwa_path = "/usr/gitc/"
        String gatk_path = "/gatk/gatk"
        String samtools_docker = "erictdawson/samtools:2020-Nov-19"
        String gatk_docker = "broadinstitute/gatk:4.1.8.1"
        String bwa_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        # If you would like to output reference files to reuse in later runs, set need_references = true. 
        Boolean need_references
      }
      
    # Perform alignment on paired fastqs with bwa-mem to produce an aligned .bam file (see task definitions below)
	call bwaMem {
    	input:
    		sampleName = sampleName,
    		read1 = read1,
    		read2 = read2,
			ref_ftp_link = ref_ftp_link, 
    		cores = cores,
            threads = threads,
            disk_size = disk_size,
    		bwa_path = bwa_path,
            docker_path = bwa_docker
    }

    # Index reference genome 
	call SamtoolsFaidx {
    	input:
    		ref_fasta = bwaMem.ref_fasta,
        	docker_path = samtools_docker
    }

    # Mark duplicate reads in .bam with GATK MarkDuplicates (see task definitions below)
    call MarkDuplicates {
    	input:
    		sampleName = sampleName,
    		sortedBAM = bwaMem.output_bam,
    		gatk_path = gatk_path,
    		disk_size = disk_size,
            docker_path = gatk_docker
    }

  	
    # Calculate base-quality scores for reads in .bam with GATK BaseRecalibrator, also creates reference dictionary (see task definitions below)
    call BaseRecalibrator {
    	input:
    		sampleName = sampleName,
    		sortedDedupBAM = MarkDuplicates.output_bam,
    		gatk_path = gatk_path,
    		ref_fasta = bwaMem.ref_fasta,
    		ref_index = SamtoolsFaidx.ref_index,
    		known_vcf = known_vcf,
    		known_vcf_idx = known_vcf_idx,
    		known_vcf_tbi = known_vcf_tbi,
    		disk_size = disk_size,
            docker_path = gatk_docker
    }
    
	# Deduplicate and recalibrate reads in .bam with GATK ApplyBQSR based on duplicate labeling and quality scoring (see task definitions below)
    call ApplyBQSR {
    	input:
    		sampleName = sampleName,
    		sortedDedupBAM = MarkDuplicates.output_bam,
    		gatk_path = gatk_path,
    		ref_fasta = bwaMem.ref_fasta,
    		ref_index = SamtoolsFaidx.ref_index,
    		ref_dict = BaseRecalibrator.ref_dict,
    		recal_table = BaseRecalibrator.output_table,
    		disk_size = disk_size,
            docker_path = gatk_docker
	}

    # Gathers output reference files to reuse in later runs
    if (need_references) {
    	call GatherReferences {
        	input:
    			ref_fasta = bwaMem.ref_fasta,
    			ref_index = SamtoolsFaidx.ref_index,
    			ref_dict = BaseRecalibrator.ref_dict
		}	
    }
	# Outputs that will be retained when execution is complete  
	output {
		File duplication_metrics = MarkDuplicates.duplicate_metrics
		File bqsr_report = BaseRecalibrator.output_table
		File analysis_ready_bam = ApplyBQSR.output_bam
        File? ref_fasta = GatherReferences.output_ref_fasta
		File? ref_index = GatherReferences.output_ref_index
		File? ref_dict = GatherReferences.output_ref_dict
	} 
}

task bwaMem {
	input {
		String sampleName
		File read1
		File read2
        Int cores
        Int threads
        String ref_ftp_link
		String? bwa_commandline = "bwa mem -M -R \"@RG\\tID:\"\"~{sampleName}\"\"\\tSM:\"\"~{sampleName}\"\"\\tPL:ILLUMINA\" -t ~{threads} -K 10000000" 
		String bwa_path
        String ref_basename = basename(ref_ftp_link)
        String ref_basename2 = basename(ref_ftp_link, ".gz")
        String docker_path
		Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
	}

	command <<<
		# set -euxo pipefail
		curl -O ~{ref_ftp_link} && \
		gunzip ~{ref_basename} && \
		~{bwa_path}bwa index ~{ref_basename2} && \
		~{bwa_path}~{bwa_commandline} ~{ref_basename2} ~{read1} ~{read2} | samtools view -S -b | samtools sort > ~{sampleName}_SORTED.bam
	>>>
	runtime {
    	# continueOnReturnCode allows the pipeline to push through "soft" fails.
        # Please check the log files for errors accumulated during the run. 
		docker : "~{docker_path}"
		preemptible: "${preemptible}"
		cpu: cores
        disks: "local-disk " + disk_size + " HDD"
		memory: "${memoryGB}GB"
		maxRetries: "${max_retries}"
        continueOnReturnCode: [0, 1]
	}
	output {
		File output_bam = "${sampleName}_SORTED.bam"
        File ref_fasta = "~{ref_basename2}"
	}
}

task MarkDuplicates {
	input {
		String sampleName
		File sortedBAM
		String gatk_path
        String docker_path
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
	}
	# How to specify multiple bam inputs: "--INPUT ~{sep=' --INPUT ' input_bams} \"
	#  1-1-Preprocessing-For-Variant-Discovery-HG38 uses the following before directly calling MarkDuplicates
	# ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G" \
	command <<<
		set -euxo pipefail
		~{gatk_path} MarkDuplicates -I ~{sortedBAM} -M ~{sampleName}_SORTED_DEDUP_METRIC.txt -O ~{sampleName}_SORTED_DEDUP.bam 
	>>>
	runtime {
    	# continueOnReturnCode allows the pipeline to push through "soft" fails.
        # Please check the log files for errors accumulated during the run. 
		docker : "~{docker_path}"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: "${preemptible}"
		memory: "${memoryGB}GB"
		maxRetries: "${max_retries}"
		continueOnReturnCode: [0, 1]
	}
	output {
		File output_bam = "~{sampleName}_SORTED_DEDUP.bam"
		File duplicate_metrics = "~{sampleName}_SORTED_DEDUP_METRIC.txt"
	}
}

task BaseRecalibrator {
	input {
		String sampleName
        String gatk_path
        String docker_path
        String ref_basename = basename(ref_fasta, ".fa")
		File sortedDedupBAM
		File ref_fasta
        File ref_index
		File known_vcf
        File known_vcf_idx
        File known_vcf_tbi
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
	}
	command <<<
		~{gatk_path} CreateSequenceDictionary -R ~{ref_fasta} -O "~{ref_basename}.dict" && \
		mv "~{ref_basename}.dict" $(dirname ~{ref_fasta}) && \
		mv ~{ref_index} $(dirname ~{ref_fasta}) && \
		~{gatk_path} BaseRecalibrator \
		-I ~{sortedDedupBAM} \
		-R ~{ref_fasta} \
		--known-sites ~{known_vcf} \
		-O "~{sampleName}_SORTED_DEDUP_RECAL_DATA.table"
		mv $(dirname ~{ref_fasta})/~{ref_basename}.dict /cromwell_root
	>>>
	runtime {
    	# continueOnReturnCode allows the pipeline to push through "soft" fails.
        # Please check the log files for errors accumulated during the run. 
		docker : "~{docker_path}"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: "${preemptible}"
		memory: "${memoryGB}GB"
		maxRetries: "${max_retries}"
		continueOnReturnCode: [0, 1]
	}
	output {
		File output_table = "~{sampleName}_SORTED_DEDUP_RECAL_DATA.table"
        File ref_dict = "~{ref_basename}.dict"
	}
}

task ApplyBQSR {
	input {
		String sampleName
        String gatk_path
        String docker_path
		File sortedDedupBAM
		File recal_table
		File ref_fasta
        File ref_index
        File ref_dict
        Int disk_size
        # Optional runtime variables (>= 16 GB RAM recommended)
        Int? preemptible = 3
        Int? memoryGB = 16
        Int? max_retries = 3
	}
	command <<<
		set -euxo pipefail
		mv ~{ref_dict} $(dirname ~{ref_fasta}) && \
		mv ~{ref_index} $(dirname ~{ref_fasta}) && \
		~{gatk_path} ApplyBQSR \
		-R ~{ref_fasta} \
		-I ~{sortedDedupBAM} \
		--bqsr-recal-file ~{recal_table} \
		-O ~{sampleName}_SORTED_DEDUP_RECAL.bam
	>>>
	runtime {
    	# continueOnReturnCode allows the pipeline to push through "soft" fails.
        # Please check the log files for errors accumulated during the run.
		docker : "~{docker_path}"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: "${preemptible}"
		memory: "${memoryGB}GB"
		maxRetries: "${max_retries}"
		continueOnReturnCode: [0, 1]
	}
	output {
		File output_bam = "~{sampleName}_SORTED_DEDUP_RECAL.bam"
	}
}

task SamtoolsFaidx {
	input {
		String docker_path
		File ref_fasta
        String sample_basename = basename(ref_fasta)
        Int? preemptible = 1
		Int? max_retries = 1
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		samtools faidx ~{ref_fasta} -o "~{sample_basename}.fai"
	>>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File ref_index = "~{sample_basename}.fai"
	}
}

task GatherReferences {
	input {
		File ref_fasta
        File ref_index
        File ref_dict
        String ref_basename = basename(ref_fasta)
        Int? preemptible = 1
		Int? max_retries = 1
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		echo "Gathering reference fasta, index, and dictionary for output associated with: ~{ref_basename}"
	>>>
    runtime {
    	docker : "ubuntu:18.04"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
		File output_ref_fasta = "~{ref_fasta}"
    	File output_ref_index = "~{ref_index}"
        File output_ref_dict = "~{ref_dict}"
	}
}