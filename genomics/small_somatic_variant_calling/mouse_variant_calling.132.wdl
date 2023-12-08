## About:
##	
## This WDL workflow is the second workflow in accepts the preprocessed BAM files produced by the Mouse WES/WGS Pipeline Preprocess   
## workflow and performs variant calling.
## 
## MAINTAINER (please contact for inquiries): 
## Neil Ruthen
## Associate Computational Biologist
## Translational Immunogenomics Lab 
## Dana-Farber Cancer Institute 
## neilg_ruthen@dfci.harvard.edu 
##
## NOTE: Please see the description below and the comments in the workflow definition for details on the usage. 
##
## Workflow defaults:
##
## - Allocated memory = 16 GB
##
## - Attempts with preemptible CPUs = 3
##
## - Supplied reference file (via mouse_preprocess workflow): Ensembl 102 GRCm38 primary assembly build
##
## Workflow usage/recommendations:
##
## - You will need provide a reference genome fasta file, index file, and dictionary as inputs to run this workflow.
##
## - You will also need a panel of normals variant calling file (e.g. pon.vcf file from mouse_pon workflow).
##
## - For WES samples, you will need a bed file containing the genomic intervals to scan over for variants. The 
##   genomic intervals file is usually a list of target regions provided by the sequencing company (e.g. Agilent). 
##   However, one can also use a list of exonic regions for the Ensembl 102 GRCm38 primary assembly build.
##
## - To output the union and/or the intersection of variant calls from Mutect2 and Strelka2, please
##   set the "consensus" boolean to "true". We recommend using the consensus_vcf for secondary analyses
##   (see the final comment at the end of the workflow definition for a list of outputs).
##
## - Set the "disk_multiplier" to 2 or 3. This ensures that the "disk_size" parameter exceeds (2x or 3x) the 
##   size of your largest sample file. 
##
## - WGS files with >30x coverage will not be processed by Mutect2 within the 24 hour limit of a preemptible CPU. 
##   In such cases, we recommend setting the preemptible parameter to 0 for the Mutect2 task.
##
## - To output Strelka2 indels in a separate MAF file, set the "separate_indels" boolean to "true".
##
## - Set the "reformat_intervals" parameter to "true" if the dictionary and intervals file have different 
##   contig formatting (Mutect2 will throw an error if they do). See comment above the ReformatBed task. 
##
## - Set the "reformat_pon" to "true" only if the samples are aligned with mm10 rather than GRCm38 and if you are
##   also using the default PON for this pipeline.
##
## Software versions:
## - GATK 4.1.8.1
## - Strelka2 2.9.3
## - Samtools-1.3.1
## - BCFtools-1.9.1
## - VCF2MAF 1.6.19

version 1.0 

workflow VariantCalling {
	input {
        String tumor_name
        String normal_name
        String pon_name
        Boolean reformat_intervals
        Boolean reformat_pon
        Boolean mutect
        Boolean strelka
        Boolean consensus
        Boolean separate_indels
        File ref_fasta
        File ref_index
        File ref_dict 
        File pon_vcf
        File normal_bam
        File tumor_bam
        Int disk_multiplier
        Int additional_space = 32
        Int disk_size = ceil(size(tumor_bam, "GB"))*disk_multiplier + additional_space
        String gatk_docker = "broadinstitute/gatk:4.1.8.1"
        String vcf2maf_docker = "nruthen/vcf2maf:0.0.2"
        String gatk_path = "/gatk/gatk"
        String samtools_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        String samtools_path = "/usr/gitc/samtools-1.3.1/"
        String curl_docker = "ellerbrock/alpine-bash-curl-ssl:0.3.0"
        String strelka_docker = "obenauflab/strelka:latest"
    	String bedtools_docker = "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
        String bcftools_docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
        String pysam_docker = "nruthen/reformvcf:0.0.1"
	}
	
    # Indexs tumor bam file with samtools
    call SamtoolsIndex as TumorIndex {
    	input:
    		docker_path = samtools_docker,
    		samtools_path = samtools_path,
    		sample_name = tumor_name,
    		sample_bam = tumor_bam,
            diskGB = disk_size
    }
    
    # Indexs normal bam file with samtools
    call SamtoolsIndex as NormalIndex {
    	input:
    		docker_path = samtools_docker,
    		samtools_path = samtools_path,
    		sample_name = normal_name,
    		sample_bam = normal_bam,
            diskGB = disk_size
    }    

    # Validates that tumor bams are free from improper formatting, incorrect alignment, corruption, etc.
    call ValidateBAM as ValidateTumor {
    	input: 
			docker_path = gatk_docker,
			gatk_path = gatk_path,
			sample_bam = TumorIndex.output_bam,
            diskGB = disk_size
    }
    
    # Validates that normal bams are free from improper formatting, incorrect alignment, corruption, etc.
    call ValidateBAM as ValidateNormal {
    	input: 
			docker_path = gatk_docker,
			gatk_path = gatk_path,
			sample_bam = NormalIndex.output_bam,
            diskGB = disk_size
    }
    
    # Optional task to remove the "chr" prefix from the chromosome column of .bed files 
    # (in order to match the expected GATK-style .bed formats) 
	if(reformat_intervals) {
    	call ReformatBed {
    		input:
            	reformat_intervals=reformat_intervals
    	}
    }

    # Optional task to add the "chr" prefix from the chromosome column of PON .vcf file
	if(reformat_pon) {
    	call ReformatPON {
    		input: 
        		pon_vcf = pon_vcf
    	}
    }
    
    # Performs variant calling with Manta (for detecting small indel candidates) and Strelka2 on the preprocessed bam files.
	call Strelka2 {
		input:
			docker_path = strelka_docker,
			ref_fasta = ref_fasta,
			ref_index = ref_index,
			ref_dict = ref_dict,
			tumor_bam = TumorIndex.output_bam,
			tumor_index = TumorIndex.output_bai,
			normal_bam = NormalIndex.output_bam,
			normal_index = NormalIndex.output_bai,
            diskGB = disk_size
	}
	
    # Sub-workflow to call variants with Mutect2 as well with Strelka2 + Manta, and take the union 
    # and intersection of the resulting variant calls, and annotate them with Ensembl VEP (recommended)
	if (consensus) {
	  # Performs variant calling with Mutect2 on the preprocessed bam files.
      call Mutect2 {
          input:
              docker_path = gatk_docker,
              gatk_path = gatk_path,
              sample_name = tumor_name, 
              ref_fasta = ref_fasta,
              ref_index = ref_index,
              ref_dict = ref_dict,
              tumor_bam = TumorIndex.output_bam,
              tumor_index = TumorIndex.output_bai,
              normal_bam = NormalIndex.output_bam,
              normal_index = NormalIndex.output_bai,
              normal_name = normal_name,
              reformatted_intervals = ReformatBed.intervals_bed,
              pon_vcf = pon_vcf,
              diskGB = disk_size
      }

	  # Performs filtering of raw variant calling with Mutect2.
      call FilterMutectCalls {
          input:
              docker_path = gatk_docker,
              gatk_path = gatk_path,
              sample_name = tumor_name, 
              ref_fasta = ref_fasta,
              ref_index = ref_index,
              ref_dict = ref_dict,
              mutect2_calls = Mutect2.output_vcf,
              mutect2_calls_stats = Mutect2.output_vcf_stats,
              diskGB = disk_size
      }

	  # Takes the consensus (intersection) of the variant calls from Strelka2 and Mutect2.
      call ConsensusVariantCalling {
          input:
              docker_path = bcftools_docker,
              ref_fasta = ref_fasta,
              ref_index = ref_index,
              ref_dict = ref_dict,
              indels_vcf = Strelka2.indels_vcf,
              snvs_vcf = Strelka2.snvs_vcf,
              mutect_vcf = FilterMutectCalls.output_vcf,
              normal_name = normal_name,
              tumor_name = tumor_name
      }
      
      # Adds caller annotation (whether a mutation was called by strelka2, mutect2, or both) to 
      # the FORMAT fields of the outputs of the consensus calling task. Also reformats the  
      # variants uniquely called by Strelka2 so that they include genotypes in the FORMAT field. 
      call ReformatVCFs {
          input:
              docker_path = pysam_docker,
              unique_mutect_vcf = ConsensusVariantCalling.unique_mutect_vcf,
              unique_strelka_vcf = ConsensusVariantCalling.unique_strelka_vcf,
              intersect_vcf = ConsensusVariantCalling.intersect_vcf,
              normal_name = normal_name,
              tumor_name = tumor_name
      }

      # Takes the union of the reformatted variant calls from Strelka2 and Mutect2.
      call UnionizeVCFs {
          input:
              docker_path = bcftools_docker,
              unique_mutect_vcf = ReformatVCFs.mutect_vcf,
              unique_strelka_vcf = ReformatVCFs.strelka_vcf,
              intersect_vcf = ReformatVCFs.consensus_vcf,
              normal_name = normal_name,
              tumor_name = tumor_name
      }
	  
      # Two custom docker images (parent image: nruthen/ensembl-vep, child image: nruthen/vcf2maf)
	  # provide support for running vcf2maf with VEP annotation tools for murine samples. Annotates 
      # all variant calls (union of variant calls) with Ensembl VEP and formats them into a MAF.
      call VCF2MAF as MAFForAll {
          input:
              docker_path = vcf2maf_docker,
              sample_name = tumor_name,
              normal_name = normal_name,
              ref_fasta = ref_fasta,
              ref_index = ref_index,
              ref_dict = ref_dict,
              sample_vcf = UnionizeVCFs.union_vcf
      }
	}

	# Sub-workflow for outputting Strelka2 indels in a separate MAF.
    if (separate_indels) {
      call ReformatSingleVCF {
          input:
              docker_path = pysam_docker,
              unique_vcf = Strelka2.indels_vcf,
              normal_name = normal_name,
              tumor_name = tumor_name
      }

      # Two custom docker images (parent image: nruthen/ensembl-vep, child image: nruthen/vcf2maf)
	  # provide support for running vcf2maf with VEP annotation tools for murine samples. Annotates 
      # Strelka2 indels with Ensembl VEP and formats them into a MAF.
      call VCF2MAF as MAFForSingle {
          input:
              docker_path = vcf2maf_docker,
              sample_name = tumor_name,
              normal_name = normal_name,
              ref_fasta = ref_fasta,
              ref_index = ref_index,
              ref_dict = ref_dict,
              sample_vcf = ReformatSingleVCF.reform_vcf
      }
    }
	
    # Outputs may include (in the order shown below) the .bam file validation summaries, 
    # variant calls unique to Mutect2, variant calls unique to Strelka2, the intersection
    # of variant calls from Mutect2 and Strelka2, a MAF file including all variant calls,
    # vcf of just the Strelka2 indels, a MAF of just the Strelka2 indels, and a vcf with
    # all variant calls from Strelka2 and Mutect2. 
    output {
    	File validation_tumor_bam = ValidateTumor.validation_summary
        File validation_normal_bam = ValidateTumor.validation_summary
        File? raw_mutect2_vcf = Mutect2.output_vcf
    	File? mutect2_vcf = ReformatVCFs.mutect_vcf
        File? strelka_vcf = ReformatVCFs.strelka_vcf
        File? consensus_vcf = ReformatVCFs.consensus_vcf
        File? output_maf = MAFForAll.output_maf
        File? single_vcf = ReformatSingleVCF.reform_vcf
        File? single_maf = MAFForSingle.output_maf
        File? union_vcf = UnionizeVCFs.union_vcf
    }
}

task Mutect2 {
	input {
        String docker_path
        String gatk_path
        String sample_name
        String normal_name
        File normal_bam
        File tumor_bam
        File tumor_index
        File normal_index
        File ref_fasta
        File ref_index
        File ref_dict
        File? reformatted_intervals
        File? unchanged_intervals
        File pon_vcf
		Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int diskGB
	}

	command <<<
		mv ~{ref_dict} $(dirname ~{ref_fasta})
		mv ~{ref_index} $(dirname ~{ref_fasta})
		~{gatk_path} IndexFeatureFile \
			-I ~{pon_vcf} && \
		if [ -f "~{reformatted_intervals}" ];
		then ~{gatk_path} Mutect2 \
			-R ~{ref_fasta} \
			-I ~{tumor_bam} \
			-I ~{normal_bam} \
			-tumor ~{sample_name} \
			-normal ~{normal_name} \
			-L ~{reformatted_intervals} \
			--panel-of-normals ~{pon_vcf} \
			-O "~{sample_name}.vcf";
		elif [ -f "~{unchanged_intervals}" ];
		then ~{gatk_path} Mutect2 \
			-R ~{ref_fasta} \
			-I ~{tumor_bam} \
			-I ~{normal_bam} \
			-tumor ~{sample_name} \
			-normal ~{normal_name} \
			-L ~{unchanged_intervals} \
			--panel-of-normals ~{pon_vcf} \
			-O "~{sample_name}.vcf";
		else ~{gatk_path} Mutect2 \
			-R ~{ref_fasta} \
			-I ~{tumor_bam} \
			-I ~{normal_bam} \
			-tumor ~{sample_name} \
			-normal ~{normal_name} \
			--panel-of-normals ~{pon_vcf} \
			-O "~{sample_name}.vcf";
		fi
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
        File output_vcf_stats = "~{sample_name}.vcf.stats"
	}
}

task FilterMutectCalls {
	input {
        String docker_path
        String gatk_path
        String sample_name
        File mutect2_calls
        File mutect2_calls_stats
        File ref_fasta
        File ref_index
        File ref_dict
		Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int diskGB
	}

	command <<<
		mv ~{ref_dict} $(dirname ~{ref_fasta})
		mv ~{ref_index} $(dirname ~{ref_fasta})
		~{gatk_path} FilterMutectCalls \
			-V ~{mutect2_calls} \
			-R ~{ref_fasta} \
			-O ~{sample_name}_filtered.vcf
 
	>>>

	runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}

	output {
		File output_vcf = "~{sample_name}_filtered.vcf"
	}
}

task VCF2MAF {
	input {
        String docker_path
        String vep_path = "/opt/vep/src/ensembl-vep"
        String vep_data = "/opt/vep/.vep/"
        String species = "mus_musculus"
        String version = "102"
        String build = "GRCm38"
        String filter_vcf = "0" 
        String sample_name
        String normal_name
        String? info_fields = "TLOD,NLOD,QSS_NT,QSI_NT,MQ,MBQ,MMQ"
        String? format_fields = "GT,AF,F1R2,F2R1,SB,CALLER"
        File sample_vcf
        File ref_fasta
        File ref_index
        File ref_dict
		Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 32
        Int? diskGB = 100
	}
	command <<<
		mv ~{ref_dict} $(dirname ~{ref_fasta})
		mv ~{ref_index} $(dirname ~{ref_fasta})
		perl /opt/vcf2maf.pl --input-vcf ~{sample_vcf} --output-maf ~{sample_name}.maf --ref-fasta ~{ref_fasta} --tumor-id ~{sample_name} --normal-id ~{normal_name} --vep-path ~{vep_path} --species ~{species} --cache-version ~{version} --ncbi-build ~{build} --vep-data ~{vep_data} --retain-fmt ~{format_fields} --retain-info ~{info_fields}
		MAF_FILE="~{sample_name}.maf"
		OUTPUT_FILE="~{sample_name}_complete.maf"
		HEADER=$(sed -n $'s/\t/\\\n/gp' ${MAF_FILE}| grep -n 'Hugo_Symbol' | cut -d: -f1)
		T_AF=$(sed -n $'s/\t/\\\n/gp' ${MAF_FILE}| grep -nx 't_AF' | cut -d: -f1)
		T_ALT=$(sed -n $'s/\t/\\\n/gp' ${MAF_FILE}| grep -nx 't_alt_count' | cut -d: -f1)
		T_DEPTH=$(sed -n $'s/\t/\\\n/gp' ${MAF_FILE}| grep -nx 't_depth' | cut -d: -f1)
		N_AF=$(sed -n $'s/\t/\\\n/gp' ${MAF_FILE}| grep -nx 'n_AF' | cut -d: -f1)
		N_ALT=$(sed -n $'s/\t/\\\n/gp' ${MAF_FILE}| grep -nx 'n_alt_count' | cut -d: -f1)
		N_DEPTH=$(sed -n $'s/\t/\\\n/gp' ${MAF_FILE}| grep -nx 'n_depth' | cut -d: -f1)
		head -n ${HEADER} ${MAF_FILE} >> ${OUTPUT_FILE}
		eval "awk -F '\t' 'BEGIN{OFS=\"\t\";} {if (NR>${HEADER}) {if (length(\$${T_AF})==0){\$${T_AF}=\$${T_ALT}/\$${T_DEPTH}; if (\$${N_DEPTH}==0){\$${N_AF}=(\$${N_ALT}+1)/(\$${N_DEPTH}+1)} else {\$${N_AF}=(\$${N_ALT}+1)/(\$${N_DEPTH})}; print \$0} else {print \$0}}}'" ${MAF_FILE} >> ${OUTPUT_FILE}
	>>>

	runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}

	output {
		File output_maf = "~{sample_name}_complete.maf"
	}
}

task SamtoolsIndex {
	input {
		String docker_path
        String sample_name
        String samtools_path
		File sample_bam
        Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int diskGB
    }
	command <<<
		samtools index ~{sample_bam} ~{sample_name}.bam.bai && \
		mv ~{sample_bam} ~{sample_name}.bam
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

task ReformatPON {
	input {
		File pon_vcf
        Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 4
        Int? diskGB = 5
    }
	command <<<
		apt-get update -y && \
		apt-get install -y gawk && \
		gunzip -c ~{pon_vcf} | cat | awk '{print(gensub(/^(2[0-2]|1[0-9]|[1-9]|[XY])/, "chr\\1", "g"))}' | gzip -c > "reformatted_pon.vcf.gz"
	>>>
    runtime {
		docker : "ubuntu:latest"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File reformatted_pon = "reformatted_pon.vcf.gz"
	}
}

task ReformatBed {
	input {
		File? intervals
        Boolean reformat_intervals
        Int? preemptible = 1
		Int? max_retries = 1
		Int? memoryGB = 4
        Int? diskGB = 5
    }
	command <<<
		awk '{gsub("chr", "")} 1' ~{intervals} > "intervals.bed"
	>>>
    runtime {
		docker : "ubuntu:latest"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File intervals_bed = "intervals.bed"
	}
}

task ValidateBAM {
	input {
		File sample_bam
        String gatk_path
        String docker_path
        Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int diskGB
    }
	command <<<
		~{gatk_path} ValidateSamFile \
			-I ~{sample_bam} \
			-MODE SUMMARY
	>>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File validation_summary = stdout()
	}
}

task Strelka2 {
    input {
        String docker_path
        String? manta_install_path = "/manta-1.5.0.centos6_x86_64"
        String? manta_analysis_path = "/manta-1.5.0.centos6_x86_64/bin"
        String? strelka_install_path = "/strelka-2.9.3.centos6_x86_64"
        String? strelka_analysis_path = "/strelka-2.9.3.centos6_x86_64/bin"
        File tumor_bam
        File normal_bam
        File tumor_index
        File normal_index
        File ref_fasta
        File ref_index
        File ref_dict
        Int? preemptible = 3
        Int? max_retries = 1
        Int? cores = 4
        Int? memoryGB = 24
        Int? allocation = memoryGB-4
        Int diskGB
    }

    command <<<
		echo "Running Manta!" && \
		~{manta_install_path}/bin/configManta.py \
		--normalBam ~{normal_bam} \
		--tumorBam ~{tumor_bam} \
		--referenceFasta ~{ref_fasta} \
		--runDir ~{manta_analysis_path} && \
		~{manta_analysis_path}/runWorkflow.py -m local -j ~{cores} -g ~{allocation} && \
		echo "Running Strelka!" && \
		~{strelka_install_path}/bin/configureStrelkaSomaticWorkflow.py \
		--normalBam ~{normal_bam} \
		--tumorBam ~{tumor_bam} \
		--referenceFasta ~{ref_fasta} \
		--indelCandidates ~{manta_analysis_path}/results/variants/candidateSmallIndels.vcf.gz \
		--runDir ~{strelka_analysis_path} && \
		~{strelka_analysis_path}/runWorkflow.py -h && \
		~{strelka_analysis_path}/runWorkflow.py -m local -j ~{cores} -g ~{allocation} && \
		mv ~{strelka_analysis_path}/results/variants/*.vcf.gz /cromwell_root
    >>>

    runtime {
        docker : "~{docker_path}" 
        disks: "local-disk ${diskGB} HDD"
        cpus: "~{cores}"
        preemptible: "~{preemptible}"
        memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
    }

    output {
        File snvs_vcf = "somatic.snvs.vcf.gz"
        File indels_vcf = "somatic.indels.vcf.gz"
    }
}

task BedtoolsConsensus {
    input {
        String docker_path
        String tumor_name
        File strelka_vcf
        File mutect2_vcf
        Int? preemptible = 3
        Int? max_retries = 1
        Int? memoryGB = 16
        Int? diskGB = 100
    }

    command <<<
        bedtools intersect -wa -wb -header -a ~{strelka_vcf} -b ~{mutect2_vcf} > "/cromwell_root/~{tumor_name}_consensus.vcf"             
    >>>

    runtime {
        docker : "~{docker_path}" 
        disks: "local-disk ${diskGB} HDD"
        preemptible: "~{preemptible}"
        memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
    }

    output {
        File consensus_vcf = "~{tumor_name}_consensus.vcf"
    }
}

task ConsensusVariantCalling {
	input {
		String docker_path
		File indels_vcf
        File snvs_vcf
        File mutect_vcf
        File ref_fasta
        File ref_index
        File ref_dict
        String normal_name
        String tumor_name
        String snvs_basename = basename(snvs_vcf)
        String indels_basename = basename(indels_vcf)
        String mutect_basename = basename(mutect_vcf)
        Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		bcftools sort -Oz -o "sorted_~{snvs_basename}" ~{snvs_vcf} && \
		bcftools sort -Oz -o "sorted_~{indels_basename}" ~{indels_vcf} && \
		bcftools sort -Oz -o "sorted_~{mutect_basename}" ~{mutect_vcf} && \
		bcftools index --force --tbi "sorted_~{snvs_basename}" && \
		bcftools index --force --tbi "sorted_~{indels_basename}" && \
		bcftools concat "sorted_~{snvs_basename}" "sorted_~{indels_basename}" -Ov --allow-overlaps -o combined_snvs_indels.vcf && \
		bcftools norm -m-any "combined_snvs_indels.vcf" | bcftools norm -Ov --check-ref w -f ~{ref_fasta} > "normalized_combined_snvs_indels.vcf" && \
		bcftools norm -m-any "sorted_~{mutect_basename}" | bcftools norm -Oz --check-ref w -f ~{ref_fasta} > "normalized_~{mutect_basename}" && \
		echo "TUMOR ~{tumor_name}" >> sample_names.txt && \
		echo "NORMAL ~{normal_name}" >> sample_names.txt && \
		bcftools reheader -s sample_names.txt -o "reheaded_combined_snvs_indels.vcf" "normalized_combined_snvs_indels.vcf" && \
		echo "~{tumor_name}" > sample_names.txt && \
		echo "~{normal_name}" >> sample_names.txt && \
		bcftools view -S sample_names.txt -o "recombined_snvs_indels.vcf" "reheaded_combined_snvs_indels.vcf"&& \
		bcftools view -Oz -o "recombined_snvs_indels.vcf.gz" "recombined_snvs_indels.vcf" && \
		bcftools view -S sample_names.txt -Oz -o "normalized_~{mutect_basename}.gz" "normalized_~{mutect_basename}" && \
		bcftools index --force --tbi "recombined_snvs_indels.vcf.gz" && \
		bcftools index --force --tbi "normalized_~{mutect_basename}.gz" && \
		bcftools isec -c snps -c indels -p union "normalized_~{mutect_basename}.gz" "recombined_snvs_indels.vcf.gz"
	>>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File output_vcf = "combined_snvs_indels.vcf"
		File unique_mutect_vcf = "/cromwell_root/union/0000.vcf"
        File unique_strelka_vcf = "/cromwell_root/union/0001.vcf"
        File intersect_vcf = "/cromwell_root/union/0002.vcf"
	}
}

task UnionizeVCFs {
	input {
		String docker_path
		File unique_mutect_vcf
        File unique_strelka_vcf
        File intersect_vcf
        String normal_name
        String tumor_name
        String mutect_basename = basename(unique_mutect_vcf)
        String strelka_basename = basename(unique_strelka_vcf)
        String intersect_basename = basename(intersect_vcf)
        Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		bcftools view -Oz -o "/cromwell_root/~{mutect_basename}.gz" ~{unique_mutect_vcf} && \
		bcftools view -Oz -o "/cromwell_root/~{strelka_basename}.gz" ~{unique_strelka_vcf} && \
		bcftools view -Oz -o "/cromwell_root/~{intersect_basename}.gz" ~{intersect_vcf} && \
		bcftools index --force --tbi "/cromwell_root/~{mutect_basename}.gz" && \
		bcftools index --force --tbi "/cromwell_root/~{strelka_basename}.gz" && \
		bcftools index --force --tbi "/cromwell_root/~{intersect_basename}.gz" && \
		bcftools concat /cromwell_root/*.vcf.gz -Ov --allow-overlaps -o mutect_strelka_union.vcf
	>>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
        File union_vcf = "mutect_strelka_union.vcf"
	}
}

task ReformatVCFs {
	input {
		String docker_path
		File unique_mutect_vcf
        File unique_strelka_vcf
        File intersect_vcf
        String normal_name
        String tumor_name
        String mutect_basename = basename(unique_mutect_vcf)
        String strelka_basename = basename(unique_strelka_vcf)
        String intersect_basename = basename(intersect_vcf)
        Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		python3 /opt/reformat_strelka.py ~{unique_mutect_vcf} "/cromwell_root/reform_~{mutect_basename}" ~{normal_name} ~{tumor_name} "mutect2" && \
		python3 /opt/reformat_strelka.py ~{unique_strelka_vcf} "/cromwell_root/reform_~{strelka_basename}" ~{normal_name} ~{tumor_name} "strelka" && \
		sed -i 's/GT1/GT/g' "/cromwell_root/reform_~{strelka_basename}" && \
		python3 /opt/reformat_strelka.py ~{intersect_vcf} "/cromwell_root/reform_~{intersect_basename}" ~{normal_name} ~{tumor_name} "both"
	>>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
		File mutect_vcf = "reform_~{mutect_basename}"
        File strelka_vcf = "reform_~{strelka_basename}"
        File consensus_vcf = "reform_~{intersect_basename}"
	}
}


task ReformatSingleVCF {
	input {
		String docker_path
		File unique_vcf
        String normal_name
        String tumor_name
        String unique_basename = basename(unique_vcf, ".gz")
        String? caller = "strelka"
        Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		gzip -cd ~{unique_vcf} > "/cromwell_root/~{unique_basename}" && \
		sed -i 's/TUMOR/~{tumor_name}/g' "/cromwell_root/~{unique_basename}" && \
		sed -i 's/NORMAL/~{normal_name}/g' "/cromwell_root/~{unique_basename}" && \
		python3 /opt/reformat_strelka.py "/cromwell_root/~{unique_basename}" "/cromwell_root/reform_~{unique_basename}" ~{normal_name} ~{tumor_name} ~{caller} && \
		sed -i 's/GT1/GT/g' "/cromwell_root/reform_~{unique_basename}"
	>>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
        File reform_vcf = "reform_~{unique_basename}"
	}
}