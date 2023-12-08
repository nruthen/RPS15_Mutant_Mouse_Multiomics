## About:
##	This WDL workflow computes gene/transcipt expression from mouse RNAseq data and merges
##  the resulting expression data with the variant calling files from the Variant Calling Workflow.

## SNAPSHOT 90 is a temp version of the pipeline for single-end reads.

version 1.0 

workflow mouse_RNAseq {
	input {
        File read1
        File read2
        File salmon_index
        File sample_maf
        File ref_fasta
        File ref_index
        File ref_dict
        String sample_name
        String ensembl_version
        # "stranded" or "unstranded"
        String stranded
        # "inward", "outward", or "matching" 
        String relative_orientation
        # "forward", "reverse", or "na"
        String directionality 
        String curl_docker = "ellerbrock/alpine-bash-curl-ssl:0.3.0"
        String fastqc_docker = "biocontainers/fastqc:v0.11.9_cv7"
        String perl_docker = "nruthen/perl-ftp:0.0.1"
        String python_docker = "python:3.4"
        String salmon_docker = "combinelab/salmon:1.4.0"
        String salmontools_docker = "nruthen/salmon-tools:0.0.1"
        String refgenie_docker = "databio/refgenie:latest"
        String star_docker = "nruthen/star:2.7.7a"
        String pandas_docker = "amancevice/pandas:1.1.4"
        String bcftools_docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
        # String samtools_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        #String samtools_docker = "biocontainers/samtools:v1.9-4-deb_cv1"
		String samtools_docker = "erictdawson/samtools:2020-Nov-19"
        String rsem_docker = "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        #String samtools_path = "/usr/gitc/samtools-1.3.1/"
        Boolean need_aligned_bam
	}
    
  	call FASTQC as FASTQC_before {
    	input:
			docker_path = fastqc_docker,
        	read1 = read1, 
            read2 = read2
	}
    
    call getTranscriptome {
    	input:
        	docker_path = curl_docker
    }
    
    call Single_End_getLibType {
    	input:
        	docker_path = python_docker,
            stranded = stranded,
            relative_orientation = relative_orientation,
            directionality = directionality
    
	}

    #call SalmonTools {
    	#input:
        	#docker_path = salmontools_docker,
        	#txome_fasta = getTranscriptome.transcripts_ftp_file,
            #gtf_file = getTranscriptome.transcripts_gtf_file,
            #genome_fasta = ref_fasta
    #}

    call Single_End_Salmon {
    	input:
        	docker_path = salmon_docker,
        	read1 = read1,
            salmon_index = salmon_index,
            lib_type = Single_End_getLibType.lib_type
    }
    
    #call STAR {
    	#input:
        	#docker_path = star_docker,
        	#read1 = read1,
			#read2 = read2,
        	#genome_fasta = getTranscriptome.genome_file,
			#sjdbGTFFile = getTranscriptome.transcripts_gtf_file
    #}
    
	call STAR_Single_End {
    	input:
        	docker_path = star_docker,
        	read1 = read1,
        	genome_fasta = getTranscriptome.genome_file,
			sjdbGTFFile = getTranscriptome.transcripts_gtf_file
    }
    
	#call RSEM {
    	#input:
        	#docker_path = rsem_docker,
			#input_bam = STAR.transcriptome_bam[0],
            #ref_fasta = getTranscriptome.genome_file,
            #gtf_file = getTranscriptome.transcripts_gtf_file
    #}

	call RSEM {
    	input:
        	docker_path = rsem_docker,
			input_bam = STAR_Single_End.transcriptome_bam[0],
            ref_fasta = getTranscriptome.genome_file,
            gtf_file = getTranscriptome.transcripts_gtf_file
    }
    #call GenerateTargets {
            #input:
                #sample_maf = sample_maf,
                #docker_path = pandas_docker
    #}
    
	#call SamtoolsIndex {
    	#input:
    		#docker_path = samtools_docker,
    		#sample_bam = STAR.output_bam[0]
    #}
    
    
    #call GatherReferences {
        #input:
    		#ref_fasta = ref_fasta,
    		#ref_index = ref_index,
    		#ref_dict = ref_dict
	#}	

    #call Mpileup {
            #input:
                #docker_path = bcftools_docker,
                #sample_name = sample_name,
		        #sample_bam = SamtoolsIndex.output_bam,
                #sample_bai = SamtoolsIndex.output_bai,
				#ref_fasta = GatherReferences.output_ref_fasta,
				#ref_index = GatherReferences.output_ref_index,
				#ref_dict = GatherReferences.output_ref_dict,
                #target_regions = GenerateTargets.output_bed
    #}

    #call CoverageReport {
    	#input:
			#sample_maf = sample_maf,
			#docker_path = pandas_docker,
            #rna_vcf = Mpileup.output_vcf
    #}
    
    if (need_aligned_bam) {
    	call GatherAlignedBAM {
    		input:
        		sample_bam=STAR_Single_End.transcriptome_bam[0]
		}
    }
    
    #if (need_aligned_bam) {
    	#call GatherAlignedBAM {
    		#input:
        		#sample_bam=STAR.output_bam[0]
		#}
    #}
   
    output {
		Array[File] first_fastqc_report = FASTQC_before.fastqc_report
    	#File mpileup_vcf = Mpileup.output_vcf
        File salmon_expression_file = Single_End_Salmon.salmon_output
        File rsem_expression_file = RSEM.output_genes[0]
        File rsem_isoforms_file = RSEM.output_isoforms[0]
        Array[File] rsem_stats = RSEM.output_stats
        Array[File] star_logs = STAR_Single_End.log_files_out
    	#File? output_bam = GatherAlignedBAM.output_bam
    }
}

task getLibType {
	input {
    	String docker_path
    	String stranded
        String relative_orientation
        String directionality
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 4
        Int? diskGB = 1
    }
    
    command <<<
		python <<CODE
		strandcode = {
			"stranded" : "S",
			"unstranded" : "U"
		}
		orientcode = {
			"inward" : "I",
			"outward" : "O",
			"matching" : "M"
		}
		directcode = {
			"forward" : "F",
			"reverse" : "R",
			"na" : " "
		}
		print((orientcode["~{relative_orientation}"]+strandcode["~{stranded}"]+directcode["~{directionality}"]).strip())
		CODE
	>>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	String lib_type = read_string(stdout())
	}    
}

task Single_End_getLibType {
	input {
    	String docker_path
    	String stranded
        String relative_orientation
        String directionality
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 4
        Int? diskGB = 1
    }
    
    command <<<
		python <<CODE
		strandcode = {
			"stranded" : "S",
			"unstranded" : "U"
		}
		orientcode = {
			"inward" : "I",
			"outward" : "O",
			"matching" : "M",
			"na" : " "
		}
		directcode = {
			"forward" : "F",
			"reverse" : "R",
			"na" : " "
		}
		print((orientcode["~{relative_orientation}"]+strandcode["~{stranded}"]+directcode["~{directionality}"]).strip())
		CODE
	>>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	String lib_type = read_string(stdout())
	}    
}

task SalmonTools {
	input {
    	File txome_fasta
        File genome_fasta
        File gtf_file
        String? output_path = "/cromwell_root"
    	String docker_path
        Int? preemptible = 3
		Int? max_retries = 3
        Int? cores = 16
        Int threads = cores*2
		Int? memoryGB = 128
        Int? diskGB = 100
    }
    
    command <<<
		bash /opt/conda/bin/SalmonTools/scripts/generateDecoyTranscriptome.sh -a ~{gtf_file} -g ~{genome_fasta} -t ~{txome_fasta} -o ~{output_path} -j ~{threads}
	>>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
		cpu: "${cores}"
	}
    
	output {
    	File aware_txome = "gentrome.fa"
        File decoys = "decoys.txt"
	}    
}

task Salmon {
    input {
		File salmon_index
        String lib_type
        File read1
        File read2
        String docker_path
        String? out_dir = "quant"
        String? index_dir = "/cromwell_root/salmon_index/"
        String base_name = basename(salmon_index)
        Int? cores = 16
        Int threads = cores*2
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 128
        Int? diskGB = 200
    }

	command <<<
		# In case the PATH variable doesnt work, the salmon executable is located in: /home/salmon-1.4.0/bin
		# LIBTYPE needs to be defined as IU, MU, or OU. Need to ask Mohamed about the library prep. See the following link: https://salmon.readthedocs.io/en/latest/library_type.html
		# Apparently, the -i or "index" argument specifies the directory in which to store the salmon index file. For the index executable, -i is the directory you want to output to 
		# (I believe...). If the directory doesn't exist, the BuildSalmonIndex.cpp script will make the directory for you. This make cause problems as we may lose track of the output 
		# directory. For the quant executable, -i is the directory where your salmon index is stored and -o is the output directory...hopefully we can keep track of it. 
		# May want to specify -p / --threads for quantification.
		# quantmerge might be a useful option if we need to merge multiple reads (or duplicates)
		mkdir salmon_index
		mv ~{salmon_index} ~{index_dir}
		tar -zxvf ~{index_dir}~{base_name}
		cd ~{index_dir}
		ls
		cd /cromwell_root
		salmon quant -i ~{index_dir} -l ~{lib_type} -1 ~{read1} -2 ~{read2} --validateMappings -o ~{out_dir}
    >>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File salmon_output = "/cromwell_root/quant/quant.sf"
	}    
}


task Single_End_Salmon {
    input {
		File salmon_index
        String lib_type
        File read1
        String docker_path
        String? out_dir = "quant"
        String? index_dir = "/cromwell_root/salmon_index/"
        String base_name = basename(salmon_index)
        Int? cores = 16
        Int threads = cores*2
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 128
        Int? diskGB = 200
    }

	command <<<
		# In case the PATH variable doesnt work, the salmon executable is located in: /home/salmon-1.4.0/bin
		# LIBTYPE needs to be defined as IU, MU, or OU. Need to ask Mohamed about the library prep. See the following link: https://salmon.readthedocs.io/en/latest/library_type.html
		# Apparently, the -i or "index" argument specifies the directory in which to store the salmon index file. For the index executable, -i is the directory you want to output to 
		# (I believe...). If the directory doesn't exist, the BuildSalmonIndex.cpp script will make the directory for you. This make cause problems as we may lose track of the output 
		# directory. For the quant executable, -i is the directory where your salmon index is stored and -o is the output directory...hopefully we can keep track of it. 
		# May want to specify -p / --threads for quantification.
		# quantmerge might be a useful option if we need to merge multiple reads (or duplicates)
		mkdir salmon_index
		mv ~{salmon_index} ~{index_dir}
		tar -zxvf ~{index_dir}~{base_name}
		cd ~{index_dir}
		ls
		cd /cromwell_root
		salmon quant -i ~{index_dir} -l ~{lib_type} -r ~{read1} --validateMappings -o ~{out_dir}
    >>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File salmon_output = "/cromwell_root/quant/quant.sf"
	}    
}


task getTranscriptome {
	input {
    	String docker_path
        String? transcripts_ftp = "ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
        String? genome_ftp = "ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
        String? transcripts_gtf = "ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
        String? ftp_base = "Mus_musculus.GRCm38.cdna.all.fa.gz"
        String? gtf_base = "Mus_musculus.GRCm38.102.gtf.gz"
        String? genome_base = "Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 8
        Int? diskGB = 64
    }
    command <<<
		#curl -O ~{transcripts_ftp} &&\
		curl -O ~{transcripts_gtf} &&\
		curl -O ~{genome_ftp} &&\
		ls
    >>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
    }
    output {
    	File transcripts_gtf_file = "~{gtf_base}"
        File genome_file = "~{genome_base}"
    }
}

task prepTranscriptome {
    input {
        String ensembl_version
        String docker_path
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 16
        Int? diskGB = 100
    }

	command <<<
		perl5.32.0 <<CODE
		$FASTA_URL = "ftp://ftp.ensembl.org/pub/release-${ensembl_version}/fasta/mus_musculus/cdna/";
		if($FASTA_URL =~ /^ftp/i) {
			$FASTA_URL =~ m/(ftp:\/\/)?(.+?)\/(.+)/;
			$ftp = Net::FTP->new($2, Passive => 1) or die "ERROR: Could not connect to FTP host $2\n$@\n";
			$ftp->login("anonymous") or die "ERROR: Could not login as anonymous\n$@\n";
			$ftp->binary();

			foreach my ($sub(split /\//, $3)) {
				$ftp->cwd($sub) or die "ERROR: Could not change directory to $sub\n$@\n";
			}
			$file = grep {$_ =~ /.all.fa.gz/} $ftp->ls;
			$ftp->get($file) or die "get failed ", $ftp->message;
			$ftp->quit;
		}
		else {
			print "Must select a valid ensembl-vep ftp link to extract the appropriate transciptome fasta from the ftp server.\n"
		}
		CODE
	>>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	Array[File] ftp_file = glob("*.all.fa.gz")
        
	}  
}

task getTranscriptome2 {
	input {
    	String docker_path
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 16
        Int? diskGB = 100
    }
    command <<<
		#refgenie pull mm10_cdna/salmon_index
		#tar -cxvf salmon_index.tar.gz /cromwell_root/mm10_cdna/salmon_index
		curl -o "salmon_index.tgz" "http://refgenomes.databio.org/v2/asset/hg19_cdna/salmon_index/archive?tag=default" 
		ls
    >>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
    }
    output {
		File salmon_index = "salmon_index.tgz"
    }
}

task STAR {
		input {
    	String docker_path
        String? genomeDirPath = "/cromwell_root/star_genome"
        String? genome_base = "Mus_musculus.GRCm38.dna.primary_assembly.fa"
        String? gtf_base = "Mus_musculus.GRCm38.102.gtf"
        File genome_fasta
        File sjdbGTFFile
        File read1
        File read2
        Int? preemptible = 3
		Int? max_retries = 3
        Int? cores = 16
		Int? memoryGB = 512
        Int? diskGB = 300
    }
    command <<<
	mkdir star_genome
	gunzip ~{genome_fasta}
	gunzip ~{sjdbGTFFile}
	genomepath=$(dirname ~{genome_fasta})
	gtfpath=$(dirname ~{sjdbGTFFile})
	echo $genomepath
	echo $gtfpath
	STAR \
		--runThreadN ~{cores} \
		--runMode genomeGenerate \
		--genomeDir "/cromwell_root/star_genome" \
		--genomeFastaFiles "${genomepath}/~{genome_base}" \
		--sjdbGTFfile "${gtfpath}/~{gtf_base}" && \
	echo "Finished generating genome, now performing alignment!" && \
	STAR \
		--readFilesIn ~{read1} ~{read2} \
		--readFilesCommand "zcat" \
		--genomeDir "~{genomeDirPath}" \
		--runThreadN ~{cores} \
		--sjdbGTFfile "${gtfpath}/~{gtf_base}" \
		--outFilterMismatchNoverLmax 0.05 \
		--outFilterMatchNmin 15 \
		--outFilterScoreMinOverLread 0 \
		--outFilterMatchNminOverLread 0 \
		--outReadsUnmapped Fastx \
		--quantMode TranscriptomeSAM \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix ./ > star_log.txt 2>&1
    >>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
    }
    output {
		File output_log = "star_log.txt"
        
        Array[File] output_bam = glob("*Aligned.sortedByCoord.out.bam*")
        Array[File] transcriptome_bam = glob("*Aligned.toTranscriptome.out.bam*")
        Array[File] log_files_out = glob("*.out")
    }
}

task STAR_Single_End {
		input {
    	String docker_path
        String? genomeDirPath = "/cromwell_root/star_genome"
        String? genome_base = "Mus_musculus.GRCm38.dna.primary_assembly.fa"
        String? gtf_base = "Mus_musculus.GRCm38.102.gtf"
        File genome_fasta
        File sjdbGTFFile
        File read1
        Int? preemptible = 3
		Int? max_retries = 3
        Int? cores = 16
		Int? memoryGB = 512
        Int? diskGB = 300
    }
    command <<<
	mkdir star_genome
	gunzip ~{genome_fasta}
	gunzip ~{sjdbGTFFile}
	genomepath=$(dirname ~{genome_fasta})
	gtfpath=$(dirname ~{sjdbGTFFile})
	echo $genomepath
	echo $gtfpath
	STAR \
		--runThreadN ~{cores} \
		--runMode genomeGenerate \
		--genomeDir "/cromwell_root/star_genome" \
		--genomeFastaFiles "${genomepath}/~{genome_base}" \
		--sjdbGTFfile "${gtfpath}/~{gtf_base}" && \
	echo "Finished generating genome, now performing alignment!" && \
	STAR \
		--readFilesIn ~{read1} \
		--readFilesCommand "zcat" \
		--genomeDir "~{genomeDirPath}" \
		--runThreadN ~{cores} \
		--sjdbGTFfile "${gtfpath}/~{gtf_base}" \
		--outFilterMismatchNoverLmax 0.05 \
		--outFilterMatchNmin 15 \
		--outFilterScoreMinOverLread 0 \
		--outFilterMatchNminOverLread 0 \
		--outReadsUnmapped Fastx \
		--quantMode TranscriptomeSAM \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix ./ > star_log.txt 2>&1
    >>>
    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
    }
    output {
		File output_log = "star_log.txt"
        
        Array[File] output_bam = glob("*Aligned.sortedByCoord.out.bam*")
        Array[File] transcriptome_bam = glob("*Aligned.toTranscriptome.out.bam*")
        Array[File] log_files_out = glob("*.out")
    }
}

task GenerateTargets {
    input {
    	String docker_path
        File sample_maf
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 4
        Int? diskGB = 10
    }

	command <<<
		python <<CODE
		import pandas
		dataframe = pandas.read_csv("~{sample_maf}", delimiter='\t', header=1)
		dataframe = dataframe[(dataframe['Variant_Classification'] != "Silent") | (dataframe['Variant_Classification'] != "Intron")]
		dataframe = dataframe.loc[:, dataframe.columns.intersection(['Chromosome', 'Start_Position', 'End_Position'])]
		dataframe.to_csv('target_regions.tsv', sep='\t', header = False, index=False)
		CODE
	>>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}
    
	output {
    	File output_bed = "target_regions.tsv"
	}    
}

task Mpileup {
	input {
		String docker_path
        String sample_name
		File sample_bam
        File sample_bai
        File ref_fasta
        File ref_index
        File ref_dict
        File target_regions
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		#fasta_dir=$(dirname ~{ref_fasta})
		#mv ~{ref_dict} "$fasta_dir" && \
		#mv ~{ref_index} "$fasta_dir" && \
		bcftools mpileup -Ou -a INFO/AD,FORMAT/AD,FORMAT/SP -f ~{ref_fasta} --regions-file ~{target_regions} ~{sample_bam} | bcftools call -m -Ov -o ~{sample_name}_filtered.vcf 
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

task CoverageReport {
	input {
    	String docker_path
        File sample_maf
        File rna_vcf
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 4
        Int? diskGB = 10
    }

	command <<<
		python <<CODE
		import io
		import os
		import pandas
		dataframe1 = pandas.read_csv("~{sample_maf}", delimiter='\t', header=1)
		dataframe1 = dataframe1[(dataframe1['Variant_Classification'] != "Silent") | (dataframe1['Variant_Classification'] != "Intron")]
		print("Filtered MAF file.")
		with open("~{rna_vcf}", 'r') as f:
			lines = [l for l in f if not l.startswith("##")]
		dataframe2 = pandas.read_csv(io.StringIO(''.join(lines)), dtype={"#CHROM": str, "POS": int, "ID": str, "REF": str, "ALT": str,"QUAL": str, "FILTER": str, "INFO": str}, sep='\t').rename(columns={"#CHROM": "CHROM"})
		print("Imported VCF file as dataframe.")
		total_percent_captured = (dataframe2.size/dataframe1.size)*100
		print("Calculating capture.")
		position_counts = dataframe2["POS"].value_counts()
		print("Calculating number of reads per position.")
		percent_missense_captured = (dataframe2[dataframe2["POS"].isin(position_counts.index[position_counts.lt(2)])].size/dataframe1[dataframe1['Variant_Classification'] == "Missense_Mutation"].size)*100
		print("Percentage of total WES variants present in RNA: " + str(total_percent_captured))
		print("Percentage of WES missense mutations present in RNA: " + str(percent_missense_captured))
		CODE
	>>>

    runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
        continueOnReturnCode: [0, 1]
	}
    
	output {
    	File output_txt = stdout()
	}   
}

task SamtoolsIndex {
	input {
		String docker_path
		File sample_bam
        String sample_name = basename(sample_bam, ".bam")
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 16
		#Int? diskGB_buffer = 10
		#Int diskGB = ceil(size(fastq1, 'G')*2 + size(fastq2, 'G') + diskGB_buffer)
        Int? diskGB = 100
    }
	command <<<
		perl -v && \
		samtools --version && \
		#samtools reheader -c 'perl -pe "s/^(@SQ.*)(\tSN:)(\d+|X|Y|MT)(\s|\$)/\$1\$2Chr\$3\$4/"' ~{sample_bam} > ~{sample_name}_out.bam && \
		#samtools view -H ~{sample_name}_out.bam && \
		samtools view -H ~{sample_bam} && \
		samtools index ~{sample_bam} ~{sample_name}.bam.bai && \
		mv ~{sample_bam} ~{sample_name}.bam && \
		samtools idxstats ~{sample_name}.bam | cut -f 1 | head -3
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

task GatherReferences {
	input {
		File ref_fasta
        File ref_index
        File ref_dict
        String ref_basename = basename(ref_fasta)
        String index_basename = basename(ref_index)
        String dict_basename = basename(ref_dict)
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 16
        Int? diskGB = 100
    }
	command <<<
		echo "Gathering reference fasta, index, and dictionary for output associated with: ~{ref_basename}"
		mv ~{ref_fasta} "/cromwell_root/~{ref_basename}"
		mv ~{ref_index} "/cromwell_root/~{index_basename}"
		mv ~{ref_dict} "/cromwell_root/~{dict_basename}"
	>>>
    runtime {
    	docker : "ubuntu:18.04"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}

	output {
		File output_ref_fasta = "~{ref_basename}"
    	File output_ref_index = "~{index_basename}"
        File output_ref_dict = "~{dict_basename}"
	}
}

task GatherAlignedBAM{
	input {
		File sample_bam
        String sample_basename = basename(sample_bam)
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 16
        Int? diskGB = 32
    }
	command <<<
		echo "Outputting aligned BAM file from STAR."
		mv ~{sample_bam} "/cromwell_root/~{sample_basename}"
	>>>
    runtime {
    	docker : "ubuntu:18.04"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}

	output {
		File output_bam = "~{sample_basename}"
	}
}

task RSEM {
    input {
        File input_bam
        File gtf_file
        File ref_fasta 
        String ref_basename = basename(ref_fasta, ".fa.gz")
        String sample_basename = basename(input_bam, "bam")
        String gtf_basename = basename(gtf_file, ".gz")
        String? star_path = "/opt/STAR-2.7.9a/bin/Linux_x86_64_static"
        String docker_path
        Int? max_frag_length = 1000
        Int? cores = 4
        Int? threads = cores
        Int? preemptible = 3
		Int? max_retries = 3
		Int? memoryGB = 32
        Int? diskGB = 200
    }

	command <<<
		mkdir /cromwell_root/ref && \
		mkdir /cromwell_root/out && \
		gzip -cd ~{gtf_file} > "/cromwell_root/~{gtf_basename}" && \
		gzip -cd ~{ref_fasta} > "/cromwell_root/~{ref_basename}" && \
		rsem-prepare-reference --gtf "/cromwell_root/~{gtf_basename}" \
			--star --star-path ~{star_path} \
			--num-threads ~{threads} "/cromwell_root/~{ref_basename}" "/cromwell_root/ref/~{ref_basename}" && \
		rsem-calculate-expression --num-threads ~{threads} \
			--fragment-length-max ~{max_frag_length} \
			--no-bam-output --paired-end \
			--estimate-rspd --bam \
			~{input_bam} "/cromwell_root/ref/~{ref_basename}" "/cromwell_root/out/~{sample_basename}"
	>>>

	runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
        maxRetries: "~{max_retries}"
	}

	output {
		Array[File] output_genes = glob("/cromwell_root/out/~{sample_basename}.genes.results*")
		Array[File] output_isoforms = glob("/cromwell_root/out/~{sample_basename}.isoforms.results*")
        Array[File] output_stats = glob("/cromwell_root/out/~{sample_basename}.stat/*")
	}
}


task FASTQC {
	input {
		File read1
		File read2
        String docker_path
		Int? preemptible = 1
		Int? max_retries = 1
		Int? memoryGB = 32
		Int? cores = 8
        Int threads = cores
		Int? diskGB = 100
	}
	command <<<
		fastqc --threads ~{threads} --outdir `pwd` ~{read1}
		fastqc --threads ~{threads} --outdir `pwd` ~{read2}
	>>>
	runtime {
		docker : "~{docker_path}"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "~{preemptible}"
		memory: "~{memoryGB}GB"
		cpus: "~{cores}"
		maxRetries: "~{max_retries}"
	}
	output {
		Array[File] fastqc_report = glob("*fastqc.*")
	}
}