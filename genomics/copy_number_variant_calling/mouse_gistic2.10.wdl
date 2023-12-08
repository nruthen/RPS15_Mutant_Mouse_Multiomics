version 1.0

workflow gistic2 {
	input {
		String gistic_docker = "nruthen/gistic2:latest"
        String group_name
        Array[File] segmentation_files
        File gistic_reference
	}

	call gistic2 {
    	input:
        	docker_path = gistic_docker,
        	group_name=group_name,
        	segmentation_files=segmentation_files,
        	gistic_reference=gistic_reference
	}
    
	output {
		Array[File] gistic_outputs = gistic2.output_files
        File merged_segmentation = gistic2.merged_seg_file
    }
}
task gistic2 {
    input {
        String group_name
        String? output_dir = "/cromwell_root/gistic_output"
        String docker_path
        String gistic2_path = "/home/gistic"
        Array[File] segmentation_files
        File gistic_reference
        Int? maxspace_between_segs = 10000
        Float? amp_threshold = 0.1
        Float? del_threshold = 0.1
        # Boolean argument to use gene-level algorithm (1 = true, 0 = false)
        Int? use_gene_level = 1
        # Set the verbosity level for the execution log, 0 is the lowest level {0,10,20,30} 
        Int? verbosity_level = 10
        # Boolean arguments to reduce memory/disk usage (1 = true, 0 = false)
        Int? use_less_mem = 0 
        Int? use_less_disk = 0
        # Boolean argument to remove sex chromosomes before analysis (1 = true, 0 = false)
        Int? remove_xy = 0
        Float? cap = 1.5
        # Boolean argument to perform identify arm-level sCNAs (1 = true, 0 = false)
        Int? arm_level_analysis = 1
        # Fraction of the chromosome affected by a sCNA required to be classified as an arm-level event
        Float? arm_threshold = 0.5 
        # Confidence level for classifying driver mutations
        Float? confidence_level = 0.99
        Int? max_seg = 10000
        Float? res = 0.05 
        Int? do_arbitration = 1
        Int? arm_peel = 1
        Int? min_seg = 4 
		Int? preemptible = 3
		Int? max_retries = 1
		Int? memoryGB = 32
        Int disk_size = 100
    }
    command <<<
		sed -i.bak '/^@/d' ~{sep=" " segmentation_files} && \
		sed -i.bak '/^CONTIG/d' ~{sep=" " segmentation_files} && \
		seg_files=(~{sep=" " segmentation_files}) && \
		for seg in ${seg_files[@]}; do
			seg_basename="$(basename $seg T.called.seg)"
			awk -v var1="$seg_basename" 'BEGIN { OFS = "\t" } {print var1,$0}' $seg | awk 'BEGIN { OFS = "\t" } NF{NF--};1' | awk '$2 ~ /^([0-9]|X|Y)/' > /cromwell_root/temp_"$seg_basename".seg
			cat /cromwell_root/temp_"$seg_basename".seg > $seg
		done
		cat ~{sep=" " segmentation_files} > "/cromwell_root/segmentation_file.txt" && \
		mkdir ~{output_dir} && \
		cd ~{gistic2_path} && \
		~{gistic2_path}/gistic2 -b ~{output_dir} \
			-seg "/cromwell_root/segmentation_file.txt" \
			-refgene ~{gistic_reference} \
			-maxspace ~{maxspace_between_segs} -ta ~{amp_threshold} -td ~{del_threshold} -fname ~{group_name} \
			-rx ~{remove_xy} -cap ~{cap} -broad ~{arm_level_analysis} -brlen ~{arm_threshold} -maxseg ~{max_seg} \
			-res ~{res} -conf ~{confidence_level} -genegistic ~{use_gene_level} -arb ~{do_arbitration} \
			-armpeel ~{arm_peel} -savegene 1 -smalldisk ~{use_less_disk} -smallmem ~{use_less_mem} \
			-gcm extreme -v ~{verbosity_level} -js ~{min_seg}
		cd "/cromwell_root"
    >>>
    runtime {
        docker : "~{docker_path}"
        preemptible: "${preemptible}"
        disks: "local-disk " + disk_size + " HDD"
        memory: "${memoryGB}GB"
        maxRetries: "${max_retries}"
    }
    output {
        Array[File] output_files  = glob("~{output_dir}/*")
        File merged_seg_file = "/cromwell_root/segmentation_file.txt"
    }
} 