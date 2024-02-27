This README file outlines how to perform stop codon read through analysis from ribosomal sequencing data.

Required Inputs:

Ribosome profiling datasets in Fastq format 
Genome assembly file in Fasta format (GRCm38.primary_assembly.genome.fa)
Ribosomal RNA (rRNA) sequence file in Fasta format (mouse_rRNA_tRNA_U1snRNA_MN537.fasta)
Transcript annotation file in genePred or GTF format (gencode.vM24.annotation.gtf)


Requirements for Compute Environment:

Linux high performance computing cluster
Perl 5.8
R 3.3
BEDtools 2.29.0
SAMtools 1.8
Read mapping software (such as Bowtie(Langmead and Salzberg, 2012) and Tophat(Kim et al., 2013))

Analysis Steps: 

1. Perform the RibORF1.0 analysis as described in the following GitHub repository:

https://github.com/zhejilab/RibORF
	
The resulting candidate ORFs (candidateORF.fa) and offset-corrected bam files (*offset.bam) are the inputs for the stop codon read through analysis.

2. With the RibORF pipeline outputs candidateORF.fa and candidateORF.genepred.txt, select for entries with 'canonical' regular expression.

3. Convert RibORF pipeline output canonical_candidateORF.genepred.txt to BED-12 file using the genePredToBed tool from UCSC Genome Browser's Utility Tools.

3. Subset annotation to include only those within candidate ORF's regions found by the RibORF pipeline:

bedtools intersect -a gencode.vM24.annotation.gtf -b canonical_candidateORF.genepred.bed > riborf_bed_12_file_of_canonical_orfs.bed

4. Subset gencode.vM24.annotation.gtf for entries with 'stop_codon' regular expression.

5. Create bed file with regions upstream and downstream of each stop codon identified above. 

python stopcodon_readthrough_upstream_downstream_bed_generator.py -bed riborf_bed_12_file_of_canonical_orfs.bed -sub_gtf stop_codon_gencode.vM24.annotation.gtf -out stopCodon_updownstream.bed 

6. Filling in $bam_file with each offset corrected bam file path and $bed_file with the stopCodon_updownstream.bed file path, along with appropriate values for the rest of the bash variables, run the following two commands in parallel for all offset corrected bam files:

samtools sort $bam_file -o ${output_dir}${bam_basename}".offset.sorted.bam" 

bedtools coverage -a $bed_file -b ${output_dir}${bam_basename}".offset.sorted.bam" -d -s -split >> ${output_dir}${bam_basename}"_perbase.txt"

7. Shift indices of coverage files to indices relative to the identified stop codons with:

python shift_coverage.py -f Het_sample1_4_perbase.txt,HET_sample2_5_perbase.txt,HET_sample3_6_perbase.txt -o ${output_dir}

8. Create periodicity plots and coverage matrices for region centered around the first nucleotide of each stop codon +/- 28 bp.

Rscript plot_stopCodon_readThrough.R

9. Perform comparison upstream vs stop codon, upstream vs downstream coverage counts across experimental conditions to determine differential ribosomal stalling and differential stop codon read through, respectively.

Rscript DeSeq2_Pipeline_Wrapper_StopCodon_ReadThrough.R







