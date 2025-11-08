#! /bin/bash

# MOUSE version of pipeline
# Perl-5.8
# R-3.3
# samtools-1.8
# anaconda3-5.0.1
# bedtools 2.29.0

# Obtain demultiplexed fastq files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293841 as input
# Run pipeline separately for each pair of fastq files associated with a sample

set -e

############################################################################################

SAMPLE1=$1 # input GZ file (demultiplexed, trimmed, filtered for length and 5' nt, as in)
SAMPLE2=$2 # input GZ file (demultiplexed, trimmed, filtered for length and 5' nt, as in)
OUT=$3 # dir for outputs
NAME=$4 # name to prefix all files with
UMI_TOOLS_PATH=$5 # Path to UMI tools installation directory
BOWTIE_PATH=$6 # Path to Bowtie installation directory
STAR_PATH=$7 # Path to STAR aligner installation directory
FASTQC_PATH=$8 # Path to FASTQC installation directory
SAMTOOLS_PATH=$9 # Path to samtools installation directory
PERIODICTY_TOOL_PATH=$10 # Path to directory containing readDist_modified2.pl script from RibORF 1.0 toolkit
BOWTIE_REF=$11 # Path to file with blacklist sequences for rRNA, tRNA, and snRNA: "~/bowtie_index/mouse_rRNA_tRNA_U1snRNA_MN537"
GENE_PRED=$12 # "~/reference/Mouse/gencode.vM24.protein_coding.RibORF.genePred"
star_genome=$13 # STAR aligner index (with genome version that matches annotation files) "~/reference/Mouse/STAR_index_mouse"

mkdir -p $OUT/umi_trimmed/
mkdir -p $OUT/Bowtie_rRNA_Removed/ || echo "Bowtie_rRNA_Removed out already exists"

# MERGE PAIRED FASTQ FILES
SAMPLE1_BASE=$(basename "$SAMPLE1")
gunzip -cd $SAMPLE1 > $OUT/$SAMPLE1_BASE".fastq"

SAMPLE2_BASE=$(basename "$SAMPLE2")
gunzip -cd $SAMPLE2 > $OUT/$SAMPLE2_BASE".fastq"

cat $OUT/$SAMPLE1_BASE".fastq" $OUT/$SAMPLE2_BASE".fastq" | awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) > 20) {print header, seq, qheader, qseq}}' > $OUT/$NAME"_merged.fastq"

rm $OUT/$SAMPLE1_BASE".fastq"
rm $OUT/$SAMPLE2_BASE".fastq"

INPUT_FQGZ=$OUT/$NAME"_merged.fastq"

# Extract UMI
NO_UMI=$OUT/umi_trimmed/$NAME".fq"
$UMI_TOOLS_PATH/umi_tools extract \
--stdin=$INPUT_FQGZ --extract-method=string \
--bc-pattern=NNNNNNNNNNNN --3prime \
--log=$OUT/umi_trimmed/$NAME".log" --stdout $NO_UMI

rm $INPUT_FQGZ

# Output directory construction 
mkdir -p $OUT/Bowtie_rRNA_Removed/Aligned/ || echo "Aligned already exists"
mkdir -p $OUT/Bowtie_rRNA_Removed/Unaligned/ || echo "Unaligned already exists"
mkdir -p $OUT/Bowtie_rRNA_Removed/Stats/ || echo "Stats already exists"
mkdir -p $OUT/Bowtie_rRNA_Removed/STDOUT/ || echo "STDOUT already exists"

ALIGNED=$OUT/Bowtie_rRNA_Removed/Aligned/$NAME".aligned.out"
ALIGNED="${ALIGNED%%.*}.aligned"
UNALIGNED=$OUT/Bowtie_rRNA_Removed/Unaligned/$NAME".unaligned.out"
UNALIGNED="${UNALIGNED%%.*}.unaligned"
STATS=$OUT/Bowtie_rRNA_Removed/Stats/$NAME".stats.out"
STATS="${STATS%%.*}.stats"
STDOUT=$OUT/Bowtie_rRNA_Removed/STDOUT/$NAME".out"
STDOUT="${STDOUT%%.*}.stdout"

# Removes rRNA and tRNA (could consider another program)

$BOWTIE_PATH/bowtie --threads 4 \
-l 20 \
--al ${ALIGNED} \
--un ${UNALIGNED} \
${BOWTIE_REF} \
${NO_UMI} \
1>>${STDOUT} \
2>>${STATS}

echo "Bowtie alignment is completed"
# rm $ALIGNED $STDOUT $TRIMMED

#STAR align

mkdir -p $OUT/STAR/

STDERR=$OUT/STAR/stderror.stderr
STDO=$OUT/STAR/stdout.stdout
OUTPUT=$OUT/STAR/file_

# Taking unaligned input from Bowtie (filtered for tRNA) as input

$STAR_PATH/STAR \
--runThreadN 4 \
--genomeDir $star_genome \
--readFilesIn $UNALIGNED --alignIntronMin 20 --alignIntronMax 100000 \
--outFilterMismatchNmax 1 --outFilterType BySJout \
--outFilterMismatchNoverLmax 0.04 \
--twopassMode Basic \
--outFileNamePrefix $OUTPUT \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
1>>${STDO} \
2>>${STDERR}

echo "STAR alignment is done"

STAR_BAM=$OUT/STAR/file_Aligned.sortedByCoord.out.bam
samtools index $STAR_BAM

# Run FastQC
mkdir -p $OUT/FastQC/
$FASTQC_PATH/fastqc -o $OUT/FastQC -f bam_mapped $STAR_BAM

mkdir -p $OUT/Fastq_FastQC/
$FASTQC_PATH/fastqc -o $OUT/Fastq_FastQC -f fastq $UNALIGNED

# RibORF plot trinucleotide periodicity

mkdir -p $OUT/RibORF/Sam/

SAM=$OUT/RibORF/Sam/$(basename $UNALIGNED)".sam"
SAM="${SAM%%.*}.sam"
HEADER=$OUT/RibORF/Sam/$(basename $UNALIGNED)".header"
HEADER="${HEADER%%.*}.header"

samtools view $STAR_BAM -o $SAM
samtools view -H $STAR_BAM -o $HEADER

for i in 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
	do
		mkdir -p $OUT/RibORF/Trinuc_plots/$i

		OUTPUT=$OUT/RibORF/Trinuc_plots/$i

		cd $OUTPUT

		perl $PERIODICTY_TOOL_PATH/readDist_modified2.pl \
		${SAM} \
		${GENE_PRED} \
		$OUTPUT $i 30 50 \

	done

perl $PERIODICTY_TOOL_PATH/offsetCorrect_samflag_corrected.pl \
$SAM_FILE \
$OFFSET \
$OUT

# Convert the resulting SAM file to BAM again. Need a header for this.
# THIS IS THE OUTPUT FOR QUANTITATION WITH STAR AND RSEM
BAM=$OUT'/'$NAME'.offset.bam'

cat $SAM_HEADER $OUT'/corrected.sam' | samtools view -b -o $BAM

# Delete the resulting SAM file once certain that BAM file has been successfully generated:
#rm $SAM_FILE
#rm $OUT'/header.sam'
rm $OUT'/corrected.sam'

all_reads=$(samtools view -c $BAM_FILE)
echo "all reads count "$all_reads
inframe_reads=$(samtools view -c $BAM)
echo "in frame reads count "$inframe_reads
