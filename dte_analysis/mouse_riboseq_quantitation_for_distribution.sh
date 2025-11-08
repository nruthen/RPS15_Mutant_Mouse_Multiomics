#! /bin/bash
#$ -cwd
#$ -l h_vmem=16g
#$ -l h_rt=24:00:00
#$ -pe smp 4 -R y -binding linear:4

source /broad/software/scripts/useuse

reuse .star-2.7.0a
reuse Rsem
reuse Samtools

SAMPLESHEET=$1 # sample sheet with one sample name per line
OUT=$2 # where output files will be created
STAR_GENOME=$3 # /cga/wu/riboseq/reference/Mouse/STAR_Riboseq_Mouse_Corrected
STAR_TRANSCRIPTOME=$4 # /cga/wu/riboseq/reference/Mouse/STAR_Riboseq_Mouse_Transcriptome
RSEM_REFERENCE=$5 #/cga/wu/riboseq/reference/Mouse/RSEM_Riboseq_Mouse_Corrected/RSEM_Riboseq_Mouse_Corrected

SAMPLE=$(awk "NR==${SGE_TASK_ID}" ${SAMPLESHEET})

SAMPLE_ID=$(basename $SAMPLE)
SAMPLE_ID=${SAMPLE_ID%%.*}

OUTDIR=$OUT/$SAMPLE_ID

#mkdir -p $OUTDIR/STAR_1
#mkdir -p $OUTDIR/STAR_2

OUTPUT_1=$OUTDIR"/STAR_1/"$SAMPLE_ID"_"

OUTPUT_2=$OUTDIR"/STAR_2/"$SAMPLE_ID"_"

# Align reads agnostically to genome to filter out reads that should not align to the transcriptome 
STAR \
--runThreadN 4 \
--genomeDir $STAR_GENOME \
--readFilesIn $SAMPLE --alignIntronMin 20 --alignIntronMax 100000 \
--outFilterMismatchNmax 1 --outFilterType BySJout \
--outFilterMismatchNoverLmax 0.04 \
--twopassMode Basic \
--outFileNamePrefix $OUTPUT_1 \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000

# Extract the reads aligning to the transcriptome and convert into fastq format for re-alignment
samtools bam2fq $OUTPUT_1'Aligned.toTranscriptome.out.bam' > $OUTPUT_1'genome_aligned.fastq'

# Extract reads that align to the pseudo-transcriptome (padded canonical orfeome)
STAR \
--runThreadN 4 \
--genomeDir $STAR_TRANSCRIPTOME \
--readFilesIn $OUTPUT_1'genome_aligned.fastq' --alignEndsType EndToEnd --alignIntronMax 1 --alignIntronMin 2 --scoreDelOpen -10000 --scoreInsOpen -10000 \
--outFilterMismatchNmax 1 --outFilterType BySJout \
--outFilterMismatchNoverLmax 0.04 \
--twopassMode Basic \
--outFileNamePrefix $OUTPUT_2 \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000
mkdir -p $OUTDIR/RSEM
cd $OUTDIR/RSEM

# Quantify expression values for the canonical orfeome
rsem-calculate-expression \
--bam \
$OUTPUT_2'Aligned.sortedByCoord.out.bam' \
$RSEM_REFERENCE \
$SAMPLE_ID



#qsub -t 1-10 STAR_mouse_RNA_align_quant_array.sh star_sheet.txt ~/out ~/STAR_Riboseq_Mouse_Corrected ~/STAR_Riboseq_Mouse_Transcriptome ~/RSEM_Riboseq_Mouse_Corrected
