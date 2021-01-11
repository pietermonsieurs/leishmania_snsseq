#!/bin/bash -l

#PBS -l walltime=08:00:00
#PBS -L tasks=1:lprocs=28

# load modules
module load BWA/0.7.17-GCCcore-8.3.0
module load SAMtools/1.9-intel-2019b
module load Java

# define variables
export THREADS=28
export SEED=100
export mapq_cutoff=30
export BWA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/
export REF_GENOME=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta
export PICARD_JAR=/user/antwerpen/205/vsc20587/data/software/picard/picard.jar
export PICARD_BIN="java -jar ${PICARD_JAR}"

### sequencing data test
# fastq_file_1=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/20210104/10_S61_L002_R1_001.fastq.gz
echo $fastq_file_1



## first create the prefix to be used to create fastq_file R2
file_prefix_full=${fastq_file_1%_L002_R1_001.fastq.gz}
fastq_file_2=${file_prefix_full}_L002_R2_001.fastq.gz

## create bwa bam output file name
file_prefix=${file_prefix_full##*/}
bam_file_prefix=${BWA_DIR}/${file_prefix}

### index reference genomes
# bwa index $REF_GENOME


# run BWA aligner
bwa mem -k $SEED -t $THREADS $REF_GENOME ${fastq_file_1} ${fastq_file_2} | samtools sort -@$THREADS -o ${bam_file_prefix}.bam

### selection on mapping quality ###
samtools view -@ $THREADS -bq $mapq_cutoff ${bam_file_prefix}.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.bam
samtools index ${bam_file_prefix}.mapq${mapq_cutoff}.bam
samtools flagstat ${bam_file_prefix}.mapq${mapq_cutoff}.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.flagstat

### remove duplicate reads ###
$PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=${bam_file_prefix}.mapq${mapq_cutoff}.bam O=${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam M=${bam_file_prefix}.mapq${mapq_cutoff}.markdups_metrics.txt
samtools index ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam
samtools flagstat ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.flagstat


### select only proper paired reads ###
samtools view -@$THREADS -bf 0x2 ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.bam
samtools index ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.bam
samtools flagstat ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.bam > ${bam_file_prefix}.mapq${mapq_cutoff}.removedups.proper_paired.flagstat





## run for all files
## instead of using the worker module, all sampels are run by specifying the fastq-file1 
## in the command line when calling qsub. Those commandline files can be created using
## the Pyhton script: bwa_snsseq_newdata_create_qsub_commands.py
