#!/bin/bash -l
#########################################################
### needed time still to be estimated. Run with different commands
### so if interrupt, start againg from last file.
### now set to 2:00:00 but probably need much less
### update: job killed after 2:00:00. 5:00:00 is enough for
### running BWA
### update2: split over different process. Adapt timing accordingly
### depending on task
#########################################################
#PBS -l walltime=23:00:00
#PBS -L tasks=1:lprocs=20

module load BWA/0.7.17-GCCcore-8.3.0
module load SAMtools/1.9-intel-2019b
module load Java

### STEP 0: set somy variables
export DATA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/
export BIN_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/
export BWA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/
export PICARD_JAR=/user/antwerpen/205/vsc20587/software/picard/picard.jar
export PICARD_BIN="java -jar ${PICARD_JAR}"
export MIN_L=100
export THREADS=20
export SEED=100

### STEP 1: filter out short reads ###
# reads_filtering_on_length.sh: this script is based on the idea 
# suggested here: http://seqanswers.com/forums/showthread.php?t=31845 
# and write all line number of *paired* data and next filters them
# adpated so that lines are written to unqiue file instaed of a 
# generic file
# $BIN_DIR/reads_filtering_on_length.sh $DATA_DIR/1_S1_L001_R1_001.fastq $DATA_DIR/1_S1_L001_R2_001.fastq $MIN_L
# $BIN_DIR/reads_filtering_on_length.sh $DATA_DIR/2_S2_L001_R1_001.fastq $DATA_DIR/2_S2_L001_R2_001.fastq $MIN_L
# $BIN_DIR/reads_filtering_on_length.sh $DATA_DIR/3_S3_L001_R1_001.fastq $DATA_DIR/3_S3_L001_R2_001.fastq $MIN_L
# $BIN_DIR/reads_filtering_on_length.sh $DATA_DIR/4_S4_L001_R1_001.fastq $DATA_DIR/4_S4_L001_R2_001.fastq $MIN_L

### STEP 2: rename and gzip ###
## rename all fastq-files and gzip them 
# mv $DATA_DIR/1_S1_L001_R1_001.fastq.100 $DATA_DIR/1_S1_L001_R1_001.filter.fastq 
# mv $DATA_DIR/1_S1_L001_R2_001.fastq.100 $DATA_DIR/1_S1_L001_R2_001.filter.fastq 
# mv $DATA_DIR/2_S2_L001_R1_001.fastq.100 $DATA_DIR/2_S2_L001_R1_001.filter.fastq 
# mv $DATA_DIR/2_S2_L001_R2_001.fastq.100 $DATA_DIR/2_S2_L001_R2_001.filter.fastq 
# mv $DATA_DIR/3_S3_L001_R1_001.fastq.100 $DATA_DIR/3_S3_L001_R1_001.filter.fastq 
# mv $DATA_DIR/3_S3_L001_R2_001.fastq.100 $DATA_DIR/3_S3_L001_R2_001.filter.fastq 
# mv $DATA_DIR/4_S4_L001_R1_001.fastq.100 $DATA_DIR/4_S4_L001_R1_001.filter.fastq 
# mv $DATA_DIR/4_S4_L001_R2_001.fastq.100 $DATA_DIR/4_S4_L001_R2_001.filter.fastq 

## if needed, first do gzip of the fastq files
# gzip $DATA_DIR/1_S1_L001_R1_001.filter.fastq 
# gzip $DATA_DIR/1_S1_L001_R2_001.filter.fastq 
# gzip $DATA_DIR/2_S2_L001_R1_001.filter.fastq 
# gzip $DATA_DIR/2_S2_L001_R2_001.filter.fastq 
# gzip $DATA_DIR/3_S3_L001_R1_001.filter.fastq 
# gzip $DATA_DIR/3_S3_L001_R2_001.filter.fastq 
# gzip $DATA_DIR/4_S4_L001_R1_001.filter.fastq 
# gzip $DATA_DIR/4_S4_L001_R2_001.filter.fastq 


#### STEP 3: run BWA #### 

# run indexing on the genome
# bwa index /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta

# # run BWA on different fastq files
bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/1_S1_L001_R1_001.filter.fastq.gz $DATA_DIR/1_S1_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/1_S1.bam -

bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/2_S2_L001_R1_001.filter.fastq.gz $DATA_DIR/2_S2_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/2_S2.bam -

bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/3_S3_L001_R1_001.filter.fastq.gz $DATA_DIR/3_S3_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/3_S3.bam -

bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/4_S4_L001_R1_001.filter.fastq.gz $DATA_DIR/4_S4_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/4_S4.bam -

# index all the bam files
samtools index $BWA_DIR/1_S1.bam
samtools index $BWA_DIR/2_S2.bam
samtools index $BWA_DIR/3_S3.bam
samtools index $BWA_DIR/4_S4.bam

#### STEP 4: remove duplicates ####
$PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/1_S1.bam O=$BWA_DIR/1_S1.removedups.bam M=$BWA_DIR/1_S1.markdups_metrics.txt
$PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/2_S2.bam O=$BWA_DIR/2_S2.removedups.bam M=$BWA_DIR/2_S2.markdups_metrics.txt
$PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/3_S3.bam O=$BWA_DIR/3_S3.removedups.bam M=$BWA_DIR/3_S3.markdups_metrics.txt
$PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/4_S4.bam O=$BWA_DIR/4_S4.removedups.bam M=$BWA_DIR/4_S4.markdups_metrics.txt

# index all the bam files
samtools index $BWA_DIR/1_S1.removedups.bam
samtools index $BWA_DIR/2_S2.removedups.bam
samtools index $BWA_DIR/3_S3.removedups.bam
samtools index $BWA_DIR/4_S4.removedups.bam

#### STEP 5: select only the proper paired reads
samtools view -@$THREADS -bf 0x2 $BWA_DIR/1_S1.removedups.bam > $BWA_DIR/1_S1.removedups.proper_paired.bam 
samtools view -@$THREADS -bf 0x2 $BWA_DIR/2_S2.removedups.bam > $BWA_DIR/2_S2.removedups.proper_paired.bam 
samtools view -@$THREADS -bf 0x2 $BWA_DIR/3_S3.removedups.bam > $BWA_DIR/3_S3.removedups.proper_paired.bam
samtools view -@$THREADS -bf 0x2 $BWA_DIR/4_S4.removedups.bam > $BWA_DIR/4_S4.removedups.proper_paired.bam 

# index all the bam files
samtools index $BWA_DIR/1_S1.removedups.proper_paired.bam
samtools index $BWA_DIR/2_S2.removedups.proper_paired.bam
samtools index $BWA_DIR/3_S3.removedups.proper_paired.bam
samtools index $BWA_DIR/4_S4.removedups.proper_paired.bam

#### STEP 6: do subsampling towards lowest coverage
# this step was previously performed in a separate script (subsampling.sh), but 
# now integrated in main script. 

## 6.1 run flagstat
samtools flagstat $BWA_DIR/1_S1.removedups.proper_paired.bam > $BWA_DIR/1_S1.removedups.proper_paired.flagstat
samtools flagstat $BWA_DIR/2_S2.removedups.proper_paired.bam > $BWA_DIR/2_S2.removedups.proper_paired.flagstat
samtools flagstat $BWA_DIR/3_S3.removedups.proper_paired.bam > $BWA_DIR/3_S3.removedups.proper_paired.flagstat
samtools flagstat $BWA_DIR/4_S4.removedups.proper_paired.bam > $BWA_DIR/4_S4.removedups.proper_paired.flagstat

## 6.2 run subsampling
# first do calculations before you can do subsampling.





#### OLD COMMANDS ###

# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/3_S3_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/3_S3_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/3_S3.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/1_S1_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/1_S1_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/1_S1.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/4_S4_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/4_S4_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/4_S4.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/2_S2_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/2_S2_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/2_S2.bam -
