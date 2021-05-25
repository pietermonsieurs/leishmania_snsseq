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
#PBS -l walltime=0:59:00
#PBS -L tasks=1:lprocs=20

module load BWA/0.7.17-GCCcore-8.3.0
module load SAMtools/1.9-intel-2019b
module load BEDTools/2.27.1-intel-2018b
module load Java
module load Python/3

### STEP 0: set somy variables
export DATA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/
export BIN_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/
export BWA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/
export SICER_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/sicer/
export MACS_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/macs/
export OVERLAP_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/overlap/
export PICARD_JAR=/user/antwerpen/205/vsc20587/software/picard/picard.jar
export PICARD_BIN="java -jar ${PICARD_JAR}"
export MIN_L=100
export THREADS=20
export SEED=100
export mapq_cutoff=30


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
# bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/1_S1_L001_R1_001.filter.fastq.gz $DATA_DIR/1_S1_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/1_S1.bam -

# bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/2_S2_L001_R1_001.filter.fastq.gz $DATA_DIR/2_S2_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/2_S2.bam -

# bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/3_S3_L001_R1_001.filter.fastq.gz $DATA_DIR/3_S3_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/3_S3.bam -

# bwa mem -t $THREADS -k $SEED $DATA_DIR/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta $DATA_DIR/4_S4_L001_R1_001.filter.fastq.gz $DATA_DIR/4_S4_L001_R2_001.filter.fastq.gz | samtools sort -@$THREADS -o $BWA_DIR/4_S4.bam -

# index all the bam files
# samtools index $BWA_DIR/1_S1.bam
# samtools index $BWA_DIR/2_S2.bam
# samtools index $BWA_DIR/3_S3.bam
# samtools index $BWA_DIR/4_S4.bam

#### STEP 4: filtering on Mapping quality ####
# samtools view -@ $THREADS -bq $mapq_cutoff $BWA_DIR/1_S1.bam > $BWA_DIR/1_S1.mapq${mapq_cutoff}.bam
# samtools view -@ $THREADS -bq $mapq_cutoff $BWA_DIR/2_S2.bam > $BWA_DIR/2_S2.mapq${mapq_cutoff}.bam
# samtools view -@ $THREADS -bq $mapq_cutoff $BWA_DIR/3_S3.bam > $BWA_DIR/3_S3.mapq${mapq_cutoff}.bam
# samtools view -@ $THREADS -bq $mapq_cutoff $BWA_DIR/4_S4.bam > $BWA_DIR/4_S4.mapq${mapq_cutoff}.bam




#### STEP 5: remove duplicates ####
# $PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/1_S1.mapq${mapq_cutoff}.bam O=$BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.bam M=$BWA_DIR/1_S1.markdups_metrics.txt
# $PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/2_S2.mapq${mapq_cutoff}.bam O=$BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.bam M=$BWA_DIR/2_S2.markdups_metrics.txt
# $PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/3_S3.mapq${mapq_cutoff}.bam O=$BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.bam M=$BWA_DIR/3_S3.markdups_metrics.txt
# $PICARD_BIN MarkDuplicates REMOVE_DUPLICATES=true I=$BWA_DIR/4_S4.mapq${mapq_cutoff}.bam O=$BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.bam M=$BWA_DIR/4_S4.markdups_metrics.txt

# index all the bam files
# samtools index $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.bam
# samtools index $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.bam
# samtools index $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.bam
# samtools index $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.bam

#### STEP 6: select only the proper paired reads
# samtools view -@$THREADS -bf 0x2 $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.bam > $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.bam 
# samtools view -@$THREADS -bf 0x2 $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.bam > $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.bam 
# samtools view -@$THREADS -bf 0x2 $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.bam > $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.bam
# samtools view -@$THREADS -bf 0x2 $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.bam > $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.bam 

# # index all the bam files
# samtools index $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.bam
# samtools index $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.bam
# samtools index $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.bam
# samtools index $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.bam

#### STEP 7: do subsampling towards lowest coverage
# this step was previously performed in a separate script (subsampling.sh), but 
# now integrated in main script. 

## 7.1 run flagstat
# samtools flagstat $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.bam > $BWA_DIR/1_S1.removedups.proper_paired.flagstat
# samtools flagstat $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.bam > $BWA_DIR/2_S2.removedups.proper_paired.flagstat
# samtools flagstat $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.bam > $BWA_DIR/3_S3.removedups.proper_paired.flagstat
# samtools flagstat $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.bam > $BWA_DIR/4_S4.removedups.proper_paired.flagstat

## 7.2 run subsampling
# first do calculations before you can do subsampling. Those fractions should
# be taken from the excel file. Should be optimized in a later stage, where all
# # bam files are normalized. 

# # sample S1
# ln -s $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.bam $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

# # sample S2
# samtools view -@$THREADS -bs 10.4927 $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.bam > $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

# # sample S3
# samtools view -@$THREADS -bs 10.7776 $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.bam > $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

# # sample S4
# samtools view -@$THREADS -bs 10.7493 $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.bam > $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

# samtools index $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

# samtools index $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

# samtools index $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

# samtools index $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam

#### STEP 8: peak detection


# python path should be adapted to the local python library so that the MACS and SICER package can be found
# export PYTHONPATH=/user/antwerpen/205/vsc20587/software/python_lib/lib/python3.7/site-packages/:$PYTHONPATH
export PYTHONPATH=/data/antwerpen/205/vsc20587/software_from_home/python_lib/lib/python3.7/site-packages/:$PYTHONPATH


# 8.1 run MACS2 

# sample S1 - control S2
/user/antwerpen/205/vsc20587/software/python_lib/bin/macs2 callpeak -f BAM -g 38e6 --nomodel --extsize 150 -n S1 -t $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam -c $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam --outdir $MACS_DIR

# sample S3 - control S4
/user/antwerpen/205/vsc20587/software/python_lib/bin/macs2 callpeak -f BAM -g 38e6 --nomodel --extsize 150  -n S3 -t $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam -c $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam --outdir $MACS_DIR

# sample S1 - control S2 - BROAD peak
/user/antwerpen/205/vsc20587/software/python_lib/bin/macs2 callpeak -f BAM -g 38e6 --nomodel --extsize 150 --broad -n S1 -t $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam -c $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam --outdir $MACS_DIR

# sample S3 - control S4 - BROAD peak
/user/antwerpen/205/vsc20587/software/python_lib/bin/macs2 callpeak -f BAM -g 38e6 --nomodel --extsize 150 --broad -n S3 -t $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam -c $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam --outdir $MACS_DIR


# 8.2 convert bam files to bed (bed-files = input for Sicer)
# converting the bam files to bed is needed to run SICER. The code below is more 
# or less a copy paste from the bash script "bedtools_bamtobed.sh"
bedtools bamtobed -i $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam > $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed

bedtools bamtobed -i $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam > $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed

bedtools bamtobed -i $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam > $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed

bedtools bamtobed -i $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bam > $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed


# 8.3 run SICER 2

/user/antwerpen/205/vsc20587/software/python_lib/bin/sicer -t $BWA_DIR/1_S1.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed -c $BWA_DIR/2_S2.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed -o $SICER_DIR -s tbruc --cpu $THREADS

/user/antwerpen/205/vsc20587/software/python_lib/bin/sicer -t $BWA_DIR/3_S3.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed -c $BWA_DIR/4_S4.mapq${mapq_cutoff}.removedups.proper_paired.subsample.bed -o $SICER_DIR -s tbruc --cpu $THREADS


#### STEP 9: overlap between peaks

# step 9.1: convert output of both algorithms to three-column bed-files.

# only absolutely needed for SICER, not needed for MACS

# sicer peaks
awk '{print $1"\t"$2"\t"$3}' $SICER_DIR/1_S1.mapq30.removedups.proper_paired.subsample-W200-G600-FDR0.01-island.bed > $SICER_DIR/S1.bed

awk '{print $1"\t"$2"\t"$3}' $SICER_DIR/3_S3.mapq30.removedups.proper_paired.subsample-W200-G600-FDR0.01-island.bed > $SICER_DIR/S3.bed

# macs narrow peaks - not needed / you can use the default file 
awk '{print $1"\t"$2"\t"$3}' $MACS_DIR/S1_peaks.narrowPeak > $MACS_DIR/S1_peaks.narrowPeak.bed
awk '{print $1"\t"$2"\t"$3}' $MACS_DIR/S3_peaks.narrowPeak > $MACS_DIR/S3_peaks.narrowPeak.bed

# macs broad peaks - not needed / you can use the default file 
awk '{print $1"\t"$2"\t"$3}' $MACS_DIR/S1_peaks.broadPeak > $MACS_DIR/S1_peaks.broadPeak.bed
awk '{print $1"\t"$2"\t"$3}' $MACS_DIR/S3_peaks.broadPeak > $MACS_DIR/S3_peaks.broadPeak.bed


# step 9.2: look for overlap of MACS with SICER peaks 

# S1 - setting -f and -F has no influence with current files
bedtools intersect -wa -f 1e-6 -F 1e-6 -a $MACS_DIR/S1_peaks.narrowPeak -b $SICER_DIR/S1.bed > $OVERLAP_DIR/S1.narrowPeak.bed
bedtools intersect -wa -f 1e-6 -F 1e-6 -a $MACS_DIR/S1_peaks.broadPeak -b $SICER_DIR/S1.bed > $OVERLAP_DIR/S1.broadPeak.bed

# S3
bedtools intersect -wa -f 1e-6 -F 1e-6 -a $MACS_DIR/S3_peaks.narrowPeak -b $SICER_DIR/S3.bed > $OVERLAP_DIR/S3.narrowPeak.bed
bedtools intersect -wa -f 1e-6 -F 1e-6 -a $MACS_DIR/S3_peaks.broadPeak -b $SICER_DIR/S3.bed > $OVERLAP_DIR/S3.broadPeak.bed



# step 9.3: select for unique peaks in both samples
bedtools intersect -v -f 1e-6 -F 1e-6 -a $OVERLAP_DIR/S1.narrowPeak.bed -b $OVERLAP_DIR/S3.narrowPeak.bed > $OVERLAP_DIR/S1.narrowPeak.unique.bed

bedtools intersect -v -f 1e-6 -F 1e-6 -b $OVERLAP_DIR/S1.narrowPeak.bed -a $OVERLAP_DIR/S3.narrowPeak.bed > $OVERLAP_DIR/S3.narrowPeak.unique.bed

bedtools intersect -v -f 1e-6 -F 1e-6 -a $OVERLAP_DIR/S1.broadPeak.bed -b $OVERLAP_DIR/S3.broadPeak.bed > $OVERLAP_DIR/S1.broadPeak.unique.bed

bedtools intersect -v -f 1e-6 -F 1e-6 -b $OVERLAP_DIR/S1.broadPeak.bed -a $OVERLAP_DIR/S3.broadPeak.bed > $OVERLAP_DIR/S3.broadPeak.unique.bed




#### OLD COMMANDS ###

# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/3_S3_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/3_S3_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/3_S3.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/1_S1_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/1_S1_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/1_S1.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/4_S4_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/4_S4_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/4_S4.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/2_S2_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/2_S2_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/2_S2.bam -
