#!/bin/bash -l
#########################################################
### needed time still to be estimated. Run with different commands
### so if interrupt, start againg from last file.
### now set to 2:00:00 but probably need much less
### update: job killed after 2:00:00
#########################################################
#PBS -l walltime=0:59:00
#PBS -L tasks=1:lprocs=28

# load necessary modules
module load SAMtools/1.9-intel-2019b
module load BEDTools/2.27.1-intel-2018b

# set necessary variables
export BAM_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/

cd $BAM_DIR

bedtools bamtobed -i 1_S1.mapq30.removedups.proper_paired.subsample.bam > 1_S1.mapq30.removedups.proper_paired.subsample.bed

bedtools bamtobed -i 2_S2.mapq30.removedups.proper_paired.subsample.bam > 2_S2.mapq30.removedups.proper_paired.subsample.bed

bedtools bamtobed -i 3_S3.mapq30.removedups.proper_paired.subsample.bam > 3_S3.mapq30.removedups.proper_paired.subsample.bed

bedtools bamtobed -i 4_S4.mapq30.removedups.proper_paired.subsample.bam > 4_S4.mapq30.removedups.proper_paired.subsample.bed
