#!/bin/bash -l
#########################################################
### needed time still to be estimated. Run with different commands
### so if interrupt, start againg from last file.
### now set to 2:00:00 but probably need much less
### update: job killed after 2:00:00
#########################################################
#PBS -l walltime=3:00:00
#PBS -L tasks=1:lprocs=20

module load BWA/0.7.17-GCCcore-8.3.0
module load SAMtools/1.9-intel-2019b

# run indexing on the genome
# bwa index /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta

# run BWA on different fastq files
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/3_S3_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/3_S3_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/3_S3.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/1_S1_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/1_S1_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/1_S1.bam -
# bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/4_S4_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/4_S4_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/4_S4.bam -
bwa mem -t 20 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/2_S2_L001_R1_001.fastq.gz /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/2_S2_L001_R2_001.fastq.gz | samtools sort -@20 -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/2_S2.bam -