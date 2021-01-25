#!/bin/bash -l

# do subsampling using samtools. First do the 
# flagstat option 

#PBS -l walltime=0:59:00
#PBS -L tasks=1:lprocs=28

module load SAMtools


ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/2_S53.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/2_S53.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.72 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/1_S52.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/1_S52.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/4_S55.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/4_S55.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.87 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/3_S54.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/3_S54.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/6_S57.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/6_S57.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.69 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/5_S56.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/5_S56.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/8_S59.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/8_S59.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.87 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/7_S58.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/7_S58.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/9_S60.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/9_S60.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.66 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/10_S61.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/10_S61.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/11_S62.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/11_S62.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.73 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/12_S63.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/12_S63.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/14_S65.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/14_S65.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.45 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/13_S64.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/13_S64.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/16_S67.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/16_S67.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.28 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/15_S66.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/15_S66.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/18_S69.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/18_S69.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.82 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/17_S68.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/17_S68.mapq30.removedups.proper_paired.subsample.bam
ln -s /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/20_S71.mapq30.removedups.proper_paired.bam /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/20_S71.mapq30.removedups.proper_paired.subsample.bam
samtools view --threads 26 -bs 10.93 /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/19_S70.mapq30.removedups.proper_paired.bam > /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/19_S70.mapq30.removedups.proper_paired.subsample.bam
