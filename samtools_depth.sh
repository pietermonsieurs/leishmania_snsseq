#!/bin/bash -l

# takes only a very short time. For four 2.2Gb bam files,
# with many contigs (> 200), only a few minutes. q1h should
# be sufficient

#PBS -l walltime=0:20:00
#PBS -L tasks=20:lprocs=1

# load necessary modules
module load SAMtools/1.9-intel-2019b

# separate file made with all chromosome names and all 4 samples
# in order to get the depth for the complete contigs. So the worker
# module will loop over 2 variables (sample and chrom)

export BAM_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/
export DEPTH_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/depth/

cd $DEPTH_DIR/

### example data 
# export sample=1
# export chrom=Chr11_3B_Tb427v10

# samtools depth -a -r $chrom $BAM_DIR/${sample}_S${sample}.proper_paired.subsample.bam > $DEPTH_DIR/${sample}_S${sample}.chrom_${chrom}.subsample.depth.csv

samtools index $chrom $BAM_DIR/${sample}_S${sample}.mapq30.removedups.proper_paired.subsample.bam

samtools depth -a -r $chrom $BAM_DIR/${sample}_S${sample}.mapq30.removedups.proper_paired.subsample.bam > $DEPTH_DIR/${sample}_S${sample}.chrom_${chrom}.subsample.depth.csv



# wsub commands
# create csv file in excel with all samples and all contigs
# /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/depth/input_sns_samtools_depth.csv
# module load worker
# wsub -data /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/depth/input_sns_samtools_depth.csv -batch /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/samtools_depth.sh

