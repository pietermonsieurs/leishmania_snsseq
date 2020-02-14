#!/bin/bash -l

# do subsampling using samtools. First do the 
# flagstat option 

#PBS -l walltime=2:00:00
#PBS -L tasks=1:lprocs=20


# load necessary modules
module load SAMtools/1.9-intel-2019b

# set necessary variables
export BAM_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa/


# running flagstats to know the amount of reads
for sample in 1_S1 2_S2 3_S3 4_S4
do
    echo $BAM_DIR/$sample.bam
    # samtools flagstat $BAM_DIR/$sample.bam > $BAM_DIR/$sample.flagstat
done

# calculate relative percentages versus the smallest bam file, which is in 
# this case 4_S4_bam
# | Sample | All_reads | Percentage |
# |--------|-----------|------------|
# | S1     | 65364709  | 68%        |
# | S2     | 91426202  | 49%        |
# | S3     | 250365416 | 18%        |
# | S4     | 44428519  | 100%       |

# percentage is .68. Seed is 10. ==> s 10.68
# samtools view -bs 10.68 $BAM_DIR/1_S1.bam > $BAM_DIR/1_S1.subsample.bam
# samtools view -bs 10.49 $BAM_DIR/2_S2.bam > $BAM_DIR/2_S2.subsample.bam
samtools view -bs 10.18 $BAM_DIR/3_S3.bam > $BAM_DIR/3_S3.subsample.bam
# ln -s $BAM_DIR/4_S4.bam $BAM_DIR/4_S4.subsample.bam


