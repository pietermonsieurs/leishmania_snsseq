#!/bin/bash -l

#PBS -l walltime=0:59:00
#PBS -L tasks=6:lprocs=4


# fg=19_S70
# bg=20_S71
# bg=gDNA_9_S9

module load Python/3.7.4-intel-2019b

## parameter settings
export PYTHONPATH=/data/antwerpen/205/vsc20587/software_from_home/python_lib/lib/python3.7/site-packages/:$PYTHONPATH
BWA_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/
OUT_DIR=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/macs_newdata/

/data/antwerpen/205/vsc20587/software_from_home/python_lib/bin/macs2 callpeak \
     -f BAM \
     -g 38e6 \
     --nomodel --extsize 150 \
     -n ${fg}_vs_${bg} \
     -t $BWA_DIR/${fg}.mapq30.removedups.proper_paired.subsample.bam \
     -c $BWA_DIR/${bg}.mapq30.removedups.proper_paired.subsample.bam \
     --outdir $OUT_DIR \
     --scale-to small


## module worker
# module load worker
# wsub -data /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/parameters_macs.csv -batch /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/macs_newdata.sh

