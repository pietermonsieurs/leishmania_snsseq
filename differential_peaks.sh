#!/bin/bash -l
#########################################################
### needed time still to be estimated. Run with different commands
#########################################################
#PBS -l walltime=0:20:00
#PBS -L tasks=20:lprocs=1

module load Python/3

export plotting_program=/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/differential_peaks.py



$plotting_program -c $chrom -s $start -e $end -d $subdir 


###### create input information 

### S1 - narrow peaks ###
# echo 'chrom,start,end,subdir' > S1.narrowPeak.unique.input_plotting.csv
# awk '{print $1","$2","$3",S1_narrowPeak"}' S1.narrowPeak.unique.bed >>  S1.narrowPeak.unique.input_plotting.csv
# module load worker
# wsub -data /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/overlap/S1.narrowPeak.unique.input_plotting.csv -batch /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/differential_peaks.sh

### S1 - broad peaks ###
# echo 'chrom,start,end,subdir' > S1.broadPeak.unique.input_plotting.csv
# awk '{print $1","$2","$3",S1_broadPeak"}' S1.broadPeak.unique.bed >>  S1.broadPeak.unique.input_plotting.csv
# module load worker
# wsub -data /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/overlap/S1.broadPeak.unique.input_plotting.csv -batch /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/differential_peaks.sh


### S3 - narrow peaks ###
# echo 'chrom,start,end,subdir' > S3.narrowPeak.unique.input_plotting.csv
# awk '{print $1","$2","$3",S3_narrowPeak"}' S3.narrowPeak.unique.bed >>  S3.narrowPeak.unique.input_plotting.csv
# module load worker
# wsub -data /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/overlap/S3.narrowPeak.unique.input_plotting.csv -batch /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/differential_peaks.sh

### S3 - broad peaks ###
# echo 'chrom,start,end,subdir' > S3.broadPeak.unique.input_plotting.csv
# awk '{print $1","$2","$3",S3_broadPeak"}' S3.broadPeak.unique.bed >>  S3.broadPeak.unique.input_plotting.csv
# module load worker
# wsub -data /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/overlap/S3.broadPeak.unique.input_plotting.csv -batch /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/differential_peaks.sh


### Centromers ###
# module load worker
# wsub -data /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/centromers/centromers.csv -batch /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/bin/differential_peaks.sh





