#!/usr/bin/env python3


import os

my_debug = 1
src_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/'
data_dir = src_dir + 'data/20210104/'
bin_dir = src_dir + 'bin/'
bwa_bin = bin_dir + 'bwa_snsseq_newdata.sh'


for fq1 in os.listdir(data_dir):

    # only select the first reads, if not a file containing the 
    # first read, skip to next one
    if not fq1.endswith("R1_001.fastq.gz"):
        continue

    qsub_command = f"qsub -v fastq_file_1={data_dir}{fq1} {bwa_bin}"

    # fq2 = fq1.replace("R1.fastq.gz", "R2.fastq.gz")

    # fq1 = data_dir + fq1
    # fq2 = data_dir + fq2

    # qsub_command = f"qsub -v fastq_file_1={fq1},fastq_file_3={fq2} {bwa_bin}"
    print(qsub_command)
    
