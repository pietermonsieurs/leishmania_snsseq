#!/usr/bin/env python3

import os
import re


data_dir =  '/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/bwa_newdata/'
threads = 26
seed = 10
my_debug = 0

flagstat_files = [f for f in os.listdir(data_dir) if f.endswith('proper_paired.flagstat')]
my_debug and print(flagstat_files)

def find_number_of_reads (flagstat_file):
    flagstat_fh = open(flagstat_file, 'r')
    for line in flagstat_fh:
        if "properly paired" in line:
            nr_of_reads = line.split()[0]
            my_debug and print(f"flagstat_file {flagstat_file} --> {nr_of_reads}")
            return(int(nr_of_reads))

for sample_id in range(1,20,2):
    ## the first ID is always the sample, the second id is the control ID which is paired
    ## with the sample_id
    control_id = sample_id + 1
    
    file_sample = [file_name for file_name in flagstat_files if file_name.startswith(f"{sample_id}_")][0]
    file_control = [file_name for file_name in flagstat_files if file_name.startswith(f"{control_id}_")][0]

    file_sample = data_dir + file_sample
    file_control = data_dir + file_control

    my_debug and print(f"sample file {file_sample} -- control file {file_control}")
    nr_of_reads_sample = find_number_of_reads(file_sample)
    nr_of_reads_control = find_number_of_reads(file_control)

    if nr_of_reads_sample > nr_of_reads_control:
        bam_file_keep = file_control.replace(".flagstat", ".bam")
        bam_file_keep_ss = file_control.replace(".flagstat", ".subsample.bam")

        bam_file_reduce = file_sample.replace(".flagstat", ".bam")
        bam_file_reduce_ss = file_sample.replace(".flagstat", ".subsample.bam")

        ratio = nr_of_reads_control/nr_of_reads_sample
        ratio = round(100*ratio)
    elif nr_of_reads_sample < nr_of_reads_control:
        bam_file_keep = file_sample.replace(".flagstat", ".bam")
        bam_file_keep_ss = file_sample.replace(".flagstat", ".subsample.bam")

        bam_file_reduce = file_control.replace(".flagstat", ".bam")
        bam_file_reduce_ss = file_control.replace(".flagstat", ".subsample.bam")

        ratio = nr_of_reads_sample/nr_of_reads_control
        ratio = round(100*ratio)

    my_debug and print(ratio)

    link_command = f"ln -s {bam_file_keep} {bam_file_keep_ss}"
    subsample_command = f"samtools view --threads {threads} -bs {seed}.{ratio} {bam_file_reduce} > {bam_file_reduce_ss}"
    # samtools view -bs 10.84 $BAM_DIR/4_S4.proper_paired.bam > $BAM_DIR/4_S4.proper_paired.subsample.bam

    print(link_command)
    print(subsample_command)