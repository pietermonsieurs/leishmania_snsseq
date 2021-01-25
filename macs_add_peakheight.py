#!/usr/bin/env python3

import pandas as pd
import os

### module load SAMtools # --> should be integrated into the bash script for qsub

src_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/'
results_dir = src_dir + 'results/'
macs_dir = results_dir + 'macs_switch/'
bwa_dir = results_dir + 'bwa/'
temp_dir = results_dir + 'temp/'

# sample = '2_S2'
# bg = '1_S1'

sample = '4_S4'
bg = '3_S3'


peaks_file = macs_dir + sample + '_peaks.narrowPeak'
peaks_file_out = macs_dir + sample + '_peaks_withmaxheigth.narrowPeak'
bam_file = bwa_dir + sample + '.mapq30.removedups.proper_paired.subsample.bam'
bam_file_bg = bwa_dir + bg + '.mapq30.removedups.proper_paired.subsample.bam'
temp_file = temp_dir + 'samtools_depth_temp.csv'
temp_file_bg = temp_dir + 'samtools_depth_bg_temp.csv'


peaks_df = pd.read_csv(peaks_file, sep="\t", header=None)
print(peaks_df.head)

max_height = []
for index, row in peaks_df.iterrows():
    chrom = row[0]
    start = row[1]
    end = row[2]
    print(f"{chrom} :: {start} -> {end}")


    depth_command = f"samtools depth {bam_file} -r {chrom}:{start}-{end} > {temp_file}"
    os.system(depth_command)

    depth_command_bg = f"samtools depth {bam_file_bg} -r {chrom}:{start}-{end} > {temp_file_bg}"
    os.system(depth_command_bg)

    depth_data = pd.read_csv(temp_file, sep="\t", header=None)
    depth_data_bg = pd.read_csv(temp_file_bg, sep="\t", header=None)
    diff_depth = depth_data.iloc[:,2] - depth_data_bg.iloc[:,2]
    max_peak = max(diff_depth)
    print(f"max_peak = {max_peak}")

    max_height.append(max_peak)
    
    #if (index > 10):
    #    break
    print(index)

print(max_height)
peaks_df['max_height'] = max_height
peaks_df.to_csv(peaks_file_out)