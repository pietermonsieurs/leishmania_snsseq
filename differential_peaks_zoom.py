#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt

my_debug = 0
src_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/'
results_dir = src_dir + 'results/'
macs_dir = results_dir + 'macs/'
# figure_dir = macs_dir + "figures/"
bam_dir = results_dir + 'bwa/'

diff_type = 'only_present_in_S3'
# diff_type = 'only_present_in_S1'
# diff_type = 'only_present_in_S3_broad'
# diff_type = 'only_present_in_S1_broad'


diffpeaks_out_file = macs_dir + diff_type + '.csv'
figure_dir = macs_dir + "figures/" + diff_type + "/"

# program directory
samtools_bin = '/Users/pmonsieurs/programming/software/samtools-1.9/samtools'

# parameter settings. window = how many flink nucleotides need to
# be taken.
window = 200
samples = ['1_S1', '2_S2', '3_S3', '4_S4']
temp_depth_file = macs_dir + 'temp.csv'

data = pd.read_csv(diffpeaks_out_file, sep="\t", header=None)
# ['chrom', 'start', 'stop', 'peak_id'])
print(data.head)

max_count = 0
counter = 0

df_list = {}

for index, row in data.iterrows():

    counter = counter + 1
    
    # extract useful parameter from the csv file
    chrom = row[0]
    start = row[1]
    stop = row[2]
    print(f"making plot for chrom {chrom} :: {start} -> {stop}")
    figure_code = f"{chrom}__{start}_{stop}"
    figure_file = figure_dir + figure_code + ".png"

    # extra coverage data for region of interest. 
    start_region = start - window
    stop_region = stop + window

    # define line types and colors
    line_styles = ['-', ':', '-', ':']
    line_colors = ['blue', 'blue', 'red', 'red']

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, 
                                        figsize=(12, 4))

    for index in range(0,len(samples)):
        sample = samples[index]
        line_style = line_styles[index]
        line_color = line_colors[index]

        bam_file = bam_dir + sample + ".proper_paired.subsample.bam"
        depth_command = f"{samtools_bin} depth -a -r {chrom}:{start_region}-{stop_region} {bam_file} > {temp_depth_file}"
        my_debug and print(depth_command)
        os.system(depth_command)
        depth_data = pd.read_csv(temp_depth_file, sep="\t", header=None)
        my_debug and print(depth_data.head)

        df_list[sample] = depth_data

        my_debug and print(depth_data.iloc[:,1])
        ax1.plot(depth_data.iloc[:,1], depth_data.iloc[:,2], label=sample, linestyle=line_style, color=line_color)    
    
    my_debug and print(df_list["1_S1"].shape)
    my_debug and print(df_list["2_S2"].shape)
    my_debug and print(df_list["3_S3"].shape)
    my_debug and print(df_list["4_S4"].shape)

    my_debug and print(df_list["1_S1"])
    my_debug and print(df_list["2_S2"])
        

    S1_diff = df_list["1_S1"].iloc[:,2] - df_list["2_S2"].iloc[:,2]
    S3_diff = df_list["3_S3"].iloc[:,2] - df_list["4_S4"].iloc[:,2]

    S1_diff_positive = S1_diff.copy()
    S3_diff_positive = S3_diff.copy()

    S1_diff_positive[S1_diff_positive < 0 ] = 0
    S3_diff_positive[S3_diff_positive < 0 ] = 0

    my_debug and print(f"S1_diff ==> {S1_diff}\n")
    my_debug and print(f"S3_diff ==> {S3_diff}\n")
    my_debug and print(df_list["1_S1"].iloc[:,2])
    


    ax1.set_xlabel('Position in Genome')
    ax1.set_ylabel('Depth of Coverage')  
    ax1.legend()      

    ax2.plot(df_list["1_S1"].iloc[:,1], S1_diff, label="S1_diff", color=line_colors[0])
    ax2.plot(df_list["1_S1"].iloc[:,1], S3_diff, label="S3_diff", color=line_colors[2])
    ax2.set_xlabel('Position in Genome')
    ax2.set_ylabel('Depth of Coverage')  
    ax2.legend()      

    ax3.plot(df_list["1_S1"].iloc[:,1], S1_diff_positive, label="S1_diff", color=line_colors[0])
    ax3.plot(df_list["1_S1"].iloc[:,1], S3_diff_positive, label="S3_diff", color=line_colors[2])
    ax3.set_xlabel('Position in Genome')
    ax3.set_ylabel('Depth of Coverage')  
    ax3.legend() 

    # plt.show()
    plt.savefig(figure_file)
    plt.close()

    #if counter > max_counter:
    #    break





