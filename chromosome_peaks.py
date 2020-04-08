#!/usr/bin/env python3

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

my_debug = 0
# src_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/'
src_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/'
results_dir = src_dir + 'results/'
depth_dir = results_dir + 'depth/'
bam_dir = results_dir + 'bwa/'

# file with all chromosome names / contigs:
chrom_file = depth_dir + 'list_of_chromosomes.csv'
chrom_data = pd.read_csv(chrom_file, sep="\t", header=None, usecols=[0])
chrom_list = chrom_data[0].to_list()
chrom_list = chrom_list[:-1]
my_debug and print(chrom_list)

# diffpeaks_out_file = macs_dir + diff_type + '.csv'
figure_dir = depth_dir + "figures/"



df_list = {}
counter = 0

for chrom in chrom_list:
    print(chrom)
    counter = counter + 1

    if counter < 259:
        continue
    
    # define line types and colors
    line_styles = ['-', ':', '-', ':']
    line_colors = ['blue', 'blue', 'red', 'red']

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, 
                                        figsize=(16, 6))


    depth_S1_file = depth_dir + "1_S1.chrom_" + chrom + ".subsample.depth.csv"
    depth_S2_file = depth_dir + "2_S2.chrom_" + chrom + ".subsample.depth.csv"
    depth_S3_file = depth_dir + "3_S3.chrom_" + chrom + ".subsample.depth.csv"
    depth_S4_file = depth_dir + "4_S4.chrom_" + chrom + ".subsample.depth.csv"

    depth_S1 = pd.read_csv(depth_S1_file, sep="\t", header=None)
    depth_S2 = pd.read_csv(depth_S2_file, sep="\t", header=None)
    depth_S3 = pd.read_csv(depth_S3_file, sep="\t", header=None)
    depth_S4 = pd.read_csv(depth_S4_file, sep="\t", header=None)

    # calculate the difference between sample and corresponding control to have 
    # background corrected samples. 
    diff_S1 = depth_S1.iloc[:,2] - depth_S2.iloc[:,2]
    diff_S3 = depth_S3.iloc[:,2] - depth_S4.iloc[:,2]

    # get minimum and maximum values to set those values in the plots, so that
    # the upper and lower graph have the same extreme values.
    y_max = max(diff_S1.max(), diff_S3.max())
    y_min = min(diff_S1.min(), diff_S3.min())
    y_max = y_max*1.05
    # print(f"min y {y_min} -- max y {y_max}")

    
    # upper plot = depth of sample S1 (bg corrected)
    ax1.plot(depth_S1.iloc[:,1], diff_S1.values)
    # , label=sample, linestyle=line_style, color=line_color)
    ax1.set_xlabel('Position in Genome')
    ax1.set_ylabel('Background corrected coverage S1')  
    ax1.grid(which="both")
    ax1.set_ylim(y_min, y_max)
    # ax1.legend()


    # lower plot = depth of sample S3 (bg corrected)
    plt.ylim(y_min, y_max*1.05)
    ax2.plot(depth_S3.iloc[:,1], diff_S3.values, color="red")
    # , label=sample, linestyle=line_style, color=line_color)
    ax2.set_xlabel('Position in Genome')
    ax2.set_ylabel('Background corrected coverage S3')  
    ax2.grid(which="both")
    ax2.set_ylim(y_min, y_max)

    figure_file = figure_dir + chrom + ".png"
    plt.savefig(figure_file)
    plt.close()

    # if counter > 10:
    #     break
    # else:
    #     continue

    continue

    for index in range(0,len(samples)):
        sample = samples[index]
        line_style = line_styles[index]
        line_color = line_colors[index]

        bam_file = bam_dir + sample + ".subsample.bam"
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





