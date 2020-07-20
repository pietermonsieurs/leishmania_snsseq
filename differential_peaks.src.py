#!/usr/bin/env python3

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import argparse


my_debug = 1
src_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/'
results_dir = src_dir + 'results/'
depth_dir = results_dir + 'depth/'

plt.rcParams.update({'font.size': 50})




def make_plot (ax, position, depth, title, color):
    ax.plot(position, depth, label=title, color=color)
    # , label=sample, linestyle=line_style, color=line_color)
    ax.set_xlabel('Position in Genome')
    # ax.set_ylabel('Coverage')  
    ax.set_ylabel(title)  
    # ax.set_title(title)
    ax.grid(which="both")
    ax.fill_between(position, depth, color=color, alpha=0.50)
    # ax1.set_ylim(y_min, y_max)
    # ax1.legend()
    return ax


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description="Create an overview of the sequencing coverage for the SNS-seq samples, starting from chromosome, start and end")
    parser.add_argument('-c', '--chromosome', required=True, type=str,
                        help="chromosome ID")
    parser.add_argument('-s', '--start', required=True, type=int,
                        help="start of the region, excluding buffer")
    parser.add_argument('-e', '--end', required=True, type=int,
                        help="end of the region, excluding buffer")
    parser.add_argument('-d', '--subdir', required=True, type=str,
                        help="type of analysis where everything should be stored i.e. name of subdirectory in figure dir")
    


    args = parser.parse_args()
    chrom = args.chromosome
    start = args.start
    end = args.end
    sub_dir_suffix = args.subdir

    # create file names with depth data
    depth_S1_file = depth_dir + "1_S1.chrom_" + chrom + ".subsample.depth.csv"
    depth_S2_file = depth_dir + "2_S2.chrom_" + chrom + ".subsample.depth.csv"
    depth_S3_file = depth_dir + "3_S3.chrom_" + chrom + ".subsample.depth.csv"
    depth_S4_file = depth_dir + "4_S4.chrom_" + chrom + ".subsample.depth.csv"

    # create figure dir based on the type of analysis running e.g. S1_narrowPeak/
    figure_dir = depth_dir + sub_dir_suffix + "/"
    if not os.path.exists(figure_dir):
        mkdir_command = f"mkdir {figure_dir}"
        os.system(mkdir_command)
    


    # check if file is empty for sample S1 (and assume that 
    # the same goes up for the oterh samples: if S1 is ok, 
    # everthing is ok)
    if os.stat(depth_S1_file).st_size == 0:
        print(f"File {depth_S1_file} is empty --> {os.stat(depth_S1_file).st_size}!!")
        sys.exit()
    
    depth_S1 = pd.read_csv(depth_S1_file, sep="\t", header=None)
    depth_S2 = pd.read_csv(depth_S2_file, sep="\t", header=None)
    depth_S3 = pd.read_csv(depth_S3_file, sep="\t", header=None)
    depth_S4 = pd.read_csv(depth_S4_file, sep="\t", header=None)

    # calculate the difference between sample and corresponding control to have 
    # background corrected samples. 
    diff_S1 = depth_S1.iloc[:,2] - depth_S2.iloc[:,2]
    diff_S3 = depth_S3.iloc[:,2] - depth_S4.iloc[:,2]
    diff_S3_S1 = diff_S3 - diff_S1


    # define line types and colors. Not always needed (obsolet after 
    # update to a six-pane figure)
    line_styles = ['-', ':', '-', ':']
    line_colors = ['blue', 'blue', 'red', 'red']
      
    # create matplotlib figures / 6 plots
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(nrows=7, ncols=1,    figsize=(60, 60))

    # first four plots are the raw read counts for each sample
    # separately
    ax1 = make_plot(ax1, depth_S1.iloc[start:end,1], depth_S1.iloc[start:end,2], "cov S1", 'purple')
    ax2 = make_plot(ax2, depth_S2.iloc[start:end,1], depth_S2.iloc[start:end,2], "cov S2", 'blue')
    ax3 = make_plot(ax3, depth_S3.iloc[start:end,1], depth_S3.iloc[start:end,2], "cov S3", 'red')
    ax4 = make_plot(ax4, depth_S4.iloc[start:end,1], depth_S4.iloc[start:end,2], "cov S4", 'orange')



    # get minimum and maximum values to set those values in the plots, so that
    # the upper and lower graph have the same extreme values.
    y_max = max(diff_S1[start:end].max(), diff_S3[start:end].max())
    y_min = min(diff_S1[start:end].min(), diff_S3[start:end].min())
    y_max = y_max*1.05
    my_debug and print(f"min y {y_min} -- max y {y_max}")


    # fifth plot = depth of sample S1 (bg corrected)
    # ax5.plot(depth_S1.iloc[start:end,1], diff_S1.values[start:end])
    # # , label=sample, linestyle=line_style, color=line_color)
    # ax5.set_xlabel('Position in Genome')
    # ax5.set_ylabel('Background corrected coverage S1')  
    # ax5.grid(which="both")

    ax5 = make_plot(ax5, depth_S1.iloc[start:end,1], diff_S1.values[start:end], "cov S1 - bg S2", 'purple')
    ax5.set_ylim(y_min, y_max)
    # ax1.legend()


    # sixth plot = depth of sample S3 (bg corrected)
    # plt.ylim(y_min, y_max*1.05)
    # ax6.plot(depth_S3.iloc[start:end,1], diff_S3.values[start:end], color="red")
    # # , label=sample, linestyle=line_style, color=line_color)
    # ax6.set_xlabel('Position in Genome')
    # ax6.set_ylabel('Background corrected coverage S3')  
    # ax6.grid(which="both")

    ax6 = make_plot(ax6, depth_S3.iloc[start:end,1], diff_S3.values[start:end], "cov S3 - bg S4", 'red')
    ax6.set_ylim(y_min, y_max)

    ax7 = make_plot(ax7, depth_S3.iloc[start:end,1], diff_S3_S1.values[start:end], "cov S3 min S1", 'black')

    figure_file = figure_dir + chrom + "." + str(start) + "_" + str(end) + ".png"
    plt.savefig(figure_file)
    plt.close()




        # start = end

    # break

    # if counter > 10:
    #     break
    # else:
    #     continue

    # continue

    # for index in range(0,len(samples)):
    #     sample = samples[index]
    #     line_style = line_styles[index]
    #     line_color = line_colors[index]

    #     bam_file = bam_dir + sample + ".subsample.bam"
    #     depth_command = f"{samtools_bin} depth -a -r {chrom}:{start_region}-{stop_region} {bam_file} > {temp_depth_file}"
    #     my_debug and print(depth_command)
    #     os.system(depth_command)
    #     depth_data = pd.read_csv(temp_depth_file, sep="\t", header=None)
    #     my_debug and print(depth_data.head)

    #     df_list[sample] = depth_data

    #     my_debug and print(depth_data.iloc[:,1])
    #     ax1.plot(depth_data.iloc[:,1], depth_data.iloc[:,2], label=sample, linestyle=line_style, color=line_color)    
    
    # my_debug and print(df_list["1_S1"].shape)
    # my_debug and print(df_list["2_S2"].shape)
    # my_debug and print(df_list["3_S3"].shape)
    # my_debug and print(df_list["4_S4"].shape)

    # my_debug and print(df_list["1_S1"])
    # my_debug and print(df_list["2_S2"])
        

    # S1_diff = df_list["1_S1"].iloc[:,2] - df_list["2_S2"].iloc[:,2]
    # S3_diff = df_list["3_S3"].iloc[:,2] - df_list["4_S4"].iloc[:,2]

    # S1_diff_positive = S1_diff.copy()
    # S3_diff_positive = S3_diff.copy()

    # S1_diff_positive[S1_diff_positive < 0 ] = 0
    # S3_diff_positive[S3_diff_positive < 0 ] = 0

    # my_debug and print(f"S1_diff ==> {S1_diff}\n")
    # my_debug and print(f"S3_diff ==> {S3_diff}\n")
    # my_debug and print(df_list["1_S1"].iloc[:,2])
    


    # ax1.set_xlabel('Position in Genome')
    # ax1.set_ylabel('Depth of Coverage')  
    # ax1.legend()      

    # ax2.plot(df_list["1_S1"].iloc[:,1], S1_diff, label="S1_diff", color=line_colors[0])
    # ax2.plot(df_list["1_S1"].iloc[:,1], S3_diff, label="S3_diff", color=line_colors[2])
    # ax2.set_xlabel('Position in Genome')
    # ax2.set_ylabel('Depth of Coverage')  
    # ax2.legend()      

    # ax3.plot(df_list["1_S1"].iloc[:,1], S1_diff_positive, label="S1_diff", color=line_colors[0])
    # ax3.plot(df_list["1_S1"].iloc[:,1], S3_diff_positive, label="S3_diff", color=line_colors[2])
    # ax3.set_xlabel('Position in Genome')
    # ax3.set_ylabel('Depth of Coverage')  
    # ax3.legend() 

    # # plt.show()
    # plt.savefig(figure_file)
    # plt.close()

    # #if counter > max_counter:
    # #    break





