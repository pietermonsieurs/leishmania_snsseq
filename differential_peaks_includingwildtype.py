#!/usr/bin/env python3

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import argparse


my_debug = 1
src_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/'
results_dir = src_dir + 'results/'
# depth_dir = results_dir + 'depth/'
depth_dir = results_dir + 'depth_lessstringent/'


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
    depth_ctrl_file = depth_dir + "2_S2.chrom_" + chrom + ".subsample.depth.csv"
    depth_ctrlbg_file = depth_dir + "1_S1.chrom_" + chrom + ".subsample.depth.csv"
    depth_treat_file = depth_dir + "4_S4.chrom_" + chrom + ".subsample.depth.csv"
    depth_treatbg_file = depth_dir + "3_S3.chrom_" + chrom + ".subsample.depth.csv"

    depth_ctrl_wt_file = depth_dir + "5_S5.chrom_" + chrom + ".subsample.depth.csv"
    depth_ctrlbg_wt_file = depth_dir + "6_S6.chrom_" + chrom + ".subsample.depth.csv"
    depth_treat_wt_file = depth_dir + "7_S7.chrom_" + chrom + ".subsample.depth.csv"
    depth_treatbg_wt_file = depth_dir + "8_S8.chrom_" + chrom + ".subsample.depth.csv"


    # create figure dir based on the type of analysis running e.g. S1_narrowPeak/
    figure_dir = depth_dir + sub_dir_suffix + "/"
    if not os.path.exists(figure_dir):
        mkdir_command = f"mkdir {figure_dir}"
        os.system(mkdir_command)
    


    # check if file is empty for sample S1 (and assume that 
    # the same goes up for the oterh samples: if S1 is ok, 
    # everthing is ok)
    if os.stat(depth_ctrl_file).st_size == 0:
        print(f"File {depth_ctrl_file} is empty --> {os.stat(depth_S1_file).st_size}!!")
        sys.exit()
    
    depth_ctrl = pd.read_csv(depth_ctrl_file, sep="\t", header=None)
    depth_ctrlbg = pd.read_csv(depth_ctrlbg_file, sep="\t", header=None)
    depth_treat = pd.read_csv(depth_treat_file, sep="\t", header=None)
    depth_treatbg = pd.read_csv(depth_treatbg_file, sep="\t", header=None)

    depth_ctrl_wt = pd.read_csv(depth_ctrl_wt_file, sep="\t", header=None)
    depth_ctrlbg_wt = pd.read_csv(depth_ctrlbg_wt_file, sep="\t", header=None)
    depth_treat_wt = pd.read_csv(depth_treat_wt_file, sep="\t", header=None)
    depth_treatbg_wt = pd.read_csv(depth_treatbg_wt_file, sep="\t", header=None)

    # calculate the difference between sample and corresponding control to have 
    # background corrected samples. 
    diff_ctrl = depth_ctrl.iloc[:,2] - depth_ctrlbg.iloc[:,2]
    diff_treat = depth_treat.iloc[:,2] - depth_treatbg.iloc[:,2]
    # diff_treat_ctrl = diff_treat - diff_ctrl

    diff_ctrl_wt = depth_ctrl_wt.iloc[:,2] - depth_ctrlbg_wt.iloc[:,2]
    diff_treat_wt = depth_treat_wt.iloc[:,2] - depth_treatbg_wt.iloc[:,2]

      
    # create matplotlib figures / 4 plots
    # fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(nrows=7, ncols=1,    figsize=(60, 60))
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(60, 40))

    # # first four plots are the raw read counts for each sample
    # # separately
    ### ax1 = make_plot(ax1, depth_ctrl.iloc[start:end,1], depth_ctrl.iloc[start:end,2], "cov ctrl", 'purple')
    ### ax2 = make_plot(ax2, depth_ctrlbg.iloc[start:end,1], depth_ctrlbg.iloc[start:end,2], "cov ctrl bg", 'blue')
    ### ax3 = make_plot(ax3, depth_treat.iloc[start:end,1], depth_treat.iloc[start:end,2], "cov treat", 'red')
    ### ax4 = make_plot(ax4, depth_treatbg.iloc[start:end,1], depth_treatbg.iloc[start:end,2], "cov treat bg", 'orange')



    # get minimum and maximum values to set those values in the plots, so that
    # the upper and lower graph have the same extreme values.
    y_max = max(diff_ctrl[start:end].max(), 
                diff_treat[start:end].max(),
                diff_ctrl_wt[start:end].max(), 
                diff_treat_wt[start:end].max())
    y_min = min(diff_ctrl[start:end].min(), 
                diff_treat[start:end].min(),
                diff_ctrl_wt[start:end].min(), 
                diff_treat_wt[start:end].min())
    y_max = y_max*1.05
    my_debug and print(f"min y {y_min} -- max y {y_max}")


    # fifth plot = depth of sample S1 (bg corrected)
    # ax5.plot(depth_S1.iloc[start:end,1], diff_S1.values[start:end])
    # # , label=sample, linestyle=line_style, color=line_color)
    # ax5.set_xlabel('Position in Genome')
    # ax5.set_ylabel('Background corrected coverage S1')  
    # ax5.grid(which="both")

    ax1 = make_plot(ax1, depth_ctrl.iloc[start:end,1], diff_ctrl.values[start:end], "S2 - S1", 'purple')
    ax1.set_ylim(y_min, y_max)

    ax2 = make_plot(ax2, depth_treat.iloc[start:end,1], diff_treat.values[start:end], "S4 - S3", 'red')
    ax2.set_ylim(y_min, y_max)

    ax3 = make_plot(ax3, depth_ctrl_wt.iloc[start:end,1], diff_ctrl_wt.values[start:end], "S5 - S6", 'orange')
    ax3.set_ylim(y_min, y_max)

    ax4 = make_plot(ax4, depth_treat_wt.iloc[start:end,1], diff_treat_wt.values[start:end], "S7R - S8R", 'blue')
    ax4.set_ylim(y_min, y_max)

    figure_file = figure_dir + chrom + "." + str(start) + "_" + str(end) + ".including_wt.png"
    plt.savefig(figure_file)
    plt.close()


