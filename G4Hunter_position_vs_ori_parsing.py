#!/usr/bin/env python3

import argparse
import os
import pandas as pd



## to be run in a for loop
# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled/
# for overlap_file in Tb427*; do echo $overlap_file; done


my_debug = 0

def get_ori_position(ori_data, peak, g4_chrom):
    # print(ori_file)
    peak_hit = ori_data[ori_data.peak_id == peak]
    # print(peak_hit)
    my_debug and print(peak_hit)
    ori_start = peak_hit.start.iloc[0]
    ori_end = peak_hit.end.iloc[0]
    ori_peak = int((ori_start + ori_end)/2)
    my_debug and print(f"{ori_start} --> {ori_end}: peak at {ori_peak}")
    return ori_peak

# def convert2bed(input_file, output_file):


if __name__ == '__main__':


    # Create an ArgumentParser object including the input option
    parser = argparse.ArgumentParser(description='starting from the G4 quadruplexes overlapping with the ORI regions, find the relative position versus this ORI')
    parser.add_argument('--input', metavar='input_file', required=True, help='Path to the input file containing G4 quadruxplexes')
    # parser.add_argument('--ori', metavar='ori_file', required=True, help='Path to file with the predicted ORI locations')
    
    args = parser.parse_args()

    ## extract the input and create the ORI file based on the information
    ## stored in the input file
    input_file = args.input
    print(input_file)

    ## create output file
    output_file = input_file.replace(".bed", ".cov")

    ## extract ORI information. you have to go back here to the original ORI file,
    ## which means before the overlap with the G4 hunter positions was caclulated
    ori_dir = os.path.dirname(input_file)
    ori_info = input_file.split(".")[-2]
    ## setting for the default ORI: 
    ori_file = f"{ori_dir}/{ori_info}-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed"
    ## setting for the randomly shuffled ORI:
    
    print(ori_file)

    ori_data = pd.read_csv(ori_file, sep = "\t", header=None)
    # ori_data.columns = ['chrom', 'start', 'end', 'peak_id', 'dummy1', 'dummy2', 'dummy3']
    ori_data.columns = ['chrom', 'start', 'end', 'peak_id', 'dummy1', 'dummy2', 'dummy3']

    ## create coverage list with all zeros
    g4_locations_coverage = [0] * 2000 * 2
    g4_locations_peaks = [0] * 2000 * 2

    ## Loop over the input file and check for each of the G4 regions 
    ## the relative position versus the ORI
    input_fh = open(input_file, 'r')
    for line in input_fh:
        # print(line)
        line = line.rstrip()
        data = line.split("\t")
        
        ## extract start and end
        g4_chrom = data[0]
        g4_start = int(data[1])
        g4_end = int(data[2])

        ## get the list of ORI peaks it matches with, and check 
        ## the start and stop position
        # peak = data[3].split(",")
        # for peak in peaks:
        # peak = data[3]
        peak = data[6]
        my_debug and print(peak)
        ori_peak = get_ori_position(ori_data, peak, g4_chrom)

        ## calculate the relative position
        rel_g4_start = min(max(-2000, g4_start - ori_peak),int(len(g4_locations_coverage)/2)-1) 

        rel_g4_end = min(g4_end - ori_peak, int(len(g4_locations_coverage)/2)-1)
        rel_g4_peak = int((rel_g4_start + rel_g4_end)/2)
        print(f"G4 start {g4_start} -- G4 end {g4_end} -- relative G4 start peak {rel_g4_start} -- relative G4 end peak {rel_g4_end} -- ori peak {ori_peak}")


        ## store in list
        for i in range(rel_g4_start, rel_g4_end + 1):
            # print(i)
            g4_locations_coverage[2000 + i] = g4_locations_coverage[2000 + i] + 1

        g4_locations_peaks[2000 + rel_g4_peak] = g4_locations_peaks[2000 + rel_g4_peak] + 1

    
    
    # print(g4_locations_coverage)
    print(g4_locations_peaks)
        
    ## write to output file handle    
    output_fh = open(output_file, 'w') 
    for i in range(0,len(g4_locations_coverage)):
        line_out = f"{i-2000},{g4_locations_coverage[i]}\n"
        print(line_out)
        output_fh.write(line_out)
    