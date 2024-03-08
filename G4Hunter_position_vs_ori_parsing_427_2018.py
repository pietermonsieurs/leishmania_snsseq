#!/usr/bin/env python3

import argparse
import os
import pandas as pd



## to be run in a for loop
# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/ori/
# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled/
# for overlap_file in Tb427*; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_position_vs_ori_parsing.py --input $PWD/$overlap_file; done

## for the Mnaseq
# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/mnase_seq
# for overlap_file in G4*merged*bed; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_position_vs_ori_parsing.py --input $PWD/$overlap_file; done
# for overlap_file in G4*shuffeled_seed*bed; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_position_vs_ori_parsing.py --input $PWD/$overlap_file; done

## for the 927 genome
# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/927/
# for overlap_file in 927_G4*bed; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_position_vs_ori_parsing.py --input $PWD/$overlap_file; done

## for the 427 genome
# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/427/
# for overlap_file in 427_G4*bed; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_position_vs_ori_parsing.py --input $PWD/$overlap_file; done

## for the 427_2018 DRIP-seq data: run only for files in which you
## are interested
# /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_position_vs_ori_parsing_427_2018.py --input /Users/pmonsieurs/programming/leishmania_snsseq/results/427_2018/427-2018_DRIP-seq_36.merged_PCF-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed.bed


my_debug = 0

def get_ori_position(ori_data, peak, g4_chrom):
    my_debug and print(f"peak --> {peak}")
    peak_hit = ori_data[ori_data.peak_id == peak]
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

    ## setting for the default ORI: 
    ori_dir = os.path.dirname(input_file)
    if "seed" in input_file:
        ori_info = input_file.split("_ORIs")[0]
        ori_info = ori_info.split("-seq_36.")[1]
        ori_file = f"{ori_dir}/{ori_info}_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed"
    else:
        ori_info = input_file.split("-")[-2]
        ori_info = ori_info.split("_36.")[-1]
        ori_file = f"{ori_dir}/{ori_info}-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed"

        
    print(ori_info)
    
    ## first line for the ori file of the original data, the second line for
    ## the Mnase-seq data, where oris have been rem
    # ori_file = f"{ori_dir}/{ori_info}-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed"
    # ori_file = f"{ori_dir}/{ori_info}_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed"
    # ori_file = f"{ori_dir}/{ori_info}_ORIs_alone_union500_nonoverlap50_woStrand.extended_2000nt.bed"

    ## setting for the randomly shuffled ORI:
    # shuffeled_seed668_PCF_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed
    # ori_dir = os.path.dirname(input_file)
    # ori_info = input_file.split(".")[-4]
    # print(ori_info)
    # ori_file = f"{ori_dir}/{ori_info}.extended_2000nt.bed"
    
    print(ori_file)

    ori_data = pd.read_csv(ori_file, sep = "\t", header=None)
    # ori_data.columns = ['chrom', 'start', 'end', 'peak_id', 'dummy1', 'dummy2', 'dummy3']
    ori_data.columns = ['chrom', 'start', 'end', 'peak_id', 'dummy1', 'dummy2', 'dummy3']
    print(ori_data)

    ## create coverage list with all zeros
    g4_locations_coverage = [0] * 2000 * 2
    g4_locations_peaks = [0] * 2000 * 2

    ## Loop over the input file and check for each of the G4 regions 
    ## the relative position versus the ORI
    print(input_file)
    input_fh = open(input_file, 'r')
    for line in input_fh:
        # print(line)
        line = line.rstrip()
        data = line.split("\t")
        
        ## extract start and end
        g4_chrom = data[0]
        g4_start = int(data[1])
        g4_end = int(data[2])
        g4_coverage = float(data[3]) ## not really G4, but MNase-seq coverage

        if (g4_coverage == 0): 
            # print("skip zeroes")
            continue

        ## get the list of ORI peaks it matches with, and check 
        ## the start and stop position
        # peak = data[3].split(",")
        # for peak in peaks:
        # peak = data[3]
        my_debug and  print(line)

        ## number 6 for original data, number 8 for mnaseq data, 7 for experimental g4 (marisco)
        # peak = data[6]
        # peak = data[8]
        peak = data[7]

        my_debug and print(f"peak --> {peak}")
        ori_peak = get_ori_position(ori_data, peak, g4_chrom)

        ## calculate the relative position
        rel_g4_start = min(max(-2000, g4_start - ori_peak),int(len(g4_locations_coverage)/2)-1) 

        rel_g4_end = min(g4_end - ori_peak, int(len(g4_locations_coverage)/2)-1)
        rel_g4_peak = int((rel_g4_start + rel_g4_end)/2)
        my_debug and print(f"G4 start {g4_start} -- G4 end {g4_end} -- relative G4 start peak {rel_g4_start} -- relative G4 end peak {rel_g4_end} -- ori peak {ori_peak}")


        ## store in list
        for i in range(rel_g4_start, rel_g4_end + 1):
            # print(i)
            g4_locations_coverage[2000 + i] = g4_locations_coverage[2000 + i] + g4_coverage

        g4_locations_peaks[2000 + rel_g4_peak] = g4_locations_peaks[2000 + rel_g4_peak] + 1

    
    
    # print(g4_locations_coverage)
    print(g4_locations_peaks)
        
    ## write to output file handle    
    output_fh = open(output_file, 'w') 
    for i in range(0,len(g4_locations_coverage)):
        line_out = f"{i-2000},{g4_locations_coverage[i]}\n"
        print(line_out)
        output_fh.write(line_out)
    