#!/usr/bin/env python3

import os
import glob

my_debug = 0
g4hunter_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/ori/'
g4_bed_files = glob.glob(f"{g4hunter_dir}/Tb427*.bed")
g4_bed_files.sort()

for bed_file in g4_bed_files: 
    my_debug and print(f"{bed_file}")

    ## extract all the metadata from the file name
    bed_file_name = bed_file.split("/")
    bed_file_name = bed_file_name[len(bed_file_name)-1]
    my_debug and print(bed_file_name)

    ## extract the sample name
    sample_name = bed_file_name.split(".")
    sample_name = sample_name[len(sample_name) -2]
    sample_name = sample_name.split("_")[1]
    my_debug and print(sample_name)

    ## extract the strand
    strand = bed_file_name.split(".")
    strand = strand[len(strand)-3]
    my_debug and print(strand)

    ## extract G4 config 
    G4_config = bed_file_name.split(".")
    G4_config = f"{G4_config[0]}.{G4_config[1]}"
    my_debug and print(G4_config)


    

    ## open the file and count the number of G4 per ORI peak
    ori_count = dict()
    fh = open(bed_file, 'r')
    for line in fh: 
        data = line.split("\t")
        ori = data[6]
        # oris = data[6].split(",")
        # for ori in oris:
        if ori in ori_count:
            ori_count[ori] = ori_count[ori] + 1 
        else:
            ori_count[ori] = 1
    nr_of_oris_with_G4 = (len(ori_count))

    ## check the total number of ORIs
    ori_source_file = f"{g4hunter_dir}/merged_{sample_name}-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed"
    ori_fh = open(ori_source_file, 'r')
    nr_of_oris = len(ori_fh.readlines())

    ## print the output
    perc_with_g4 = int(100*nr_of_oris_with_G4/nr_of_oris)
    print(f"{sample_name}\t{G4_config}\t{strand}\t{nr_of_oris}\t{nr_of_oris_with_G4}\t{perc_with_g4}%")

