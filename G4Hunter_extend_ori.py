#!/usr/bin/env python3

import os

# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/ori_predictions/'
# out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/ori/'

# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/ori_predictions_shuffled/'
# out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled/'

##### ADDITIONAL SETTINGS ####

## settings for the data related to strain 927 - normal, with normal and
## shuffled ORI dir (i.e. snsseq_dir) 
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data/927_ORIs_bed/'
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data/927_shuffledORIs-bed/'
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data_set2/927_ORIs_bed/'
snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data_set2/927_shuffledORIs-bed/'
out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/927/'

## settings for the data related to strain 427 - normal, with normal and
## shuffled ORI dir (i.e. snsseq_dir) 
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data/427_ORIs_bed/'
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data/427_shuffledORIs-bed/'
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data_set2/427_ORIs_set2_bed/'
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data_set2/427_shuffledORIs_set2_bed/'
# out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/427/'

## settings for the data related to strain 427_2018 
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427-2018_data_set2/427-2018_ORIs_set2_bed/'
# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427-2018_data_set2/427-2018_shuffledORIs_set2-bed/'
# out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/427_2018/'

## settings for the data related to strain 427_2018 - normal: this was already 
## analysed as this is the reference genome also used by Bridlin to caclulate
## the SNS-seq data

## set the length of the extension in bp
# extension = 500
extension = 2000

ori_files = os.listdir(snsseq_dir)


for ori_file in ori_files:
    
    if not ori_file.endswith(".bed") or not "ORI" in ori_file:
        continue

    out_file = ori_file.replace('.bed', f'.extended_{extension}nt.bed')
    print(f"{ori_file} -> {out_file}")

    ## create file names with full path
    ori_file_full = f"{snsseq_dir}/{ori_file}"
    out_file = f"{out_dir}/{out_file}"

    ## open files for reading and writing respectively
    ori_fh = open(ori_file_full, 'r')
    out_fh = open(out_file, 'w')

    for line in ori_fh:
        line = line.rstrip()
        # print(line)
        data = line.split("\t")

        ## update start and end position
        center = int((int(data[2]) + int(data[1]))/2)
        # print(center)
        new_start = max(1, int(center) - extension)
        new_end = int(center) + extension

        data[1] = str(new_start)
        data[2] = str(new_end)

        ## create new output line
        out_line = "\t".join(data)
        out_line = f"{out_line}\n"
        # print(out_line)
        out_fh.write(out_line)

    ori_fh.close()
    out_fh.close()

