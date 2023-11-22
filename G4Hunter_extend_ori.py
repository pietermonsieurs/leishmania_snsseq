#!/usr/bin/env python3

import os

# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/ori_predictions/'
# out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/ori/'

# snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/ori_predictions_shuffled/'
# out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled/'

snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_ORIs_suffledORIs_G4-hunter_Mnase-seq/'
out_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/results/mnase_seq/'

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

