#!/usr/bin/env python3

import os

snsseq_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/data/ori_predictions/'
out_dir = '/Users/pmonsieurs/programming/trypanosoma_sofi_mining/results/ori/'

## set the length of the extension in bp
# extension = 500
extension = 2000

ori_files = os.listdir(snsseq_dir)


for ori_file in ori_files:
    
    if not ori_file.endswith(".bed"):
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
        print(line)
        data = line.split("\t")

        ## update start and end position
        new_start = max(1, int(data[1]) - extension)
        new_end = int(data[2]) + extension

        data[1] = str(new_start)
        data[2] = str(new_end)

        ## create new output line
        out_line = "\t".join(data)
        out_line = f"{out_line}\n"
        print(out_line)
        out_fh.write(out_line)

    ori_fh.close()
    out_fh.close()

