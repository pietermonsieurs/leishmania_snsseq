#!/usr/bin/env python3

import os
import glob 
import argparse
from Bio import SeqIO


## to be run for all .bed file using the find command combined
## with xargs
# input_files=$(find ./ -type f -name "*.bed")
# for input_file in $input_files; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/polynucleotide_get_sequences.py --input ${input_file}; done


## for the final files of Tb427 (NOT Tb427_2018)
# cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data
# cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data_set2/
# ref_genome=/Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-62_TbruceiLister427_Genome.fasta
# input_files=$(find ./ -type f -name "*ORI*.bed")
# for input_file in $input_files; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/polynucleotide_get_sequences.py --input ${input_file} --ref_genome=${ref_genome}; done

## for the final files of Tb927 
# cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data
# cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data_set2/
# ref_genome=/Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-63_TbruceiTREU927_Genome.fasta
# input_files=$(find ./ -type f -name "*ORI*.bed")
# for input_file in $input_files; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/polynucleotide_get_sequences.py --input ${input_file} --ref_genome=${ref_genome}; done

## for the final files of Tb427_2018
# cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427-2018_data
# cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427-2018_data_set2/
# ref_genome=/Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta
# input_files=$(find ./ -type f -name "*.bed")
# for input_file in $input_files; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/polynucleotide_get_sequences.py --input ${input_file} --ref_genome=${ref_genome}; done


def read_genome(fasta_file):
    chroms = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        chroms[record.id] = record
    return chroms   


def extract_subsequence(chroms, chromosome, start, end):
    record = chroms[chromosome]
    sequence = record.seq[start - 1:end]
    return sequence 


if __name__ == '__main__':
    # Create an ArgumentParser object including the input option
    parser = argparse.ArgumentParser(description='extract the flaking regions of the origins predicted using SNS-seq or random shuffled data')
    parser.add_argument('--input', metavar='input_file', required=True, help='file containing the ori predictions')
    parser.add_argument('--window', required=False, default=2000, help='lenght of upstream and downstream region')
    parser.add_argument('--output_dir', required=False, default="/Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/", help='output directory to write the output file')
    parser.add_argument('--ref_genome', required=True, default='/Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta', help='Genome fasta file from which you want to cut the sequences upstream and downstream of the ORIs')
 
    args = parser.parse_args()

    ## extract the input and extract up and downstream region
    input_file = args.input
    window = args.window
    output_dir = args.output_dir
    genome_fasta_file = args.ref_genome

    print(output_dir)


    ## get all chromosome sequences. This goes faster than parsing 
    ## the genome fasta file everytime again
    chroms = read_genome(genome_fasta_file)

    input_fh = open(input_file, 'r')
    
    ## create the output file
    output_file = os.path.basename(input_file)
    output_file = output_file.replace(".bed", f".window{window}.fasta")
    output_file = f"{output_dir}/{output_file}"
    out_fh = open(output_file, 'w')

    for line in input_fh:
        line = line.rstrip()
        data = line.split("\t")

        ## extract relevant information
        chrom = data[0]
        center = int((int(data[2]) + int(data[1]))/2)
        start = center - window
        end = center + window

        print(f"starting from {chrom} - {data[1]} - {data[2]} :: center {center} --> extracting {start} - {end}")

        ## get the sequence
        sequence = extract_subsequence(chroms, chrom, start, end)

        fasta_line = f">{chrom}|center:{center}|chrom:{chrom}|{start}-{end}\n{sequence}\n"
        out_fh.write(fasta_line)


    ## file handles closing
    input_fh.close()
    out_fh.close()





