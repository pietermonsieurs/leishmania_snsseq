#!/usr/bin/env python3

import os
import glob 
import argparse
from Bio import SeqIO


genome_fasta_file = '/Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta'


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
    parser.add_argument('--window', required=False, default=2500, help='lenght of upstream and downstream region')
    parser.add_argument('--output_dir', required=False, default="/Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/", help='output directory to write the output file')

 
    args = parser.parse_args()

    ## extract the input and extract up and downstream region
    input_file = args.input
    window = args.window
    output_dir = args.output_dir
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
        center = int(int(data[2]) - int(data[1])/2)
        start = center - window
        end = center + window

        ## get the sequence
        sequence = extract_subsequence(chroms, chrom, start, end)

        fasta_line = f">{chrom}|center:{center}|chrom:{chrom}|{start}-{end}\n{sequence}\n"
        out_fh.write(fasta_line)


    ## file handles closing
    input_fh.close()
    out_fh.close()





