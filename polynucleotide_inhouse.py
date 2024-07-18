#!/usr/bin/env python3

import pandas as pd
import argparse


# patterns = ["AAAA", "TTTT", "CCCC", "GGGG"]

# input_file = "/Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/merged_BSF-b_ORIs_alone_union500_nonoverlap50.window2500.fasta"

# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide
# input_files=$(find ./ -type f -name "*.window2000.fasta")
# input_files=$(find ./ -type f -name "*427-2018*.window2000.fasta")
# for input_file in $input_files; do echo $input_file; /Users/pmonsieurs/programming/leishmania_snsseq/bin/polynucleotide_inhouse.py --input ${input_file} --length 4; done

my_debug = 0

## Function to find specific patterns in a sequence and return their positions
def find_patterns(sequence):
    # patterns = ["AAAA", "TTTT", "CCCC", "GGGG"]
    positions = {pattern: [] for pattern in patterns}
    for pattern in patterns:
        start = 0
        while start < len(sequence):
            start = sequence.find(pattern, start)
            if start == -1:
                break
            positions[pattern].append(start)
            start += 1
    return positions

## Function to read FASTA file and find specific patterns in each sequence
def find_patterns_in_fasta(file_path):
    print("finding patterns in fasta")
    with open(file_path, 'r') as fasta_file:
        sequences = {}
        current_sequence = ""
        current_header = ""
        for line in fasta_file:
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = current_sequence
                current_header = line.strip()[1:]
                current_sequence = ""
            else:
                current_sequence += line.strip()
        # Process the last sequence in the file
        if current_header:
            sequences[current_header] = current_sequence
    
    pattern_positions = {}
    for header, sequence in sequences.items():
        pattern_positions[header] = find_patterns(sequence)
    return pattern_positions



if __name__ == '__main__':
    # Create an ArgumentParser object including the input option
    parser = argparse.ArgumentParser(description='extract polynucleotide (AAA, TTT, ...) of different lengths from a fasta file')
    parser.add_argument('--input', metavar='input_file', required=True, help='file containing the fasta sequences on which to do the analysis')
    parser.add_argument('--length', required=True, help='length of the poly-nucleotide')
    parser.add_argument('--window', required=False, default=2000, help='lenght of upstream and downstream region')


    ## extract the arguments
    args = parser.parse_args()
    input_file = args.input
    length = args.length
    window = args.window
    # output_dir = args.output_dir

    ## create the output file
    output_file=input_file.replace(".fasta", f".poly_{length}.csv")

    ## create the patterns to look for
    bases = ["A", "T", "G", "C"]
    patterns = [base * int(length) for base in bases]
    print(patterns)

    # Example usage
    pattern_positions = find_patterns_in_fasta(input_file)

    ## create an empty data frame that you can use to sum the 
    ## different positions
    df = pd.DataFrame(0, index=range(window*2), columns=patterns)


    count = 0
    for header, positions in pattern_positions.items():
        count = count + 1
        print(f"Header [{count}/{len(pattern_positions)}]: {header}", end = "\r")
        for pattern, positions_list in positions.items():
            my_debug and print(f"Pattern: {pattern}, Positions: {positions_list}")
            for position in positions_list:
                my_debug and print(position)
                if position < 4000: 
                    df.at[position, pattern] += 1

    print(df)
    df.to_csv(output_file)








