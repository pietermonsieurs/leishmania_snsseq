#!/usr/bin/env python3

import argparse


## to be run in a for loop
# cd /Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter
# for csv_file in *.csv; do /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_convert_csv2bed.py --input $csv_file; done

my_debug = 1

def create_outputfile(input_file):
    output_file_plus = input_file.replace(".csv", ".plus.bed")
    output_file_min = input_file.replace(".csv", ".min.bed")

    return output_file_plus, output_file_min

def convert2bed(input_file, output_file_plus, output_file_min):

    ## open files for reading and writing
    in_fh = open(input_file, 'r')
    out_fh_plus = open(output_file_plus, 'w')
    out_fh_min = open(output_file_min, 'w')


    chrom = ""
    for line in in_fh:
        line = line.rstrip()

        ## check if we begin with a new chromosome
        if line.startswith(">"):
            chrom = line.replace(">", "")
            my_debug and print(f"new chromosome {chrom} found!")
            continue
        elif line.startswith("Start"):
            ## header line so skip
            continue

        data = line.split("\t")
        data = [x.replace(" ", "") for x in data]

        ## create output line
        line_out = "\t".join(data[:2])
        bed_line_out = f"{chrom}\t{line_out}\n"

        ## split between + and - strand based on the sign (positive or
        ## negative) of the score: + = G-rich = + strand, - = C-rich = 
        ## - strand
        if float(data[4]) > 0:
            ## write to positive strand output
            out_fh_plus.write(bed_line_out)
        elif float(data[4]) < 0:
            ## write to negative strand output
            out_fh_min.write(bed_line_out)
        ## if not a line announcing a new chromosome, extract
        ## the start and end information
    
    ## close filehandles
    out_fh_plus.close()
    out_fh_min.close()
    in_fh.close()

if __name__ == '__main__':


    # Create an ArgumentParser object including the input option
    parser = argparse.ArgumentParser(description='Read an input file from the command line')
    parser.add_argument('--input', metavar='input_file', required=True, help='Path to the input file')
    args = parser.parse_args()

    ## extract the input file
    input_file = args.input

    # Call the main function with the input file path
    output_file_plus, output_file_min = create_outputfile(input_file)
    convert2bed(input_file, output_file_plus, output_file_min)
