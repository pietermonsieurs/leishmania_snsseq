#!/usr/bin/env python3

import os

my_debug = 0

src_dir = '/Users/pmonsieurs/programming/leishmania_snsseq/'
# reference_dir = src_dir + 'data/refgenome/'
# reference_genome = reference_dir + 'TriTrypDB-46_TbruceiLister427_2018_Genome.fasta'

# extract the contigs and the corresponding lenghts from the bam-file
# rather than the genome reference. Bam file is basis for bed files, and 
# is certainly corresponding with the final bed output 
bwa_dir = src_dir + 'results/bwa/'
# bam_file = bwa_dir + '1_S1.subsample.bam'
bam_file = bwa_dir + '1_S1.proper_paired.subsample.bam'

# bin file for samtools
samtools = '/Users/pmonsieurs/programming/software/samtools-1.9/samtools'

temp_header_file = src_dir + 'results/header.temp'
sam_command = f"{samtools} view -H {bam_file} > {temp_header_file}"
os.system(sam_command)

contig_ids = []
contig_lenghts = {}

header_fh = open(temp_header_file, 'r')
for line in header_fh:
    line = line.rstrip()

    # only keep lines with contig information, the other ones can be 
    # ignored.
    if not line.startswith("@SQ"):
        continue

    data = line.split("\t")
    my_debug and print(data)
    contig_id = data[1].replace("SN:", "")
    contig_length = data[2].replace("LN:", "")
    print(f"contig [{contig_id}] -> {contig_length}")
    contig_ids.append(contig_id)
    contig_lenghts[contig_id] = int(contig_length)


#contig_list = "['" + "', '".join(contig_id) + "']"

# to be added to the file "GenomeData.py"
print(f"\n\ntbruc_chroms = {contig_ids}\n\n")
print(f"tbruc_chrom_lengths = {contig_lenghts}\n\n")



