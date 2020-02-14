#!/usr/bin/env python3

import os

my_debug = 1
src_dir = '/user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/'
data_dir = src_dir + 'data/'
out_dir = src_dir + 'results/'
ref_genome = data_dir + 'refgenome/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta'
ref_genome_index = ref_genome + ".fai"
threads = 20

index_command = f"bwa index {ref_genome}"
print(index_command)

for fq1 in os.listdir(data_dir):

    # only select the first reads, if not a file containing the 
    # first read, skip to next one
    if not fq1.endswith("R1_001.fastq.gz"):
        continue

    fq2 = fq1.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
    out_file = fq1.replace("_L001_R1_001.fastq.gz", ".bam")

    fq1 = data_dir + fq1
    fq2 = data_dir + fq2
    out_file = out_dir + "bwa/" + out_file 
    
    bwa_command = f"bwa mem -t {threads} {ref_genome} {fq1} {fq2} | samtools sort -@{threads} -o {out_file} -"
    print(bwa_command)
    

