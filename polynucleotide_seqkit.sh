## run seqkit locate with a pattern contain all types
## of (poly)nucleotides
seqkit_bin=/Users/pmonsieurs/programming/software/seqkit/seqkit

## create pattern file
cd /Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/seqkit/
vi pattern.fa


## testrun seqkit
for fasta_file in *.fasta; do

    ## create output file name
    out_file_seqkit="${fasta_file/.fasta/.seqkit.txt}"
    out_file_seqkit="./seqkit/${out_file_seqkit}"
    
    echo $out_file_seqkit

    ## run seqtk
    ${seqkit_bin} locate -P -f seqkit/pattern.fa ${fasta_file} -o $out_file_seqkit

done
    

## make different subgroups
cd /Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/seqkit
for seqkit_file in *.seqkit.txt; do
    out_file_seqkit_mono="${fasta_file/.seqkit.txt/.seqkit.mono.txt}"
    out_file_seqkit_di="${fasta_file/.seqkit.txt/.seqkit.di.txt}"
    out_file_seqkit_tri="${fasta_file/.seqkit.txt/.seqkit.tri.txt}"
    out_file_seqkit_tetra="${fasta_file/.seqkit.txt/.seqkit.tetra.txt}"

    grep "mono"  $seqkit_file > ${out_file_seqkit_mono}
    grep "di"  $seqkit_file > ${out_file_seqkit_di}
    grep "tri"  $seqkit_file > ${out_file_seqkit_tri}
    grep "tetra" $seqkit_file > ${out_file_seqkit_tetra}

done