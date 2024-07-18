## run seqkit locate with a pattern contain all types
## of (poly)nucleotides
seqkit_bin=/Users/pmonsieurs/programming/software/seqkit/seqkit

## create pattern file
cd /Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/seqkit/
vi pattern.fa

cd /Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/

## testrun seqkit
for fasta_file in *.fasta; do

    ## create output file name
    out_file_seqkit="${fasta_file/.fasta/.seqkit.txt}"
    out_file_seqkit="./seqkit/${out_file_seqkit}"
    
    echo $out_file_seqkit

    ## run seqtk
    # ${seqkit_bin} locate -P -f seqkit/pattern.fa ${fasta_file} -o $out_file_seqkit
    ${seqkit_bin} locate -P -f seqkit/pattern_mono.fa ${fasta_file} -o $out_file_seqkit

done
    

## make different subgroups
cd /Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/seqkit
for seqkit_file in *.seqkit.txt; do
    echo $seqkit_file

    out_file_seqkit_mono="${seqkit_file/.seqkit.txt/.seqkit.mono.txt}"
    out_file_seqkit_di="${seqkit_file/.seqkit.txt/.seqkit.di.txt}"
    out_file_seqkit_tri="${seqkit_file/.seqkit.txt/.seqkit.tri.txt}"
    out_file_seqkit_tetra="${seqkit_file/.seqkit.txt/.seqkit.tetra.txt}"

    grep "mono" $seqkit_file > ${out_file_seqkit_mono}
    grep "di" $seqkit_file > ${out_file_seqkit_di}
    grep "tri" $seqkit_file > ${out_file_seqkit_tri}
    grep "tetra" $seqkit_file > ${out_file_seqkit_tetra}

done

## if only the mono group is run, you can simply make symbolic links
cd /Users/pmonsieurs/programming/leishmania_snsseq/results/polynucleotide/seqkit
for seqkit_file in *.seqkit.txt; do
    echo $seqkit_file
    out_file_seqkit_mono="${seqkit_file/.seqkit.txt/.seqkit.mono.txt}"
    ln -s $seqkit_file ${out_file_seqkit_mono}
done