## define the bedtools binary
bedtools_bin=/Users/pmonsieurs/programming/software/bedtools2-master/bin/bedtools

## first run the python script to extend the ORI bed file. For the 
## overlap this was only set to 500. Now, set to 2000bp in line with 
## the publication of Comoglio 2015 Cell Reports - Figure 2A
# /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_extend_ori.py


## calculate the overlap between G4Hunter output and the 
## extended bed file containing ORIs. Only report as output the 
## G4Hunter positions that are overlapping. 

## set intersect between different SNS-seq data of Bridlin
# g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter/
# snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/ori/

## set intersect between different shuffled SNS-seq data of Bridlin
g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter/
snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled/


## for Tb427 SNSseq data
g4_files=($(find $g4_dir -name "Tb427*.bed"))
snsseq_files=($(find $snsseq_dir -name "shuffeled*.extended_2000nt.bed"))

for g4_file in ${g4_files[@]}; do
    # echo $g4_file
    g4_file_short=$(basename ${g4_file})
    g4_file_short=${g4_file_short/.bed/}
    echo $g4_file_short
    for snsseq_file in ${snsseq_files[@]}; do
        # echo $snsseq_file
        snsseq_file_short=$(basename ${snsseq_file})
        snsseq_file_short=${snsseq_file_short/-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed/}
        echo $snsseq_file_short     
        output_file=${snsseq_dir}/${g4_file_short}.${snsseq_file_short}.bed
        echo $output_file   
        # overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        ${bedtools_bin} intersect -wb -a $g4_file -b $snsseq_file > ${output_file}
        
    done
done