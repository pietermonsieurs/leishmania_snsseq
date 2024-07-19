## define the bedtools binary
#bedtools_bin=/Users/pmonsieurs/programming/software/bedtools2-master/bin/bedtools
bedtools_bin=/Users/pmonsieurs/programming/software/bedtools2/bin/bedtools

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
# g4_files=($(find $g4_dir -maxdepth 1  -name "Tb427*.bed"))
# snsseq_files=($(find $snsseq_dir -name "merged*.extended_2000nt.bed"))

## set intersect between different shuffled SNS-seq data of Bridlin, but now
## with the shuffled ORI sequences
# g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter/
# snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled/
# g4_files=($(find $g4_dir -maxdepth 1  -name "Tb427*.bed"))
# snsseq_files=($(find $snsseq_dir -name "shuffeled*.extended_2000nt.bed"))

## set intersect between the G4 Hunter data with the ORI sequences
## either true ones or shuffled ones. Different genetic background
## (Tb427 vs Tb427_2018), but combining real and shuffled in one 
## analysis
# g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_ORIs_suffledORIs_G4-hunter_Mnase-seq
# snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/mnase_seq/
# g4_files=($(find $g4_dir -maxdepth 1  -name "G4*.bed"))
# snsseq_files=($(find $snsseq_dir -name "*ORI*extended*.bed"))





## for the 427_2018 genetic background on the new data set (set2). G4 Hunter data 
## need to come from your own results as they are not provided by Bridlin but have
## been run by you for the first time
g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter/
snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/427_2018/
snsseq_files=($(find $snsseq_dir -name "*ORI*extended*.bed"))
g4_files=($(find $g4_dir -maxdepth 1  -name "Tb427*.bed"))




## set intersect between the G4 experimental data with the ORI sequences
## either true ones or shuffled ones. Same information as the initial 
## picture, only difference is that now experimental data (Marisco) are
## used instead of the G4Hunter data. First need to do conversion from 
## bigwig to bed
bigwig2bed=/Users/pmonsieurs/programming/software/bigWigToBedGraph/bigWigToBedGraph
# g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data/927_G4_experimental_Marsico/
g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data_set2/927_G4_experimental_Marsico/
for bw_file in ${g4_dir}/*.bw; do
    bed_file=${bw_file/.bw/.bed}
    $bigwig2bed $bw_file $bed_file
done
g4_files=($(find $g4_dir -maxdepth 1  -name "*G4*.bed"))
snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/927/
snsseq_files=($(find $snsseq_dir -name "*ORI*extended*.bed"))


## print the file lists
echo ${g4_files[@]}
echo ${snsseq_files[@]}


## calculate the overlap
for g4_file in ${g4_files[@]}; do
    echo $g4_file
    g4_file_short=$(basename ${g4_file})
    g4_file_short=${g4_file_short/.bed/}
    echo $g4_file_short 
    for snsseq_file in ${snsseq_files[@]}; do
        echo " --> ${snsseq_file}"
        snsseq_file_short=$(basename ${snsseq_file})
        echo " --> basename ${snsseq_file_short}"
        
        ## do substitution depending on the type of file
        # snsseq_file_short=${snsseq_file_short/-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed/}
        # snsseq_file_short=${snsseq_file_short/_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed/}
        # snsseq_file_short=${snsseq_file_short/_ORIs_alone_union500_nonoverlap50_woStrand.extended_2000nt.bed/}
        snsseq_file_short=${snsseq_file_short/.extended_2000nt.bed/}
        echo " --> snsseq_file_short ${snsseq_file_short}"     
        output_file=${snsseq_dir}/${g4_file_short}.${snsseq_file_short}.bed
        echo " --> output_file ${output_file}"
        echo " --> g4_file ${g4_file}"
        # overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        ${bedtools_bin} intersect -wb -a $g4_file -b $snsseq_file > ${output_file}
        
    done
done