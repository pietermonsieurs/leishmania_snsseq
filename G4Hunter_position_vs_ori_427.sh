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
g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data/427_G4hunter_predictions_bed/
snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/427/
g4_files=($(find $g4_dir -maxdepth 1  -name "427*.bed"))
snsseq_files=($(find $snsseq_dir -name "427_merged_*.extended_2000nt.bed"))

## set intersect between different shuffled SNS-seq data of Bridlin, but now
## with the shuffled ORI sequences
g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data/427_G4hunter_predictions_bed/
snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/427/
g4_files=($(find $g4_dir -maxdepth 1  -name "427*.bed"))
snsseq_files=($(find $snsseq_dir -name "427_shuffeled*.extended_2000nt.bed"))

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
        snsseq_file_short=${snsseq_file_short/_ORIs_alone_union500_nonoverlap50_woStrand.extended_2000nt.bed/}
        echo " --> snsseq_file_short ${snsseq_file_short}"     
        output_file=${snsseq_dir}/${g4_file_short}.${snsseq_file_short}.bed
        echo " --> output_file ${output_file}"
        echo " --> g4_file ${g4_file}"
        # overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        ${bedtools_bin} intersect -wb -a $g4_file -b $snsseq_file > ${output_file}
        
    done
done




## get the MNase-seq data 
mnase_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data/427_Mnase-seq_bw/
mnase_files=($(find $mnase_dir -maxdepth 1  -name "427*.bed"))
echo ${mnase_files[@]}

## run for shuffled sequences as well as the normal ori sequences!!
snsseq_files=($(find $snsseq_dir -name "427_merged_*.extended_2000nt.bed"))
snsseq_files=($(find $snsseq_dir -name "427_shuffeled*.extended_2000nt.bed"))

## afterwards remove the file where BSF/PCF in sns and BSF/PCF 
## are mixed up in the same file. Only do comparison per sample so e.g. 
## remove 427_BSF_YT3_rep2_T_brucei_427.427_shuffeled_seed668_PCF.bed

## calculate the overlap 
for mnase_file in ${mnase_files[@]}; do
    echo $mnase_file
    mnase_file_short=$(basename ${mnase_file})
    mnase_file_short=${mnase_file_short/.bed/}
    echo $mnase_file_short 
    for snsseq_file in ${snsseq_files[@]}; do
        echo " --> ${snsseq_file}"
        snsseq_file_short=$(basename ${snsseq_file})
        echo " --> basename ${snsseq_file_short}"
        
        ## do substitution depending on the type of file
        # snsseq_file_short=${snsseq_file_short/-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed/}
        # snsseq_file_short=${snsseq_file_short/_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed/}
        snsseq_file_short=${snsseq_file_short/_ORIs_alone_union500_nonoverlap50_woStrand.extended_2000nt.bed/}
        echo " --> snsseq_file_short ${snsseq_file_short}"     
        output_file=${snsseq_dir}/${mnase_file_short}.${snsseq_file_short}.bed
        echo " --> output_file ${output_file}"
        echo " --> g4_file ${mnase_file}"
        # overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        ${bedtools_bin} intersect -wb -a $mnase_file -b $snsseq_file > ${output_file}
        
    done
done