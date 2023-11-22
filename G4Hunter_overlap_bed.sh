## define the bedtools binary
bedtools_bin=/Users/pmonsieurs/programming/software/bedtools2-master/bin/bedtools

## set intersect between different SNS-seq data of Bridlin
g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter/
snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/ori_predictions


## for Tb427 SNSseq data
g4_files=($(find $g4_dir -name "Tb427*.bed"))
snsseq_files=($(find $snsseq_dir -name "merged*.bed"))

for g4_file in ${g4_files[@]}; do
    # echo $g4_file
    g4_file_short=$(basename ${g4_file})
    for snsseq_file in ${snsseq_files[@]}; do
        # echo $snsseq_file
        snsseq_file_short=$(basename ${snsseq_file})
        overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        echo "${g4_file_short},${snsseq_file_short},${overlap}"
    done
done



## for Tb927 SNSseq data
g4_files=($(find $g4_dir -name "Tb927*.bed"))
snsseq_files=($(find $snsseq_dir -name "*TB927v66.bed"))

for g4_file in ${g4_files[@]}; do
    # echo $g4_file
    g4_file_short=$(basename ${g4_file})
    for snsseq_file in ${snsseq_files[@]}; do
        # echo $snsseq_file
        snsseq_file_short=$(basename ${snsseq_file})
        overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        echo "${g4_file_short},${snsseq_file_short},${overlap}"
    done
done




#### same procedure but now using the extended bed files


## define the bedtools binary
bedtools_bin=/Users/pmonsieurs/programming/software/bedtools2-master/bin/bedtools

## set intersect between different SNS-seq data of Bridlin
g4_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter/
snsseq_dir=/Users/pmonsieurs/programming/trypanosoma_sofi_mining/results/ori/


## for Tb427 SNSseq data
g4_files=($(find $g4_dir -name "Tb427*.bed"))
snsseq_files=($(find $snsseq_dir -name "merged*.bed"))

for g4_file in ${g4_files[@]}; do
    # echo $g4_file
    g4_file_short=$(basename ${g4_file})
    for snsseq_file in ${snsseq_files[@]}; do
        # echo $snsseq_file
        snsseq_file_short=$(basename ${snsseq_file})
        overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        echo "${g4_file_short},${snsseq_file_short},${overlap}"
    done
done



## for Tb927 SNSseq data
g4_files=($(find $g4_dir -name "Tb927*.bed"))
snsseq_files=($(find $snsseq_dir -name "*TB927v66.*.bed"))

for g4_file in ${g4_files[@]}; do
    # echo $g4_file
    g4_file_short=$(basename ${g4_file})
    for snsseq_file in ${snsseq_files[@]}; do
        # echo $snsseq_file
        snsseq_file_short=$(basename ${snsseq_file})
        overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        echo "${g4_file_short},${snsseq_file_short},${overlap}"
    done
done


