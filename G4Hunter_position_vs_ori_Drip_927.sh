## define the bedtools binary
#bedtools_bin=/Users/pmonsieurs/programming/software/bedtools2-master/bin/bedtools
bedtools_bin=/Users/pmonsieurs/programming/software/bedtools2/bin/bedtools

## first run the python script to extend the ORI bed file. For the 
## overlap this was only set to 500. Now, set to 2000bp in line with 
## the publication of Comoglio 2015 Cell Reports - Figure 2A
# /Users/pmonsieurs/programming/leishmania_snsseq/bin/G4Hunter_extend_ori.py


## set intersect between the DRIP-seq data with the ORI sequences
## either true ones or shuffled ones. First need to do conversion from 
## bigwig to bed
bigwig2bed=/Users/pmonsieurs/programming/software/bigWigToBedGraph/bigWigToBedGraph
drip_dir=/Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data/927_DRIP-seq_bw
for bw_file in ${drip_dir}/*.bw; do
    bed_file=${bw_file/.bw/.bed}
    $bigwig2bed $bw_file $bed_file
done
drip_files=($(find $drip_dir -maxdepth 1  -name "*DRIP*.bed"))
snsseq_dir=/Users/pmonsieurs/programming/leishmania_snsseq/results/927/
snsseq_files=($(find $snsseq_dir -name "*ORI*extended*.bed"))


## print the file lists
echo ${drip_files[@]}
echo ${snsseq_files[@]}


## calculate the overlap
for drip_file in ${drip_files[@]}; do
    echo $drip_file
    drip_file_short=$(basename ${drip_file})
    drip_file_short=${drip_file_short/.bed/}
    echo $drip_file_short 
    for snsseq_file in ${snsseq_files[@]}; do
        echo " --> ${snsseq_file}"
        snsseq_file_short=$(basename ${snsseq_file})
        echo " --> basename ${snsseq_file_short}"
        
        ## do substitution depending on the type of file
        # snsseq_file_short=${snsseq_file_short/-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed/}
        # snsseq_file_short=${snsseq_file_short/_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed/}
        snsseq_file_short=${snsseq_file_short/_ORIs_alone_union500_nonoverlap50_woStrand.extended_2000nt.bed/}
        echo " --> snsseq_file_short ${snsseq_file_short}"     
        output_file=${snsseq_dir}/${drip_file_short}.${snsseq_file_short}.bed
        echo " --> output_file ${output_file}"
        echo " --> g4_file ${drip_file}"
        # overlap=$(${bedtools_bin} intersect -a $snsseq_file -b $g4_file 2>/dev/null | wc | awk '{print $1}')
        ${bedtools_bin} intersect -wb -a $drip_file -b $snsseq_file > ${output_file}
        
    done
done