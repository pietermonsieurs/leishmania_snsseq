## define the binary
bigwig2bed=/Users/pmonsieurs/programming/software/bigWigToBedGraph/bigWigToBedGraph

## convert the bw to .bed
cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data/427_Mnase-seq_bw
$bigwig2bed 427_PCF_Amt_WT_rep1_T_brucei_427.bigwig 427_PCF_Amt_WT_rep1_T_brucei_427.bed
$bigwig2bed 427_BSF_YT3_rep2_T_brucei_427.bigwig 427_BSF_YT3_rep2_T_brucei_427.bed