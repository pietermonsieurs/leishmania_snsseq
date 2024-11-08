## define the binary
bigwig2bed=/Users/pmonsieurs/programming/software/bigWigToBedGraph/bigWigToBedGraph

## convert the bw to .bed
cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_427_data/427_Mnase-seq_bw
$bigwig2bed 427_PCF_Amt_WT_rep1_T_brucei_427.bigwig 427_PCF_Amt_WT_rep1_T_brucei_427.bed
$bigwig2bed 427_BSF_YT3_rep2_T_brucei_427.bigwig 427_BSF_YT3_rep2_T_brucei_427.bed

cd /Users/pmonsieurs/programming/leishmania_snsseq/data/for-Pieter_927_data_set3/927_G4_experimental_Marsico/
$bigwig2bed 927_G4_Minus_K.bw 927_G4_Minus_K.bed
$bigwig2bed 927_G4_Plus_K.bw 927_G4_Plus_K.bed
