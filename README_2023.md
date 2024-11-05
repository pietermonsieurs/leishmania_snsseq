# SNS-seq new analysis

## Input data
* Ori prediction to compare with: 
    * prediction by Bridlin on Tb427 2018
    * prediction by Bridlin on Tb927
    * prediction as reported by Marisco et al. 2019
* Run G4 Hunter software in three different modes: 
    * Tb427 reference genome - aim for 5000 sites
    * Tb427 reference genome - aim for 20000 sites
    * Tb927 reference genome - aim for 3500 site
* comparisons to be made
    * G4 Hunter with Tb427 - 5000 sites with results Bridlin on Tb427 
    * G4 Hunter with Tb427 - 20000 sites with results Bridlin on Tb427 
    * G4 Hunter with Tb927 - 5000 sites with results Bridlin on Tb927
    * G4 Hunter with Tb927 - 3500 sites with results Marisco 2019
    * results Bridlin on Tb927 with results Marisco 2019

## G4 Hunter
* different types of running
    * installation *and* running script: [G4Hunter_install.sh](G4Hunter_install.sh)
    * run online on https://bioinformatics.ibp.cz/
        * create / upload your own genome sequences
        * run online with different parameter settings
        * download results from "csv grouped"
* Run G4Hunter using settings that predicted the correct number of sites
    * Tb427 -  5000 sites: size 25 - treshold 1.85 --> 4562 hits // 1.80 --> 6283 hits
    * Tb427 - 20000 sites: size 25 - treshold 1.57 --> 12887 hits // 1.56 --> 30476 hits
    * Tb927 -  3500 sites: size 25 - treshold 2.00 --> 3471 hits
* Convert the csv download file from the online G4 hunter to a bed file: 
    * [G4Hunter_convert_csv2bed.py](G4Hunter_convert_csv2bed.py): initially converting the csv file to a bed file. 
    * G4Hunter data should be split up in positive versus negative strand, so output needs to results in two different bed files: one for the + strand, one for the - strand.


## Overlap G4Hunter versus Ori positions
* ori predictions are based on a very short region. To be able to calculate the overlap with the G4 hunter predictions, we have to extend the ORI predictions
    * [G4Hunter_extend_ori.py](G4Hunter_extend_ori.py)
* calculate the overlap of each of the four G4Hunter scenarios (threshold 1.56, 1.57, 1.8, 1.85) with the merged ori-prediction files merged_BSF and merged_PCF
    * G4Hunter files stored in /Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter: Tb427_window25_score1.56_30476hits.bed	Tb427_window25_score1.85_4562hits.bed Tb427_window25_score1.57_13409hits.bed	Tb427_window25_score1.8_6283hits.bed
    * ori files stored in /Users/pmonsieurs/programming/leishmania_snsseq/results/ori/: merged_BSF-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed and merged_PCF-b_ORIs_alone_union500_nonoverlap50.extended_2000nt.bed
    * overlap can be calculated using the script [G4Hunter_overlap_bed.sh](G4Hunter_overlap_bed.sh). This will give an overlap of the ORI positions and the G4 hunter peaks, but only the number of overlaps, NO information given on the position. To produce information that can be used to make the 'positional plots', see below.  --> SKIP! NOT needed to make the figures
* calculate the relative proportion of the G4 peaks versus the ORI predictions of Bridlin
    * use the bedtools algorithm to get the overlap and return the position of overlap: [G4Hunter_position_vs_ori.sh](G4Hunter_position_vs_ori.sh)
    * do parsing of the output of this overlapping region: [G4Hunter_position_vs_ori_parsing.py](G4Hunter_position_vs_ori_parsing.py)
    * do plotting of the results in a graph that gives the coverage by G4 hunter peaks relatively versus the ORI position (ORI = position 0)
* repeat with randomly shuffled ORI positions: 
    * extend the shuffled / random positions of Bridlin using [G4Hunter_extend_ori.py](G4Hunter_extend_ori.py)
    * use the bedtools algorithm to get the overlap and return the position of overlap: [G4Hunter_position_vs_ori.sh](G4Hunter_position_vs_ori.sh)
    * do the parsing of the random file: [G4Hunter_position_vs_ori_parsing.py](G4Hunter_position_vs_ori_parsing.py)
* visualisation: [G4Hunter_position_vs_ori_plotting.R](G4Hunter_position_vs_ori_plotting.R)
    * initially done on two different folders, one with the "ori" and one with the "shuffled ori". This has been used for the first draft of the paper. 
    * for the revision done on one folder containing both the "oris" and shuffled ones: /Users/pmonsieurs/programming/leishmania_snsseq/results/427_2018/
    * this can be used to create SuppFigure 10 (version1) / SuppFigure 11A (version2)
* polyA/N regions
    * for poly regions length 4 to 8: use inhouse tool
    * for mono regions (nucleotide frequency): use the seqkit tool for the "mono" type
        * run seqkit [polynucleotide_seqkit.sh](polynucleotide_seqkit.sh)
        * plot in R [polynucleotide_seqkit_plot.R](polynucleotide_seqkit_plot.R)


## Tb427 & Tb927 new ref genome
* repeat the analysis but now for the Tb927 reference genome, and create the same plot as before (with polyA + G4Hunter, but now experimental + ORIs as predicted in Tb927) (This is "first plot" in mail Bridlin 09/01/24) --> in the publication this is Figure 3D?
    * polyA regions
        * create the input fasta files, i.e. the fasta file containing sequencing before and after the ORI: [polynucleotide_get_sequences.py](polynucleotide_get_sequences.py). Use the ORIs as stored in 927_ORIs_bed and 927_shuffledORIs-bed
        * calculate the polyA sequences up and downstream of the ORI: [polynucleotide_inhouse.py]
    * calculate the overlap between the G4 data and the ORI, but now using the experimental G4 data instead of the experimentally predicted ones: 
        * first extend the ORIs 2000 up and downstream to be able to calculate the overlap between the ORI and the G4: 
        * use G4 data as stored in 927_G4_experimental_Marsico [G4Hunter_extend_ori.py](G4Hunter_extend_ori.py)
        * calculate overlap between G4 and predicted ORIs: [G4Hunter_position_vs_ori.sh](G4Hunter_position_vs_ori.sh)
        * do the parsing of the files: [G4Hunter_position_vs_ori_parsing.py](G4Hunter_position_vs_ori_parsing.py)
    * Figure 3D: combination of polyA with G4Marisco in 927 background [polynucleotide_G4Hunter_combined_fig3D.R](polynucleotide_G4Hunter_combined_fig3D.R)

* repeat the analysis but now for the Tb427 reference genome instaed of the Tb427_2018 ref genome, with plotting G4 hunter + MNaseSeq data relative to the ORI. (This is "second plot" in mail Bridlin 09/01/24)
    * new reference genome, so polyA tail should be created again
        * create the input fasta files, i.e. the fasta file containing sequencing before and after the ORI: [polynucleotide_get_sequences.py](polynucleotide_get_sequences.py) -> use as input the ORI for which you want to do the plots
            * this is now run for the three directory: for-Pieter_427-2018_data, for-Pieter_427_data, for-Pieter_927_data. In each of those directories are the ORI regions predicted either on the real data, or on the randomly shuffled data
            * output contains fasta sequences 2000nt up and downstream of those ORIs, and can be used to calculate the polyA/T/C/G in the neighborhood of those sequences
            * those fasta sequences are stored in the directory .../results/polynucleotide/
        * calculate the polyA sequences up and downstream of the ORI: [polynucleotide_inhouse.py](polynucleotide_inhouse.py)
    * repeat steps from previous analysis: 
        * extend the ORI: [G4Hunter_extend_ori.py](G4Hunter_extend_ori.py). The resulting bed files is stored in the .../results/427/ directory
        * calculate overlap between G4 and predicted ORIs: [G4Hunter_position_vs_ori_427.sh](G4Hunter_position_vs_ori_427.sh)
        * do the parsing of the files: [G4Hunter_position_vs_ori_parsing.py](G4Hunter_position_vs_ori_parsing.py)
    * include the Mnase-seq data by checking the overlap with the ori positions
        * include the code for this in [G4Hunter_position_vs_ori_427.sh](G4Hunter_position_vs_ori_427.sh)
            * run for normal and shuffled
        * run again the variant python script [G4Hunter_position_vs_ori_parsing_427.py](G4Hunter_position_vs_ori_parsing_427.py) but now using as input file the bed file of the MNase seq data instead of the G4Hunter data
        * visualisation using [polynucleotide_G4Hunter_combined_427.R](polynucleotide_G4Hunter_combined_427.R)



* repeat the analysis but now for the Tb927 reference genome with plotting G4 experimental data from Marisco + DRIP-seq data relative to the ORI. (This is "third A plot" = 3A in mail Bridlin 09/01/24)
    * step from previous step do not to be repeated for Tb927, as this is already done above with the Marisco data, so [G4Hunter_extend_ori.py](G4Hunter_extend_ori.py), [G4Hunter_position_vs_ori.sh] and [G4Hunter_position_vs_ori_parsing.py](G4Hunter_position_vs_ori_parsing.py) are already run
    * do the steps for the DRIP-seq data of Tb927:
        * calculate overlap between the DRIPseq data and the ori's. This is wrongly called G4Hunter however there is nog G4Hunter involved anymore (Experimentally predicted G4 sites). [G4Hunter_position_vs_ori_Drip_927.sh](G4Hunter_position_vs_ori_Drip_927.sh)
        * convert the .bed file from the overlap in previous script to a coverage file (.cov) using [G4Hunter_position_vs_ori_parsing_927.py](G4Hunter_position_vs_ori_parsing_927.py)
        * visualisation using [polynucleotide_G4Hunter_combined_927.R](polynucleotide_G4Hunter_combined_927.R)

* repeat the analysis but now for the Tb427_2018 reference genome with plotting G4 hunter + DRIP-seq data relative to the ORI. (This is "third B plot" = 3B in mail Bridlin 09/01/24)
    * The steps for the overlap between G4Hunter and the ORIs has already been performed as one of the first analyses, and the outputs are stored in: /Users/pmonsieurs/programming/leishmania_snsseq/results/ori and /Users/pmonsieurs/programming/leishmania_snsseq/results/ori_shuffled. So following scripts were used (no need to re-run)
        * [G4Hunter_extend_ori.py](G4Hunter_extend_ori.py)
        * [G4Hunter_position_vs_ori.sh][G4Hunter_position_vs_ori.sh]
        * [G4Hunter_position_vs_ori_parsing.py](G4Hunter_position_vs_ori_parsing.py) 
    * do the steps for the DRIP-seq data of Tb427_2018:
        * calculate overlap between the DRIPseq data and the ori's. This is wrongly called G4Hunter however there is nog G4 involved anymore. [G4Hunter_position_vs_ori_Drip_427_2018.sh](G4Hunter_position_vs_ori_Drip_427_2018.sh)
        * convert the .bed file from the overlap in previous script to a coverage file (.cov) using [G4Hunter_position_vs_ori_parsing_427_2018.py](G4Hunter_position_vs_ori_parsing_927.py)
        * visualisation using [polynucleotide_G4Hunter_combined_427_2018.R](polynucleotide_G4Hunter_combined_427_2018.R)


    
## overview for the different pictures from the publication
* Figure 3D: polynucleotide_G4Hunter_combined_fig3D.R - combination of polyA with G4Marisco in 927 background
* Figure 3E: polynucleotide_G4Hunter_combined_fig3E.R
* Figure 4: DripSeq and MNase-Seq sequencing data: see above
* SuppFigure 8 : polynucleotide_inhouse_plot.R
* SuppFigure 10 / SuppFigure 11A in version2: G4Hunter_position_vs_ori_plotting.R

## overview of the different pictures redone for review
•	Fig3 main manuscript and Supplementary Figure 9 -> poly_nucleotide.poly_all.tiff:  [polynucleotide_inhouse_plot.R](polynucleotide_inhouse_plot.R) (both for mono (min and max = 1) as poly (min = 4 and max = 8))
•	Supplementary Figure 11A -> G4hunter_versus_ori.tiff [G4Hunter_position_vs_ori_plotting.R](G4Hunter_position_vs_ori_plotting.R)
•	Figure 5A -> Tb927_G4_with_polyA.variant.tiff [polynucleotide_G4Hunter_combined_fig3D.R](polynucleotide_G4Hunter_combined_fig3D.R)
•	Figure 5B and Supplementary Figure 11B -> Tb427_window25_score1.85_4562hits_with_polyA.variant.tiff , Tb427_window25_score1.8_6283hits_with_polyA.variant.tiff, Tb427_window25_score1.57_13409hits_with_polyA.variant.tiff, Tb427_window25_score1.56_30476hits_with_polyA.variant.tiff [polynucleotide_G4Hunter_combined_fig3E.R](polynucleotide_G4Hunter_combined_fig3E.R)
•	Figure 5C  still need to think about this one. This would give 4 plots (BSF + BSF, BSF + BSF-PCF, PCF + BSF-PCF and PCF-PCF), correct? 
•	Figure 5D -> Tb427_2018_with_polyA.variant.dual_axes.tiff [polynucleotide_G4Hunter_combined_427_2018.R](polynucleotide_G4Hunter_combined_427_2018.R)



## New analysis - paper revision

### preliminary analyses
* extend the ORIs in all 3 genetic backgrounds with 2000bp, by setting the input directory to a different directory in [G4Hunter_extend_ori.py](G4Hunter_extend_ori.py)


### PolyA polyT (4 nt to 8 nt) centered around the origins - 427-2018 genomic background
* run the python script [polynucleotide_get_sequences.py](polynucleotide_get_sequences.py) to retrieve the sequences around the ORI or around the shuffled ORI for each dataset. 

## update Figure 5D, which used to be Figure 3D
* For the figure 5D, is comparable to the file Tb927_G4_with_polyA.variant.tiff that I sent: done with script [polynucleotide_G4Hunter_combined_fig3D.R](polynucleotide_G4Hunter_combined_fig3D.R). 
    * mail from september 30/09/2024: redo with 927_G4_MinusK.bw and 927_G4_PlusK.bw instead of 927_G4_Minus_PDS.bw and 927_G4_Minus_PDS.bw
    * file merged_BSF-PCF_ORIs_RNASE_927.bed is same file as previously
    * redo G4 extend ori process
        * convert G4 bigwig (bw) to .bed files: [bigwig_bw_to_bed.sh](bigwig_bw_to_bed.sh)
        * file containing ORI's (merged_BSF-PCF_ORIs_RNASE_927.bed) is not changed and is already extended 
        * [G4Hunter_position_vs_ori.sh][G4Hunter_position_vs_ori.sh] (contains also bigwig to bed file commands)
        * [G4Hunter_position_vs_ori_parsing.py](G4Hunter_position_vs_ori_parsing.py)
    * Figure 3D --> now becomes figure 5D: combination of polyA with G4Marisco in 927 background, but now with other G4 predictions: 
        * [polynucleotide_G4Hunter_combined_fig3D.R](polynucleotide_G4Hunter_combined_fig3D.R)
        * make sure to only select those bed files with "_K" in it

* Figure 6C - backgrond 427 - original name was (427_with_polyA.variant.dual_axes.tiff)
    * 

* For Figure 7D, script to use is [polynucleotide_G4Hunter_combined_427_2018.R] but replace only the G4 files with grep on 1.56 to grep on 1.57
    * output file name is Tb427_2018_with_polyA.variant.dual_axes.tiff
    * stored in ~/results/Tb427_2018

