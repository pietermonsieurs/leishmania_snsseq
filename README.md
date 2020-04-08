# Leishmania SNS-seq

## Input data
* reference genome: downloaded from TriTrypDB. Used version indicated with 2018 (release data december 2018). The other file of T. brucei lister is from 2010. https://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427_2018/fasta/data/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta. Total genome size is 42.25Mbp
* Four different input data sets (8 fastq-files) that should be aligned using BWA: [bwa_createcommands.py](bwa_createcommands.py)
* Two different conditions, each with their own control: 
    * Normalize and subtract sample 2 (background control) from sample 1 (SNS non-induced). Call peaks in non-induced sample.
    * Normalize and subtract sample 4 (background control) from sample 3 (SNS induced RNAi Mlp1). Call peaks in induced RNAi Mlp1 mutant sample.
    * Compare the peaks of subtracted 1 and 3 samples (Are the peaks all the same or there are some differences between subtracted 1 and 3 samples).
    * Mapping of different peaks (origins of replication) between subtracted samples 1 and 3.


* do subsampling to bring all samples to the same average read depth. Subsampling can be done using the samtools implementation. Recalculating based on the flagstat statistics mentioned below. 
    * S4 is the smallest sample, everything is calculated relatively versus this one. 
    * different stats can be used (all reads, mapped reads or proper pairs - see below)
    * for now, percentage are calculated based on the  total amount of reads
| Sample | All_reads | Percentage | Mapped    | Percentage | Properly paired | Percentage |
|--------|-----------|------------|-----------|------------|-----------------|------------|
| S1     | 65364709  | 68%        | 45203712  | 62%        | 43630068        | 63%        |
| S2     | 91426202  | 49%        | 64394674  | 44%        | 62964564        | 44%        |
| S3     | 250365416 | 18%        | 139939427 | 20%        | 134851818       | 20%        |
| S4     | 44428519  | 100%       | 28233444  | 100%       | 27518060        | 100%       |

* check the length of the reads, as many very small / narrow peaks have been detected using these input fastq files. 
  * gunzip 1_S1_L001_R1_001.fastq.gz
  * awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' 1_S1_L001_R1_001.fastq
  * visualization of the read lengths: [readlenght_stats.R](readlenght_stats.R)


* do indexing on all the bam-files using samtools
    * ~/programming/software/samtools-1.9/samtools index 1_S1.subsample.bam
    * ~/programming/software/samtools-1.9/samtools index 2_S2.subsample.bam
    * ~/programming/software/samtools-1.9/samtools index 3_S3.subsample.bam
    * ~/programming/software/samtools-1.9/samtools index 4_S4.subsample.bam
* convert bam-file to .bed files using the bedtools samtobed command.  

## Peak detection

### Data preprocessing
Sicer for broad peak detection - no sharp peaks:
* Installation on support:
    * docs available at: https://zanglab.github.io/SICER2/ 
    * installation: `pip3 install --trusted-host pypi.org --trusted-host files.pythonhosted.org SICER2`
    * Installation path: /Library/Frameworks/Python.framework/Versions/3.8/bin/sicer
    * GenomeData.py file that needs to be adopted to work with new species: /lib/python3.8/site-packages/sicer/lib/GenomeData.py
    * GenomeData.py on computer Hideo: /usr/local/lib/python3.6/dist-packages/sicer/lib/GenomeData.py
* debugging: only take the first 10,000 lines to run Sicer. Sicer crashing at certain point.
    * subsetting:
        * head -n 1000000 1_S1.subsample.bed > 1_S1.subsample.head1M.bed
        * head -n 1000000 2_S2.subsample.bed > 2_S2.subsample.head1M.bed
        * head -n 1000000 3_S3.subsample.bed > 3_S3.subsample.head1M.bed
        * head -n 1000000 4_S4.subsample.bed > 4_S4.subsample.head1M.bed
    * run sicer_df: 
        * sicer_df -t 3_S3.subsample.head10k.bed 1_S1.subsample.head10k.bed -c 4_S4.subsample.head10k.bed 2_S2.subsample.head10k.bed -s tbruc  > sicer_df.head10k.out
        * sicer_df -t 3_S3.subsample.head1M.bed 1_S1.subsample.head1M.bed -c 4_S4.subsample.head1M.bed 2_S2.subsample.head1M.bed -s tbruc  > sicer_df.head1M.out

### SICER2
* Sicer can only work with genomes which are integrated into the system. It needs to know all the chromosome IDs and the corresponding lengths.
    * For T. brucei, only a contig assembly is available. Adding all these contigs done in an automized way: [sicer_addorganism.py](sicer_addorganism.py)
    * add output of previous script to corresponding configuration file of sicer2 --> /Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/sicer/lib/GenomeData.py
* run Sicer with differential peak detection, including controls:
    * output file does not have to be specified. Sicer create its own names for output. Output directory might be beneficial to keep overview
    * sicer_df -t 3_S3.subsample.bed 1_S1.subsample.bed -c 4_S4.subsample.bed 2_S2.subsample.bed -s tbruc -o ./sicer/
    * Sicer2 tool gives always errors when running. Not clear where this comes from.
* Use older version of SICER (version 1.1)
    * installed on the linux computer of Hideo!!
    * Download link found at https://home.gwu.edu/~wpeng/Software.htm
    * download of software via: https://home.gwu.edu/~wpeng/SICER_V1.1.tgz
    * This version is not completely python based, but rather a bash script calling different python modules
        * header of the bash-script should be changed, to give the correct path to the sicer directory: SICER=/home/smrtanalysis/snsseq/software/SICER_V1.1/SICER/ in SICER-df.sh and PATHTO=/home/smrtanalysis/snsseq/software/SICER_V1.1/ in SICER.sh in directory /home/smrtanalysis/snsseq/software/SICER_V1.1/SICER
        * specis hard coded in the SICER-df script --> set to "tbruc". Also, the GenomeData.py script should be updated in the lib directory of SICER, similar as what needed to be done when using SICER2
        * also the python scripts are stored under the lib-directory of the sicer directory. However this runs python2.7 so should be run within a virtual environment
            * virtualenv -p /usr/bin/python2.7 --distribute temp-python
            * source /home/smrtanalysis/snsseq/software/SICER_V1.1/SICER/lib/temp-python/bin/activate
            * pip install numpy
            * pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org scipy
        * after activation of the virtual environment, run SICER bash script
    * Running SICER:
        * /home/smrtanalysis/snsseq/software/SICER_V1.1/SICER/SICER-df.sh  3_S3.subsample.bed 4_S4.subsample.bed 1_S1.subsample.bed 2_S2.subsample.bed 200 600 0.05 0.05
        * End results is a list of differential peaks between both background corrected samples
        * important! only differential peaks when they were identified in peaks when comparing with their background. No peak detected = not included in the differential analysis!
    
### MACS
The final aim is to find peaks in the third condition (S3) that cannot be found back in the first one (S1), and this after background correction.

* step 1: peak detection using MACS
    * installation of MACS via conda on laptop: `conda install -c bioconda macs2`. After installation, MACS can be run via macs2
    * MACS integrates background correction, so can be used to control S1 with S2 output, and S3 with S4.
    * can be run in two different modes: narrowPeaks (default) and broadPeak (set via de parameter --broad)
    * also effective genome size has to be set. For human genome, this is set to 90% of the full genome size (3Gbp --> 90% = 2.7Gbp = 2.7e9). For T. brucei lister, the genome size = 42.25 Mbp --> 90% = 38Mbp = 38e6
    * different runs:
        * default: 
            * S1 vs S2: `macs2 callpeak -f BAM -g 38e6 -n S1 -t ../bwa/1_S1.subsample.bam -c ../bwa/2_S2.subsample.bam`
            * S3 vs S4: `macs2 callpeak -f BAM -g 38e6 -n S3 -t ../bwa/3_S3.subsample.bam -c ../bwa/4_S4.subsample.bam`
        * broad peaks
            * S1 vs S2: `macs2 callpeak -f BAM -g 38e6 -n S1_broad --broad -t ../bwa/1_S1.subsample.bam -c ../bwa/2_S2.subsample.bam`
            * S3 vs S4: `macs2 callpeak -f BAM -g 38e6 -n S3_broad --broad -t ../bwa/3_S3.subsample.bam -c ../bwa/4_S4.subsample.bam`
* step 2: overlap between both files: select those peaks present in S3 and not in S1
    * use the bedtools intersect method to calculate intersect between two experiments. The amount of overlap is defined as the fraction of the whole genome. http://quinlanlab.org/tutorials/bedtools/bedtools.html#bedtools-intersect --> seems to be quite robust: varying between 0.001 and 1e-9 has no, or only limited effect
    * when comparing the tool `intersect` with `subtract`, subtract seems to be less stringent i.e. if you subtract B from A, almost all peaks of A are retained, which is not the case for intersect option.
    * detect peaks present in only one condition: 
        * only in S3: `~/programming/software/bedtools2-master/bin/bedtools intersect -v -f 1e-6 -F 1e-6 -b S1_peaks.narrowPeak -a S3_peaks.narrowPeak > only_present_in_S3.csv`
        * only in S1: `~/programming/software/bedtools2-master/bin/bedtools intersect -v -f 1e-6 -F 1e-6 -a S1_peaks.narrowPeak -b S3_peaks.narrowPeak > only_present_in_S1.csv`
        * broad - only in S3: `~/programming/software/bedtools2-master/bin/bedtools intersect -v -f 1e-6 -F 1e-6 -b S1_broad_peaks.broadPeak -a S3_broad_peaks.broadPeak > only_present_in_S3_broad.csv`
        * broad - only in S1: `~/programming/software/bedtools2-master/bin/bedtools intersect -v -f 1e-6 -F 1e-6 -a S1_broad_peaks.broadPeak -b S3_broad_peaks.broadPeak > only_present_in_S1_broad.csv`
* step 3: visualization of the peaks to detect aberrations from what is expected. Using matplotlib in combnation with samtools depth to make coverage depth plots. 
    * make sure to run samtools depth with -a option. Otherwise position with zero depth will not be reported and all data frames will not have the same length.
    * [differential_peaks_zoom.py](differential_peaks_zoom.py): creates per peak a plot with the sequencing coverage for a peak in either of the 4 samples (S1 to S4):
        * Left: sequencing depth for the peak + 200 bp flanking windows on both sides
        * Middle: coverage of S1 and S3 when subtracted from their corresponding control (S2 or S4 respectively). This might result in negative values, in case the control has a higher coverage for a position in the genome than the actual sample.
        * Same as the middle figure, but negative values have been corrected to 0. 
* step 4: convert MACS output files to bigwig files
    * cp S1_peaks.csv S1_peaks.copy
    * sed '/^#/ d' S1_peaks.csv > S1_peaks.wig
    * awk '{print $1"\t"$2"\t"$3"\t"$4}' S1_peaks.wig > S1_peak.bigwig
    


### SWEMBL
export PATH=/Users/pmonsieurs/programming/software/SWEMBL/:/Users/pmonsieurs/programming/software/samtools-1.9:$PATH

SWEMBL -m 500 -f 150 -R 0.0025 -F -i 1_S1.subsample.bam -a 2_S2.subsample.bam


### in-house developed tool
Start from the samtools depth output generated from [samtools_depth.sh](samtools_depth.sh) and use the core of the [chromosome_peaks.py](chromosome_peaks.py) script to calculate the difference. On this difference, and easy-to-use peak detection algorithm can be applied. 
    

#




https://github.com/iakerman/SNS-seq/blob/master/github_SNS-seq_pipeline_6-2019.txt


# Use this to create your temporary python "install"
# (Assuming that is the correct path to the python interpreter you want to use.)
virtualenv -p /usr/bin/python2.7 --distribute temp-python

# Type this command when you want to use your temporary python.
# While you are using your temporary python you will also have access to a temporary pip,
# which will keep all packages installed with it separate from your main python install.
# A shorter version of this command would be ". temp-python/bin/activate"
source /home/smrtanalysis/snsseq/software/SICER_V1.1/SICER/lib/temp-python/bin/activate
pip install numpy
pip install scipy
pip install --trusted-host=pypi.python.org --trusted-host=pypi.org --trusted-host=files.pythonhosted.org scipy

# When you no longer wish to use you temporary python type
deactivate
