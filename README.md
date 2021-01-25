# Leishmania SNS-seq

## Input data & alignment
* Two different conditions, each with their own control: 
    * Normalize and subtract sample 2 (background control) from sample 1 (SNS non-induced). Call peaks in non-induced sample.
    * Normalize and subtract sample 4 (background control) from sample 3 (SNS induced RNAi Mlp1). Call peaks in induced RNAi Mlp1 mutant sample.
    * Compare the peaks of subtracted 1 and 3 samples (Are the peaks all the same or there are some differences between subtracted 1 and 3 samples).
    * Mapping of different peaks (origins of replication) between subtracted samples 1 and 3.
* reference genome: downloaded from TriTrypDB. Used version indicated with 2018 (release data december 2018). The other file of T. brucei lister is from 2010. https://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427_2018/fasta/data/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta. Total genome size is 42.25Mbp

## Alignment
* After a first analyis, it was noticed that the fastq file contained very short reads, based on the trimming done by the sequencing company ==> check the length of the reads, as many very small / narrow peaks have been detected using these input fastq files. 
  * gunzip 1_S1_L001_R1_001.fastq.gz
  * awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' 1_S1_L001_R1_001.fastq
  * visualization of the read lengths: [readlenght_stats.R](readlenght_stats.R)
  * Reads are now filtered on length using the script: [reads_filtering_on_length.sh](reads_filtering_on_length.sh), which is integrated in the BWA bash script (see further below)
* Run alignment, which includeds also the filtering of the input file bases on the read length. Overall bash script is [bwa_snsseq.sh](bwa_snsseq.sh)
    * Four different input data sets (8 fastq-files) that should be aligned using BWA: [bwa_createcommands.py](bwa_createcommands.py). Those commands are inserted in the bash script mentioned above. 
        * many short peaks are detected with only a length of 19nt. This corresponds perfectly with the default seed lengths of the BWA algorithm. No idea why, but apparently 19nt seed is enough to be reported. Therefore, minimum seed length is now set to 100: redo the analysis and set the minimum seed length to 100 (-k 100) when running BWA.
        * alternative could be to use Bowtie. However, the publication of Lombrana (2016) is using bowtie version 1 instead of bowtie2, and the parameter settings mentioned there are not available anymore in bowtie2
    * remove duplicate alignments from the bam-file: included the Picard MarkDuplicates command into the bwa bash script, with the option activated to directly remove the duplicates. Between 13% and 35% of the alignments are removed from the bam-files.

* do subsampling to bring all samples to the same average read depth. Subsampling can be done using the samtools implementation. Recalculating based on the flagstat statistics mentioned below. 
    * option 1:
        * S4 is the smallest sample, everything is calculated relatively versus this one. 
        * different stats can be used (all reads, mapped reads or proper pairs - see below)
        * for now, percentage are calculated based on the  total amount of reads
            | Sample | All_reads | Percentage | Mapped    | Percentage | Properly paired | Percentage |
            |--------|-----------|------------|-----------|------------|-----------------|------------|
            | S1     | 65364709  | 68%        | 45203712  | 62%        | 43630068        | 63%        |
            | S2     | 91426202  | 49%        | 64394674  | 44%        | 62964564        | 44%        |
            | S3     | 250365416 | 18%        | 139939427 | 20%        | 134851818       | 20%        |
            | S4     | 44428519  | 100%       | 28233444  | 100%       | 27518060        | 100%       |
    * option 2:
        * first select properly paired reads, and only continue with them:
            * samtools view -bf 0x2 1_S1.bam > 1_S1.proper_paired.bam 
            * samtools view -bf 0x2 2_S2.bam > 2_S2.proper_paired.bam 
            * samtools view -bf 0x2 3_S3.bam > 3_S3.proper_paired.bam
            * samtools view -bf 0x2 4_S4.bam > 4_S4.proper_paired.bam 
        * check flagstat again for all four bam-files, and do recalculation for doing subsampling. Now, the bam files are only based on fastq-file with reads with sufficient length (length > 100nt) and where only proper paired reads are selected. This should further refine the peaks
        * overview tables
            * without increasing the seed length
                | Sample | All\_reads | Percentage |
                |--------|------------|------------|
                | S1     | 20163824   | 83%        |
                | S2     | 19440774   | 86%        |
                | S3     | 16652744   | 100%       |
                | S4     | 18672428   | 89%        |
            * with increasing the seed length to 100 (-k 100)
                | Sample | All\_reads | Percentage |
                |--------|------------|------------|
                | S1     | 15434734   | 80%        |
                | S2     | 15430228   | 80%        |
                | S3     | 12395674   | 100%       |
                | S4     | 14800620   | 84%        |
    * option 3: 
        * duplicate reads have been remove. Check flagstats again to see the percentage of proper paired reads. 
        * select proper paired reads: 
            * samtools view -bf 0x2 1_S1.removedups.bam > 1_S1.removedups.proper_paired.bam 
            * samtools view -bf 0x2 2_S2.removedups.bam > 2_S2.removedups.proper_paired.bam 
            * samtools view -bf 0x2 3_S3.removedups.bam > 3_S3.removedups.proper_paired.bam
            * samtools view -bf 0x2 4_S4.removedups.bam > 4_S4.removedups.proper_paired.bam 

    * option 4: Most stringent filtering. Redo the analysis with different filtering steps implemented, most if them are implemented in the bash script [bwa_snsseq.sh](bwa_snsseq.sh)
        * check on length of the reads (see also above: [reads_filtering_on_length.sh](reads_filtering_on_length.sh))
        * BWA: with seed parameter of 100 = very stringent
        * samtools view: remove the alignments with a mapping quality lower than 30
        * Picard: remove duplicates with the option REMOVE_DUPLICATES=TRUE
        * samtools view: select only proper paired reads. 
        * subsampling now integrated in the bwa_snsseq.sh script
    * running subsampling: [subsampling.sh](subsampling.sh) to subsample all bam-file to the same amount of reads. Output looks like: 4_S4.proper_paired.subsample.bam --> ! Update: this step is now also integrated in [bwa_snsseq.sh](bwa_snsseq.sh)

* do indexing on all the bam-files using samtools
* convert bam-file to .bed files using the bedtools samtobed command. 
    * [bedtools_bamtobed.sh](bedtools_bamtobed.sh) !! this step is *not* needed if you want to run MACS, only if you want to run SICER2.
    * afterwards copy-paste bed file to laptop for peak detection. !! Update: Sicer2 now also installed on the calcua cluster, see further below. 


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
Summary: SICER2 was initially not running on MacBook, which forced us to use SICER1.1 After mailing with authors, this seems to give the same output as SICER 2. The authors updated the SICER package, and now SICER 2 is running on MacBook (March2020). Later, the SICER2 pacakges was also installed locally on the CalcUA cluster (see further below), so it can be run on multiple clusters. Initially (when running on MacBook), the *sicer_df* package was used where in one command all differential peaks were predicted (i.e. integrating information from S1 until S4). However, in a second step, we want to use two or more prediction tools, and look for overlapping peaks either detected in S1 or detected in S3. Therefore, we will run separatley for the two conditions using the *sicer* package. 

* Sicer can only work with genomes which are integrated into the system. It needs to know all the chromosome IDs and the corresponding lengths.
    * For T. brucei, only a contig assembly is available. Adding all these contigs done in an automized way: [sicer_addorganism.py](sicer_addorganism.py)
    * add output of previous script to corresponding configuration file of sicer2 --> /Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/sicer/lib/GenomeData.py
* run Sicer with differential peak detection, including controls:
    * output file does not have to be specified. Sicer create its own names for output. Output directory might be beneficial to keep overview
    * sicer_df -t 3_S3.proper_paired.subsample.bed 1_S1.proper_paired.subsample.bed -c 4_S4.proper_paired.subsample.bed 2_S2.proper_paired.subsample.bed -s tbruc -o ./sicer/
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

* Run sicer (not sicer_df) on S1 and S3 respectively to identify peaks. This can be run on CalcUA now. 
    * Run SICERInstallation of SICER2 on CalcUA: see mail of Franky Backeljauw 20200304:
        * install in python-lib of software directory: /user/antwerpen/205/vsc20587/software/python_lib/lib/python3.7/site-packages
        * export PYTHONPATH=/user/antwerpen/205/vsc20587/software/python_lib/lib/python3.7/site-packages/:$PYTHONPATH
        * pip3 install --prefix="/user/antwerpen/205/vsc20587/software/python_lib/" sicer2
        * pip3 show sicer2
        * binaries are installed in /user/antwerpen/205/vsc20587/software/python_lib/bin
        * copy-paste GenomeData.py (see above) from MacBook to /user/antwerpen/205/vsc20587/software/python_lib/lib/python3.7/site-packages/sicer/lib: this contains genome data for T. brucei
    * Running sicer2 needs different input parameters: 
        * the parameters set in the previous experiment with sicer_df is probably the same as now. But default settings look ok: --window_size 200 and --gap_size 600. False discovery rate is by default 0.01 (--false_discovery_rate) which seems more stringent than the 0.05 in previous setting
        * run SICER for S1 and S3 with own control
            * /user/antwerpen/205/vsc20587/software/python_lib/bin/sicer -t 1_S1.mapq30.removedups.proper_paired.subsample.bed -c 2_S2.mapq30.removedups.proper_paired.subsample.bed -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/sicer/ -s tbruc --cpu 20
            * /user/antwerpen/205/vsc20587/software/python_lib/bin/sicer -t 3_S3.mapq30.removedups.proper_paired.subsample.bed -c 4_S4.mapq30.removedups.proper_paired.subsample.bed -o /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/results/sicer/ -s tbruc --cpu 20




    
### MACS
The final aim is to find peaks in the third condition (S3) that cannot be found back in the first one (S1), and this after background correction.

#### MACS on calcua
* installation of MACS2 on CalcUA server. Use same approach as done for SICER2
    * module load Python/3
    * export PYTHONPATH=/user/antwerpen/205/vsc20587/software/python_lib/lib/python3.7/site-packages/:$PYTHONPATH
    * pip install --prefix="/user/antwerpen/205/vsc20587/software/python_lib/" MACS2
    * MACS2 is installed under: /user/antwerpen/205/vsc20587/software/python_lib/bin/ --> when running, first do PYTHONPATH and do load Python3 module
* running peak detection
    * macs2 is installed on calcua with approach above in `/user/antwerpen/205/vsc20587/software/python_lib/bin`
    * add new parameters: --no-model
    * callpeak: 
        * /user/antwerpen/205/vsc20587/software/python_lib/bin/macs2 callpeak -f BAM -g 38e6 --nomodel --extsize 150 -n S1 -t ~/scratch/leishmania_snsseq/results/bwa/1_S1.mapq30.removedups.proper_paired.subsample.bam -c ~/scratch/leishmania_snsseq/results/bwa/2_S2.mapq30.removedups.proper_paired.subsample.bam
        * /user/antwerpen/205/vsc20587/software/python_lib/bin/macs2 callpeak -f BAM -g 38e6 --nomodel --extsize 150  -n S3 -t ~/scratch/leishmania_snsseq/results/bwa/3_S3.mapq30.removedups.proper_paired.subsample.bam -c ~/scratch/leishmania_snsseq/results/bwa/4_S4.mapq30.removedups.proper_paired.subsample.bam


#### approach only focussing on MACS, not using the SICER output
* step 1: peak detection using MACS
    * installation of MACS via conda on laptop: `conda install -c bioconda macs2`. After installation, MACS can be run via macs2
    * MACS integrates background correction, so can be used to control S1 with S2 output, and S3 with S4.
    * can be run in two different modes: narrowPeaks (default) and broadPeak (set via de parameter --broad)
    * also effective genome size has to be set. For human genome, this is set to 90% of the full genome size (3Gbp --> 90% = 2.7Gbp = 2.7e9). For T. brucei lister, the genome size = 42.25 Mbp --> 90% = 38Mbp = 38e6
    * different runs:
        * default: 
            * S1 vs S2: `macs2 callpeak -f BAM -g 38e6 -n S1 -t ../bwa/1_S1.subsample.bam -c ../bwa/2_S2.subsample.bam`
            * S3 vs S4: `macs2 callpeak -f BAM -g 38e6 -n S3 -t ../bwa/3_S3.subsample.bam -c ../bwa/4_S4.subsample.bam`
        * update after filtering: 
            * S1 vs S2: `macs2 callpeak -f BAM -g 38e6 -n S1 -t ../bwa/1_S1.proper_paired.subsample.bam -c ../bwa/2_S2.proper_paired.subsample.bam`
            * S3 vs S4: `macs2 callpeak -f BAM -g 38e6 -n S3 -t ../bwa/3_S3.proper_paired.subsample.bam -c ../bwa/4_S4.proper_paired.subsample.bam`
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
    * make sure to run samtools depth with -a option. Otherwise position with zero depth will not be reported and all data frames will not have the same length. Script has to be run on the supercomputer for each contig separately. 
        * depth directory contains list of all chromosomes / contigs 
        * bash script: [samtools_depth.sh](samtools_depth.sh)
    * perform filtering on the peaks that are detected, for example by only selecting the one with the highest coverage (example for S3):
        * `mv only_present_in_S3.csv only_present_in_S3.all.csv`
        * `awk '$5 > 100' only_present_in_S3.all.csv > only_present_in_S3.csv`

    * [differential_peaks_zoom.py](differential_peaks_zoom.py): creates per peak a plot with the sequencing coverage for a peak in either of the 4 samples (S1 to S4):
        * Left: sequencing depth for the peak + 200 bp flanking windows on both sides
        * Middle: coverage of S1 and S3 when subtracted from their corresponding control (S2 or S4 respectively). This might result in negative values, in case the control has a higher coverage for a position in the genome than the actual sample.
        * Same as the middle figure, but negative values have been corrected to 0. 
* step 4: convert MACS output files to bigwig files
    * cp S1_peaks.csv S1_peaks.copy
    * sed '/^#/ d' S1_peaks.csv > S1_peaks.wig
    * awk '{print $1"\t"$2"\t"$3"\t"$4}' S1_peaks.wig > S1_peak.bigwig ! wrong!! need specialized script


#### updated MACS2 with switched samples [July 2020]
Samples might have been switched. S1 and S2 might be switched, as well as S3 and S4. This has been changed in the bwa_snsseq.sh script and renamed to [bwa_snsseq_sampleswitch.sh](bwa_snsseq_sampleswitch.sh)
* S1 is now control and S2 is the treatment sample. Same for S3 (control) and S4 (treatment)
* also adapt the visualization script []   




### SWEMBL
export PATH=/Users/pmonsieurs/programming/software/SWEMBL/:/Users/pmonsieurs/programming/software/samtools-1.9:$PATH

SWEMBL -m 500 -f 150 -R 0.0025 -F -i 1_S1.proper_paired.subsample.bam -a 2_S2.proper_paired.subsample.bam

Update: swembl is giving a segmentation fault. Might be related to the bam-file which is corrupt. The segementation fault is not there when running a bam-file from the archive (~/programming/leishmania_snsseq/results/bwa/archive/2_S2.subsample.bam)

Docs:
http://www.bioconductor.org/help/course-materials/2010/EMBL2010/Chip-seq.SWEMBL.pdf




## Combining MACS and SICER2 [May 2020]
This part was used to get the most stringent conditions. Those conditions seem to be too stringent, so now switched back (temporary?) to only use MACS. 

### Data Integration [May 2020]
Integrating the peaks predicted by MACS and SICER2 (and in the future other tools like SWEMBL?), in order to only retain those peaks that are predicted with both methodologies, in order to have the highest confidence.
* use the bedtools intersect command to look for overlap between both approaches. For MACS, we have to choose between narrow and broad peaks. Broad peaks are still relatively narrow compared to sicer, but return more peaks than narrow peak (= default MACS2 setting)
    * first convert the peak / summit files of SICER and MACS into .bed-files of *only* three columns
    * next run bedtools using the bedtools intersect -wa setting
    * the default overlap fraction is 1e-9, which corresponds to 1nt with human genome. 
    * see detailed commands in [bwa_snsseq.sh](bwa_snsseq.sh)
* check the unique peaks either in S1 or S3 by using the -v option with bedtools intersect
    * commands have been added to [bwa_snsseq.sh](bwa_snsseq.sh)

### Stringent approach using MACS and SICER2 --- kind of duplicate of data integration ----  [May 2020]
Both algorithms have been run separately. Afterward, the peaks from MACS and SICER2 are combined, and only the peaks confirmed in both algorithms, are reported. Next, only the peaks uniquely predicted in S1 or uniquely predicted in S3 are reported.
* all the different steps have been combined into on bash sript: [bwa_snsseq.sh](bwa_snsseq.sh)
    * first different bwa alignment steps and filtering steps are performed
    * next overlap between both algorithms is done --> only when confirmation using both approaches, the peaks are retained (separate procedure for narrow and broad peaks)
    * select the unique peaks for either S1 or S3







## full pictures of sequencing data
* request to visualize full sequencing data set. First run depth for all 4 samples, then subtract with the corresponding background, and do the visualization of S1 on top (bg-corrected with S2), and S3 below (bg-corrected with S4): [chromosome_peaks.py](chromosome_peaks.py). Additionally, include the raw data on top of the plot for each of the 4 samples. 
    * get the depth for all chromosomes and all samples (n=4)
        * using samtools depth to get the depth for all chromosomes and all samples
        * re-run only needed when something is changed on the .bam files (e.g. additional filterin step). The name of the corresponding bam-file needs to be adapted in the script. 
    * create picture with sequencing depth for each of the differential peaks 
        * software code [chromosome_peaks.py](chromosome_peaks.py) was now used for creating pictures of 50kb spanning the whole chromosome. This has to be adpated that 1) it works with a specified chromosome + start and end position, 2) that is does not work on all chromosomes combined, but takes one parameter setting by one
        * new script is called: [differential_peaks.py](differential_peaks.py): takes as input three arguments: chrom, start, end
        



### in-house developed tool
Start from the samtools depth output generated from [samtools_depth.sh](samtools_depth.sh) and use the core of the [chromosome_peaks.py](chromosome_peaks.py) script to calculate the difference. On this difference, and easy-to-use peak detection algorithm can be applied. 




    

## New data

### Input data
New input data have been generated. More explanation on the biological background of the different samples can be found in the mail of Slavica on Jan 4th 17:09 (docx file also stored in the "docs" directory of the project directory). The first sample is always the real sample, the second sample is the control (e.g. 1 = sample, 2=background). In summary:
* samples 9 until 14: SNS-seq for Procylic forms, i.e. three replicates each with their background
* samples 15 until 20: SNS-seq for bloodstream forms, i.e. three replicates each with their background
* samples 1 until 8: replication of previous experiment of Akila, with Mlp1 gene knock-out (with background + control without knock-out and it corresponding background)
All those data have been stored in a new subdirectory: /user/antwerpen/205/vsc20587/scratch/leishmania_snsseq/data/20210104

### BWA alignment
* BWA alignment is run via bash script [bwa_snsseq_newdata.sh](bwa_snsseq_newdata.sh), which is largely based on the general BWA protocol
* creating the different qsub commands using the -v option can be done using the script [bwa_snsseq_newdata_create_qsub_commands.py](bwa_snsseq_newdata_create_qsub_commands.py)
* subsample the different paired control - sample bam files to the lowest number of reads
    * for each pair of samples, the first one (e.g. 1) is the sample, the second one (e.g. 2) is the control or background
    * calculate for each pair which one is the largest, and calculate ratio for subsampling: [subsample_bamfiles.py](subsample_bamfiles.py)
    * output of this script are direct linux commands (inlcuding samtools commands) which can be used as basis for the qsub script
        * add qsub header
        * load the SAMtools
        * for now stored in [subsample_bamfiles.sh](subsample_bamfiles.sh)
* data sent via belnet filesender








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
