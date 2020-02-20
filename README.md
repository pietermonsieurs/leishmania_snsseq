# Leishmania SNS-seq

## Input data
* reference genome: downloaded from TriTrypDB. Used version indicated with 2018 (release data december 2018). The other file of T. brucei lister is from 2010. https://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427_2018/fasta/data/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta
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



* convert bam-file to .bed files using the bedtools samtobed command.  

## Peak detection

### SICER2
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
