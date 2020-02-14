# Leishmania SNS-seq

## Input data
* reference genome: downloaded from TriTrypDB. Used version indicated with 2018 (release data december 2018). The other file of T. brucei lister is from 2010. https://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427_2018/fasta/data/TriTrypDB-46_TbruceiLister427_2018_Genome.fasta
* Four different input data sets (8 fastq-files) that should be aligned using BWA: [bwa_createcommands.py](bwa_createcommands.py)

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
* 

