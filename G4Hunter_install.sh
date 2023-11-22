## export pythonpath for python2.7
export PYTHONPATH=/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python2.7/site-packages/:${PYTHONPATH}

## python install libraries 
pip2.7 install --prefix="/user/antwerpen/205/vsc20587/data/software/python_lib/" matplotlib
pip2.7 install --prefix="/user/antwerpen/205/vsc20587/data/software/python_lib/" numpy
pip2.7 install --prefix="/user/antwerpen/205/vsc20587/data/software/python_lib/" biopython

## git clone G4Hunter package
cd programming/software
git clone https://github.com/AnimaTardeb/G4-hunter.git

## do a test run
python2.7 ./G4Hunter.py  -i ~/programming/leishmania_susl/data/refgenomes/TriTrypDB-46_LinfantumJPCM5_Genome.fasta -o /Users/pmonsieurs/programming/leishmania_snsseq/results/g4hunter -w 50 -s 1

## get TriTrypDB genomes. Two different versions should 
## be downloaded: Tb427 is used by Bridlin, Tb927 is used 
## in some other publications.
cd /Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome
wget https://tritrypdb.org/common/downloads/release-55/TbruceiLister427_2018/fasta/data/TriTrypDB-55_TbruceiLister427_2018_Genome.fasta

## run G4 Hunter on the TriTrypDB genome
python2.7 /Users/pmonsieurs/programming/software/G4-hunter/G4Hunter.py -i /Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-55_TbruceiLister427_2018_Genome.fasta -o /Users/pmonsieurs/programming/leishmania_snsseq/results/ -w 50 -s 2


python2.7 /Users/pmonsieurs/programming/software/G4-hunter/G4Hunter.py -i /Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-55_TbruceiLister427_2018_Genome.fasta -o /Users/pmonsieurs/programming/leishmania_snsseq/results/ -w 50 -s 2



## run with different parameters to produce different kind of data sets

## set 1: return ~ 20.000 G4 quadruplexed for Tb427
python2.7 /Users/pmonsieurs/programming/software/G4-hunter/G4Hunter.py -i /Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-55_TbruceiLister427_2018_Genome.fasta -o /Users/pmonsieurs/programming/leishmania_snsseq/results/ -w 25 -s 1.5

## set 1: return ~ 5.000 G4 quadruplexed for Tb427
python2.7 /Users/pmonsieurs/programming/software/G4-hunter/G4Hunter.py -i /Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-55_TbruceiLister427_2018_Genome.fasta -o /Users/pmonsieurs/programming/leishmania_snsseq/results/ -w 25 -s 1.5


python2.7 /Users/pmonsieurs/programming/software/G4-hunter/G4Hunter.py -i /Users/pmonsieurs/programming/leishmania_snsseq/data/refgenome/TriTrypDB-63_TbruceiTREU927_Genome.fasta -o /Users/pmonsieurs/programming/leishmania_snsseq/results/ -w 25 -s 2

