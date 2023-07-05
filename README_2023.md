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
* Convert the csv download file from the online G4 hunter to a bed file
    * 