# COVID_genotyping_probes

#basic protocol and order of running sh 

Run sh CoreAll.sh

Make sure to delete contents of scratch and plot directories in 
 the Matlab and Multicore directories.
  - Note, the directories need to be present, but the contents should be deleted. 

Run RunAll.m in the Matlab directory
 - you will need to adjust the variables in this file.
    - GenRec.samples = vector of sample directories within the exposure time directories
    - GenRec.is_cDNA = vector of zeros the same dimension as GenRec.samples
    - GenRec.sampletag = vector of the sample names
    - CompSampleNames = Just create a vector of same dimensions as GenRec.samples with the elements 'clade'

Run the RunRscript.sh to make the final calls.
 - Define path to the *012_R_* files in the .sh file.

Optional steps.

Run cat *.fasta > in.fa

Run compare_fasta_R.pl

Currently Under Development

Version log 1.0 basic version of protocol uploaded and shared along with necessary sh and code
Directory information provided 

Verision log 0.1 8/9/2022
Additional scripts and codes related to the SARS-CoV-2 genome tiling array and genotyping probe sets will be deposited and shared
