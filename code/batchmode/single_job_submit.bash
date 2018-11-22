#!/bin/bash

## export JOB_NAME=SUB_${start_subj}_${end_subj}
#$ -cwd
#$ -q free64
#$ -pe openmp 1
#$ -l mem_free=4G
### #$ -R y

module load MATLAB/r2017b

# ./analyzedata_batch $whsim $dotest $isHPC $start_subj $end_subj 
./analyzedata_batch_v2 $whsim $dotest $isHPC $start_subj $end_subj $logging $logfile

module unload MATLAB/r2017b