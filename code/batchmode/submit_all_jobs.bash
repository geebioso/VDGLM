#!/bin/bash

##############################################################################

# This file will iterate through sequences of subjects and submit a 
# job to the scheduler taking start subject, end subject, simulation number,
# and whether we are testing as inputs. The script dynamically assigns the
# job name, error log file, and output log file 

# INPUTS: 
# 	$1: the simulation number 
# 	$2: run in test mode? (0/1) 
#  	$3: are we running on the HPC? (0/1)
# 	$4: starting subject
# 	$5: subject interval 
# 	$6: ending subject 
# 	$7: logging level (e.g., ALL,TRACE,DEBUG,INFO,WARN,ERROR,FATAL,OFF
# 	$8: logfile 
# 	$9: job submit script 

# Example Usage 
# source submit_all_jobs.bash 26 0 1 1 1 875 INFO log.txt single_job_submit.bash


##############################################################################

# compute the last starting point
# this allows me to pass more interpretable inputs to test_for_loop.bash
# e.g., 1 5 10 makes more sense than 1 5 6 for running subjects 1-10
LAST_START=$(($6 - $5 + 1))

# iterate over subjects 
for i in `seq $4 $5 $LAST_START`
do
    # compute the end subject for the current job, start subject is subject i 
    END_SUBJ=$(($i + $5 -1))

    # compute the file suffix indicating subject for the job 
    SUFFIX=whs_${1}_${i}_${END_SUBJ}_test${2}

    # submit job 
    qsub -v whsim=$1,dotest=$2,isHPC=$3,start_subj=$i,end_subj=$END_SUBJ,logging=$7,logfile=$8 \
        -N ${SUFFIX} -e error_log/${SUFFIX}.err \
        -o output_log/${SUFFIX}.out $9
done

