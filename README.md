# Intro

The code in this repository will fit the VDGLM to the 1200 subject release 
of the Human Connectome Project (HCP) on the UCI High Performance Cluster (HPC). 

# Running on the HPC 

## Sending Code and Data to HPC 

To send code and the requisite data to the HPC, open up a terminal and type 
the following commands: 

```
source toHPC.bash
source send_data_to_HPC.bash
```

## Compiling on the HPC 

```
module load MATALB/r2017b
mcc -m analyzedata_batch.m -I /data/users/ggaut/VDGLM
```

## Submitting Jobs 

The main functions for submitting jobs on the HPC are `submit_all_jobs.bash` 
and `single_job_submit.bash`. To submit jobs, cd to the directory `/data/users/ggaut/VDGLM`
and type the command: 

```
source submit_all_jobs.bash [whsim] [dotest] [isHPC] [start_sub] [increment] [end_sub] [loglevel] [logfile] single_job_submit.bash
``` 

This command will run the analysis starting at subject `start_sub` and ending
 at subject `end_sub` incrementing by `increment` subjects. E.g., `1 1 5` will 
run subject 1,2,3,4,5 and `1 2 5` will run subjects 1,3,5. The other options 
are:    
* whsim: the simulation number (26 for CIfTI with prewhitening)     
* dotest: are we running in test mode (0/1)    
* isHPC: are we running on the HPC (0/1)     
* loglevel: one of {ALL, TRACE, DEBUG, INFO, WARN, ERROR, FATAL, OFF}     
* logfile: where to save log. is deprecated since I turned logging to file off 
in the code     

## HPC Results

The HPC Results directory is contained in `/pub/ggaut/VDGLM`. The Results 
structure will be the same as the directory structure on my local machine. 
To copy results from the HPC to my machine, use the `rsync` command. For 
example: 

```
rsync -av ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/batch_analyses/combined .
```

will copy the combined directory on the HPC to the current directory. The
results structure is as follows: 

`Results`: contains all results    
    `single_analyses`: results for analyses run all at once    
    `batch_analyses`: results for analyses run in batch    
        `single_jobs`: contains results from single jobs    
        `combined`: contains results combined from all single jobs results    
        `null_single_jobs`: contains results from single null sampling jobs    
        `null_combined`: contains results combined from all single null 
sampling jobs 

## Combining Results 

The function `combine_results.m` will combine all results from the directory 
`Results\batch_analyses\single_jobs` into single files and store them in 
`Results\batch_analyses\combined`.  

## To-do

1) Use job arrays rather than for-loop for submitting multiple jobs 

# Code Overview 
    
This code will run the analyses  **Important**: you will have to make several 
changes to get this code to run on your machine. You will have to edit the file
 `set_analysis_options.m` to set the paths to the data and the design files. 
The files that we are using for analysis 
are `tc_csf_wm_motion_out_globalMask5000_203subj.mat` for the OSUWB data and 
`all_tc_hcp.mat` for the HCP data. The design files we are using are 
`WorkingMem_mtx_3column_filter90s.mat` for the OSUWB data and 
`behav_Result_WM_LR.mat` and `behav_Result_WM_LR_comb.mat` for the HCP data. 
`behav_Result_WM_LR.mat` contains the uncombined design (10 conditions) and 
`behav_Result_WM_LR_comb.mat` contains the combined design (4 conditions). 

The key files are:    
    
`analyzedata_batch_v2.m`: This is the main function for analyzing a batch of 
subjects and the compiled version of this function is called when we submit a 
job to the HPC.     

`set_analysis_options_v2.m`: This file sets analysis parameters and is called 
by `analyzedata_batch_v2.m`. Edit this function to specify the paths to the 
data and the design matrices. The function takes an analysis number as a 
parameter and returns the options for analysis. Some important options that are
 set are all data and output file paths, the models to run,  whether to 
prewhiten, and whether to add head motion regressors. 
    
`fit_models_cv.m`: This file is the workhorse function for fitting models. It 
fitsboth GLM and VDGLM models and can fit pre-whitened or non-pre-whitened 
models. Cross validation is done inside this function. 

`plot_cerebral_cortex.m`: This function plots all results and output a 
structure `sims`. `sims` contains a field `simi` for each simulation `i` that 
was plotted.     

# Example Output    
    
If we load an output file (e.g., `paramswm_cc_whs1.mat`) it will load several variables:     

*  `allbicm [NS x R x P]`: contains the BIC scores for `NS` subjects, `R` regions, and `P` parameters     
*  `allllsm [NS x R x P]`: contains the log likelihood scores     
*  `bestmodelBIC [NS x R]`: contains the best models for each subject and region compute by BIC. The number indexes into the cell array `models` (e.g, `bestmodelBIC(s,r) = 1`, corresponds to `models{1}`)     
*  `bestmodelCV [NS x R]`: contains the best model for each subject and region computed by out-of-sample log likelihood     
*  `M`: number of models     
*  `K`: number of cross validation folds    
*  `savedwhsim`: inicates which analysis was just saved     
*  `models`: a cell array containing information about each model. The fields in models are:    
*  `code`: lists which model we are running and which columns are in the mean and variance design    
*  `meaneffects`: lists the mean effects in text    
*  `vareffects`: lists the variance effects in text    
*  `description`: text description of the model    
*  `meancols`: lists which columns of the full design are part of the mean design    
*  `varcols`:  lists which columns of the full design are part of the variance design    
*  `designlabels`: labels of the design matrix    
*  `inits`: initial values for the optimization    
*  `paramlabels`: cell of parameter descriptions    
*  `whsim`: analysis number    
*  `ismeanparam`: index into which parameters correspond to the mean    
*  `isvarparam`: index into which parameters correspond to the variance    
*  `isinterceptparam`: index into which parameters are intercepts    
*  `allparams`: all parameter values [NS x P x R]    
*  `allpredsm`: all mean predictions [T x S x R]    
*  `allpredsv`: all variance predictions [T x S x R]        
    
        For example:    
 
           models{1}    
                   
               ans =     
                   
                 struct with fields:    
                   
                               code: {{1×4 cell}  {1×2 cell}  'Var+Mean+Fix+Ins'}    
                        meaneffects: {'Intercept'  '2-back'  'Instruction'  'Fixation'}    
                         vareffects: {'Intercept'  '2-back'}    
                        description: 'Var+Mean+Fix+Ins'    
                           meancols: [1 3 4 5]    
                            varcols: [1 3]    
                       designlabels: {'Intercept'  '2-back'  'Instruction'  'Fixation'}    
                              inits: [0 0.1000 0.1000 0.1000 1 0.1000]    
                        paramlabels: {6×1 cell}    
                              whsim: 2    
                        ismeanparam: [1 1 1 1 0 0]    
                         isvarparam: [0 0 0 0 1 1]    
                   isinterceptparam: [1 0 0 0 1 0]    
                          allparams: [142×6×299 double]    
                          allpredsm: [395×142×299 single]    
                          allpredsv: [395×142×299 single]       

## Where to Modify Code

To modify the code to run on all voxels, you will have to edit the functions 
`set_analysis_options.m` and `analyzedata_cerebral_cortex.m`. In 
`set_analysis_options.m`, add another simulation with the same options that 
are output for the other simulations. In `analyzedata_cerebral_cortex.m`, you 
will have to change the data that is loaded (Line 13). The data needs to be 
size T x S x V where T is the length of the time series, S is the number of 
subjects, and V is the number of voxels. 


## Working Directory Structure

The directories on my local machine are: 

`Results`: contains all results     
    `single_analyses`: results for analyses run all at once        
    `batch_analyses`: results for analyses run in batch    
        `single_jobs`: contains results from single jobs     
        `combined`: contains results combined from all single jobs results     
        `null_single_jobs`: contains results from single null sampling jobs     
        `null_combined`: contains results combined from all single null sampling jobs     

`images`: contains figures and images 
    `contasts`: volumetric analysis contrasts 
    `wb_contrasts`: work bench contrasts (CIfTI) 
    `wb_effects`: work bench effect indicators 

`ROI2NIfTI`: contains nii files and functions for converting between vectors and nii     
    `files`: nii files used to create images in `images\contrasts` and `images\wb_contrasts`     
        `effects`: nii effect indicators used to create images in `images\wb_effects`            
    `dicm2nii`: contains code for converting between dicm and nii       

`onsetfiles`: contains the onsetfiles for HCP and OSUWB datasets 



