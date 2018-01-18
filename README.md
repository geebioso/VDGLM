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

* `Results`: contains all results    
    * `single_analyses`: results for analyses run all at once    
    * `batch_analyses`: results for analyses run in batch    
        * `single_jobs`: contains results from single jobs    
        * `combined`: contains results combined from all single jobs results    
        * `null_single_jobs`: contains results from single null sampling jobs    
        * `null_combined`: contains results combined from all single null 
sampling jobs 

## Combining Results 

The function `combine_results.m` will combine all results from the directory 
`Results\batch_analyses\single_jobs` into single files and store them in 
`Results\batch_analyses\combined`.  

## To-do

1) Use job arrays rather than for-loop for submitting multiple jobs 

# Code Overview 
    
This code will run the analyses. The directory structure of the code is as 
follows: 

* `code`: main directory 
    * `optimization`: contains VDGLM objective function, gradients, and constraint
    * `plotting`: functions for plotting results      
    * `stats`: functions for computing statistics     
    * `batchmode`: functions for fitting models in batch mode     
    * `utils`: utility functions    



**Important**: you will have to make several 
changes to get this code to run on your machine. You will have to edit the file
 `set_analysis_options.m` to set the paths to the data and the design files. 
The files that we are using for analysis 
are `tc_csf_wm_motion_out_globalMask5000_203subj.mat` for the OSUWB data and 
`all_tc_hcp.mat` for the HCP data. The design files we are using are 
`WorkingMem_mtx_3column_filter90s.mat` for the OSUWB data and 
`behav_Result_WM_LR.mat` and `behav_Result_WM_LR_comb.mat` for the HCP data. 
`behav_Result_WM_LR.mat` contains the uncombined design (10 conditions) and 
`behav_Result_WM_LR_comb.mat` contains the combined design (4 conditions). 

The key analysis files are:    
    
`batchmode/analyzedata_batch_v2.m`: This is the main function for analyzing a 
batch of subjects and the compiled version of this function is called when we 
submit a job to the HPC. This funciton performs parameter analysis and computes
 out-of-sample log likelihood,  model predictions, and BIC. 

`batchmode/set_analysis_options_v2.m`: This file sets analysis parameters and 
is called by `analyzedata_batch_v2.m`. Edit this function to specify the paths 
to the data and the design matrices. The function takes an analysis number as a
 parameter and returns the options for analysis. Some important options that 
are set are all data and output file paths, the models to run, whether to 
prewhiten, and whether to add head motion regressors. 
    
`batchmode/fit_models_cv.m`: This file is the workhorse function for fitting 
models. It fitsboth GLM and VDGLM models and can fit pre-whitened or 
non-pre-whitened models. Cross validation is done inside this function. 

`plot_cerebral_cortex.m`: This function plots all results and output a 
structure `sims`. `sims` contains a field `simi` for each simulation `i` that 
was plotted.     

The main statistical functions are: 

`stats/compute_cohens_d.m`: This computes cohen's d for each parameter of each
model. The output is stored in nii format in the directory `../ROI2NIfTI/files`
. This output is then visualized using the HCP work bench software. 

`stats/null_sample_hypothesis_test.m`: This will compute the null distribution 
by sampling datsets. Output can be used to test the expected false positive
rate for the VDGLM. This function is not ready to be used quite yet. 

The main plotting function is: 

`plotting/plot_simulations.m`: This function plots are plots of model
comparisons, brain visualizations of model comparisons, parameter histograms, 
and model predictions versus the actual data. 

# Examples

Before trying any examples, make sure to set the paths by running 
`add_all_paths.m`.

## Setting up a simulation 
Here is the example output from running the function 
`set_analysis_options_v2.m` with the parameters `whsim=26`, `isHPC=0`, 
`dotest=0`, and `LOG=LOG = log4m.getLogger(l'')`: 

```
opts = 

  struct with fields:

                whsim: 26
                  whs: 2
                    K: 10
                 seed: 1
              Tremove: 10
        doconstrained: 1
            prewhiten: 1
    var_log_transform: 0
                 TukN: 5
         multivariate: 0
              roifile: '/Users/Garren/Dropbox/FMRI/restingstatedata/tc_WM_875subj.mat'
           designfile: '/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/onsetfiles/behav_Result_WM_LR.mat'
       combdesignfile: '/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/onsetfiles/behav_Result_WM_LR_comb.mat'
            runmodels: {4×1 cell}
     output_directory: '/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/Results/batch_analyses/single_jobs'
            addmotion: 1
```

`opts.runmodels` contains information on the 4 models we are going to fit. 
`opts.runmodels{m}{1}` contains the mean design columns for model `m`. 
`opts.runmodels{m}{2}` contains the variance design columns, and 
`opts.runmodels{m}{3}` contains the model description such that 
`opts.runmodels{m}{3}{1}` contains the model name and 
`opts.rummodels{m}{3}{2}` contains the model name of the model that model `m` 
is dependent on (if applicable). 

### Model Dependencies 
The variance models are hard coded to depend on the output of a corresponding
mean model in order to speed up fitting time. This change is motivated by how 
prewhitening works. The prewhitening procedure works as follows: 
- Fit GLM0
    - prewhiten (using residuals of GLM0) 
    - Fit GLM1
This code adds the additional VDGLM fitting step so that the whole procedure 
is now:
- Fit GLM0
    - prewhiten (using residuals of GLM0) 
    - Fit GLM1
    - Fit VDGLM0
Encoding dependencies allows us to fit the mean models using OLS (which is 
faster than our optimization procedure). Model dependencies are hard-coded in 
the third entry of each entry of run_models. The entry is a cell of length two 
with the model name in the first cell element and the model that must run 
before the current model in the second cell elements. e.g. {'Var+Mean', 'Mean'}
 indicates that the 'Mean' model must have been run for the  'Var+Mean' to run.
No dependency is indicated by an empty string.  

## Fitting a single subject 
To fit subject 1 in test mode (test mode only supports subjects 1-5), run the 
following command: 

```
analyzedata_batch_v2(26, 1, 0, 1, 1, 'INFO', '')
```

The command line output will be: 

```
INFO:Running Simulation 26
INFO:    whsim: 26.000000
INFO:    whs: 2
INFO:    K: 10.000000
INFO:    seed: 1
INFO:    Tremove: 10
INFO:    doconstrained: 1
INFO:    prewhiten: 1
INFO:    var_log_transform: 0
INFO:    TukN: 5.000000
INFO:    multivariate: 0
INFO:    roifile: /Users/Garren/Dropbox/FMRI/restingstatedata/tc_WM_875subj.mat
INFO:    designfile: /Users/Garren/Dropbox/FMRI/Projects/varianceGLM/onsetfiles/behav_Result_WM_LR.mat
INFO:    combdesignfile: /Users/Garren/Dropbox/FMRI/Projects/varianceGLM/onsetfiles/behav_Result_WM_LR_comb.mat
INFO:    output_directory: /Users/Garren/Dropbox/FMRI/Projects/varianceGLM/Results/batch_analyses/single_jobs_test
INFO:    addmotion: 1
INFO:    runmodels:
INFO:        model: Var+Mean, depends on: Mean
INFO:        model: Mean, depends on: 
INFO:        model: Var, depends on: Intercept
INFO:        model: Intercept, depends on: 
INFO:Loading data: /Users/Garren/Dropbox/FMRI/restingstatedata/tc_WM_875subj.mat
INFO:Converting data and transforming to z-scores
INFO:Creating convolved design matrix for all experimental variables
INFO:    dotest = 1
INFO:Creating design for model 1 of 4 (Var+Mean)
INFO:Creating design for model 2 of 4 (Mean)
INFO:Creating design for model 3 of 4 (Var)
INFO:Creating design for model 4 of 4 (Intercept)
INFO:Fitting models ..... working on subject 1 of [1:1]

INFO:Total Run Time = 8.70
INFO:saving file /Users/Garren/Dropbox/FMRI/Projects/varianceGLM/Results/batch_analyses/single_jobs_test/paramswm_whs26_batch_1_to_1_test
```


## Output  
    
If we load an output file (e.g., i`paramswm_whs26_batch_1_to_1_test`), the
following variables  will be loaded:     

*  `allbicm [NS x R x P]`: contains the BIC scores for `NS` subjects, `R` 
regions, and `P` parameters     
*  `allllsm [NS x R x P]`: contains the log likelihood scores     
*  `bestmodelBIC [NS x R]`: contains the best models for each subject and 
region compute by BIC. The number indexes into the cell array `models` (e.g, 
`bestmodelBIC(s,r) = 1`, corresponds to `models{1}`)     
*  `bestmodelCV [NS x R]`: contains the best model for each subject and region 
computed by out-of-sample log likelihood     
*  `M`: number of models     
*  `K`: number of cross validation folds    
*  `savedwhsim`: inicates which analysis was just saved     
*  `models`: a cell array containing information about each model. The fields 
in models are:    
*  `code`: lists which model we are running and which columns are in the mean 
and variance design    
*  `meaneffects`: lists the mean effects in text    
*  `vareffects`: lists the variance effects in text    
*  `description`: text description of the model    
*  `meancols`: lists which columns of the full design are part of the mean 
design    
*  `varcols`:  lists which columns of the full design are part of the variance 
design    
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

``` 
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

```

## Where to Modify Code

To modify the code to run your machine, you will have to edit the functions
`set_analysis_options_v2.m` and `utils/set_results_directory.m`. The paths to 
the data must be updated in `set_analysis_options_v2.m`. You must update
`utils/set_results_directory.m` to indicate where you want your results to be 
saved. Make sure that the directory structure is the same as the `Results` 
directory structure outlined below.

## Working Directory Structure

The directories on my local machine are: 

* `Results`: contains all results     
    * `single_analyses`: results for analyses run all at once        
    * `batch_analyses`: results for analyses run in batch    
        * `single_jobs`: contains results from single jobs     
        * `combined`: contains results combined from all single jobs results     
        * `null_single_jobs`: contains results from single null sampling jobs     
        * `null_combined`: contains results combined from all single null
 sampling jobs     

* `images`: contains figures and images 
    * `contasts`: volumetric analysis contrasts 
    * `wb_contrasts`: work bench contrasts (CIfTI) 
    * `wb_effects`: work bench effect indicators 

* `ROI2NIfTI`: contains nii files and functions for converting between vectors 
and nii     
    * `files`: nii files used to create images in `images\contrasts` and 
`images\wb_contrasts`     
        * `effects`: nii effect indicators used to create images in 
`images\wb_effects`            
    * `dicm2nii`: contains code for converting between dicm and nii       

* `onsetfiles`: contains the onsetfiles for HCP and OSUWB datasets 

Most of these are excluded from the git repo to keep the size of the repo 
maneageable. You will need to make stand in directories for some of the 
output to save. `Results` and `images` and their subdirectories will need to 
exist for functions to run.  


