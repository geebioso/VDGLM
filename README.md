# Intro

The code in this repository is what we used fit the VDGLM to the 1200 subject release 
of the Human Connectome Project (HCP) on the UCI High Performance Cluster (HPC). 
Note that this code is not intended to be production level code and will require 
considerable hacking to run. 

## Data 

Our data files and a description of our preprocessing pipelines can be found 
on our [Open Science Foundation page](https://osf.io/4rvbz/wiki/home/). 

Our analyses begin with data from the [Minimal Preprocessing Pipeline](https://www.sciencedirect.com/science/article/pii/S1053811913005053). We run 
out anlaysis on 875 subjects that had low head motion. From the original 
grayordinate data, we extracted the time series for 333 surface regions of 
interest (ROIs) based on the [Gordon et al. atlas](https://academic.oup.com/cercor/article/26/1/288/2367115). Additionally, we 
performed scrubbing and motion correction. The data from the working memory 
task are stored in the file `tc_WM_875subj.mat.mat`. This file contains the 
following fields: 
    
    motionX: estimated mean head motion for each subject and run 
    readme: short readme 
    subjs: subject ids
    tasks: run ids 
    tc: timecourse for each subject and run 
    
The files `design_WM_LR.mat` and `combdesign_WM_LR.mat` contain the 
corresponding uncombined (10 conditions) and combined (4 conditions) design 
matrices, respectively. Each file contains the fields: 

    conditions: condition labels
    rst: a structure that contains variables describing the time series. rst.mtx is the design matrix. 


The main files we use for
analysis are `tc_WM_875subj.mat` which contains the imaging time course and  
`design_WM_LR.mat` and `combdesign_WM_LR.mat`, which contain the  
uncombined design matrix (10 conditions) and combined design matrix (4 
conditions), respectively. 

# Code Overview 
    
This code will run the analyses. The directory structure of the code is as 
follows: 

* `code`: main directory 
    * `optimization`: contains VDGLM objective function, gradients, and constraint
    * `plotting`: functions for plotting results      
    * `stats`: functions for computing statistics     
    * `batchmode`: functions for fitting models in batch mode     
    * `utils`: utility functions    
    * `HPC_scripts`: scripts for running jobs on the HPC (these must be 
located on the HPC)      

I would recommend running code on a desktop machine to save the hassle of 
sending code to and from the HPC and compiling MATLAB code. 

**Important**: you will have to make several 
changes to get this code to run on your machine. You will have to edit the file
 `batchmode/set_analysis_options.m` to set the paths to the data and the design
 files. You will have to modify the function 
`batchmode/set_results_directory.m` to set the result directory path. You can 
supply an option that will automatically create the required subdirectories. 

The key files for running the analysis on your own machine are:     

`batchmode/run_subjects_in_parallel.m`: This function is used to set up and run 
all the analyses you would like to run. There are options provided in the 
function for running the null hypothesis test or a comparison between the 
VDGLM and GLM on HCP data. 
    
`batchmode/analyzedata_batch.m`: This is the main function for analyzing a 
batch of subjects and the compiled version of this function is called when we 
submit a job to the HPC. This function performs parameter estimation and computes
 out-of-sample log likelihood, model predictions, and BIC. 

`batchmode/set_analysis_options.m`: This file sets analysis parameters and 
is called by `analyzedata_batch.m`. Edit this function to specify the paths 
to the data and the design matrices. The function takes an analysis number as a
 parameter and returns the options for analysis. Some important options that 
are set are all data and output file paths, the models to run, whether to 
prewhiten, and whether to add head motion regressors. 
    
`batchmode/fit_models_cv.m`: This file is the workhorse function for fitting 
models. It fitsboth GLM and VDGLM models and can fit pre-whitened or 
non-pre-whitened models. Cross validation is done inside this function. 

`batchmode/combine_results.m`: This file combines the results from all single
jobs. 

The main statistical functions are: 

`stats/compute_cohens_d.m`: This computes cohen's d for each parameter of each
model. The output is stored in nii format in the directory `../ROI2NIfTI/files`
. This output is then visualized using the HCP work bench software. 

`stats/null_sample_hypothesis_test.m`: This will compute the null distribution 
by sampling datsets. Output can be used to test the expected false positive
rate for the VDGLM. This function is not ready to be used quite yet. 

The main plotting function is: 

`plotting/plot_simulations.m`: This function plots of model
comparisons, brain visualizations of model comparisons, parameter histograms, 
and model predictions versus the actual data. 

# Running the Analysis in the Paper 



1) First you will have to edit directory paths for the following files:     
    * `add_all_paths.m`: change  the path to a NIFTI installation     
    * `batchmode/set_analysis_options.m`: change the paths to the time course and the 
design matrix     
    * `batchmode/set_results_directory`: specify the paths where you would like to save 
model results, generated images, and ROI brain files 
    * `plotting/wb_global_variables.bash`: change the MAIN_FILE_DIRECTORY to be
equal to the directory where you store ROI brain files, change the IMAGE_DIRECTORY
to the same image directory specified in `set_results_directory.m`     
2) Next run the function `add_all_paths.m`, which will add the important 
code subdirectories to your path.      
3)  Run the following functions in order:         
    * `batchmode/run_subjects_in_parallel`: this file will perform one analysis
for each subject in parallel    
    * `batchmode/combine_results`: this file will join all the results into group level
files     
    * `stats/compute_cohens_d`: this function will compute group-level Cohen's d for
all contrasts specified in the paper     
    * `stats/compute_color_bounds`: this function will compute the appropriate color 
bounds and create color bars for the Cohen's d values just computed (colorbars 
 are saved in the image directory)     
    * `plotting/plot_simulations`: this will plot the scatter plot and any figures
used to check our analysis (model predictions, area plots, e.g. non-brain images)     

4) UNDER CONSTRUCTION: To plot the brain images in the paper, perform the 
following:     
    a)  Get rid of all MAC newlines so that windows bash can run the scripts:     
`cd plotting`
`sed -i 's/\r$//' format_wb_scripts_for_windows.bash`
	b) There are a few important files for how to format images:     
		* `wb_cont_format.bash`: this will format an image for wb_view on a 
non-fixed scale for the VDGLM     
	    * `wb_cont_format_mdl2.bash`: this will format an image for wb_view on 
a non-fixed scale for the GLM     
		* `wb_cont_format_fixed(_mdl2).bash`: these scripts will format images 
for wb_view on a fixed scale (_mdl2 is GLM )     
		These files can be fun as follows: `source wb_cont_format.bash 0.2` to 
format an image with a non-fixed color scheme with a Cohen’s d threshold of 0.2     
	c) There are also corresponding files for writing the images (these iterate
 over the thresholds we need): 
		```wb_cont_to_image.bash
		wb_cont_to_image_fixed_mdl2.bash
		wb_cont_to_image_fixed.bash
		wb_cont_to_image_mdl2.bash```    
		These files can be run as follows: `source wb_cont_to_image.bash`
	d) For each effect size threshold, you need to MANUALLY CREATE A SCENE FILE:     
		1) For example, for a small threshold, run the command: 
			```source wb_cont_format.bash 0.2 
			wb_view contrasts.spec``` 
		2) Create a scene file with the name cohens_d_whs26_[threshsize].scene.
 Thresh size much match Cohen’s d: {small:0.2, medium:0.5, large:0.8}    
			a) click the movie click board in workbench     
			b) for each MAP in workbench add a new scene, click add window        
    e) Once you have saved all the scene files, you can print images using to 
image commands as follows: `source wb_cont_to_image.bash` 

# Examples

Before trying any examples, make sure to set the paths by running 
`add_all_paths.m`.

## Setting up a simulation 
Here is the example output from running the function 
`set_analysis_options.m` with the parameters `whsim=26`, `isHPC=0`, 
`dotest=0`, and `LOG=LOG = log4m.getLogger('test_log.txt')`: 

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
    - Fit GLM1 on prewhitened data     

This code adds the additional VDGLM fitting step so that the whole procedure 
is now:    
- Fit GLM0    
    - prewhiten (using residuals of GLM0)     
    - Fit GLM1 on prewhitened data    
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

Most of these are excluded from the git repo to keep the size of the repo 
maneageable. You will need to make stand in directories for some of the 
output to save. `Results` and `images` and their subdirectories will need to 
exist for functions to run.  

<!-- 
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

-->
