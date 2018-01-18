function [ opts, dotest ] = set_analysis_options_v2(whsim, isHPC, dotest, LOG)

% NOTE!!!!!!
%   Right now this code will not run for any sim numbers except 26 and 27.
%   I got rid of the other simulations to make the code more readable and
%   because we will most likely not run many in the future--add them later.
%   
%   This functions changes from the original function
%   set_analysis_options.m as follows: 
% 
%   I'm hard coding in a way to make the variance models depend on the
%   output of the mean models in order to speed up fitting time. The
%   motivation for this change is how prewhitening works: 
%   - Fit GLM0
%       - prewhiten
%       - Fit GLM1
%   I am going to add an additional step so that the procedure becomes: 
%   - Fit GLM0
%       - prewhiten
%       - Fit GLM1
%       - Fit VDGLM0
%   Then we can fit the mean models using OLS (which is faster than our
%   optimization procedure). 
%   Model dependencies are hard-coded in the following way: 
%       Int  <- Var
%       Mean <- Var+Mean 
%   Where <- indicates that the first model must run for the second to run.
%   
%   The way we represent this dependency is in the third entry of each entry of
%   run_models. The entry is a cell of length two with the model name in
%   the first cell element and the model that must run before the current
%   model in the second cell elements. e.g. {'Var+Mean', 'Mean'} indicates 
%   Mean <- Var+Mean. 
% 
%   The set-up is reflected here and the change in the code is reflected in
%   the function fit_models.m 


% INPUT:
%   whsim: which analysis do we want to run?
%   isHPC: switch for running on a local machine or on the UCI HPC 
%   dotest: switch for running in test mode (discards a lot of data) 

% OUTPUT:
%   dotest: switch for running in test mode (discards a lot of data) 
%   opts: structure that contains 
%       whs: which dataset
%       K: number of cross validation folds
%       seed: random seed
%       Tremove: how many TRs to remove from design matrix
%       doconstrained: flag to do constrained or unconstrained optimization
%       prewhiten: prewhiten the data?
%       var_log_transform: do we transform the variance to be on a log scale?
%       TukN: How long do we want our Tukey Filter to be?
%       roifile: data file to load
%       designfile: design file to load
%       runmodels: list of models that we are running


%% Set directories for OSUWB and HCP data
if isHPC
    % hcproifile =  '/pub/ggaut/HCP/all_tc_hpc.mat';
    hcproifile =  '/pub/ggaut/HCP/tc_WM_875subj.mat';
    
    % design files
    hcpdesignfile =  '/pub/ggaut/HCP/behav_Result_WM_LR.mat';
    hcpcombdesignfile =  '/pub/ggaut/HCP/behav_Result_WM_LR_comb.mat';
    
    results_directory = '/pub/ggaut/VDGLM/Results'; 
else
    % osu
    osuroifile = fullfile(getenv('HOME'), 'Dropbox','FMRI', 'restingstatedata',...
        'tc_motion_CSF_WM_out_WorkingMem177.mat');
    osudesignfile = fullfile(getenv('HOME'), 'Dropbox','FMRI', 'Projects', ...
        'varianceGLM', 'onsetfiles','WorkingMem_mtx_3column_filter90s.mat');
    osucombdesignfile = '';
    
    % hcp
    if ismember(whsim, [26, 27])
        hcproifile =  fullfile(getenv('HOME'), 'Dropbox','FMRI', ...
            'restingstatedata','tc_WM_875subj.mat');
    else
        hcproifile =  fullfile(getenv('HOME'), 'Dropbox','anatomical',...
            'ICAresults','all_tc_hcp.mat');
    end
    hcpdesignfile =  fullfile( getenv('HOME'), 'Dropbox','FMRI', 'Projects', ...
        'varianceGLM', 'onsetfiles', 'behav_Result_WM_LR.mat');
    hcpcombdesignfile =  fullfile( getenv('HOME'), 'Dropbox','FMRI', 'Projects', ...
        'varianceGLM', 'onsetfiles', 'behav_Result_WM_LR_comb.mat');
    
    results_directory = fullfile(getenv('HOME'), 'Dropbox', 'FMRI', ...
        'Projects', 'varianceGLM', 'Results');
end

%% Pick simulation
if (whsim==26)
    % HCP DATA BASIC CIFTI (NO TD, PW) with FIXATION IN VARIANCE
    whs     = 2; % which data set? 8= Ohio State Well Being (edge voxels removed); 0 = HCP WM
    K       = 10; % Number of folds in CV; K=0 do not CV
    seed    = 1;
    Tremove = 10; % How many time steps to remove from data?
    doconstrained = true; % 1 = Use fmincon; 0 = fminunc
    
    prewhiten = true;
    var_log_transform = false;
    TukN = 5;
    multivariate = false;
    addmotion = true; 
    
    roifile = hcproifile;
    designfile = hcpdesignfile;
    combdesignfile = hcpcombdesignfile;
    output_directory = fullfile( results_directory, 'batch_analyses', 'single_jobs'); 
    
    runmodels  = {
        {{ 'Intercept' , '0-back', '2-back' , 'Instruction' } , { 'Intercept' , '0-back', '2-back', 'Instruction' }    , {'Var+Mean', 'Mean'} },
        {{ 'Intercept' , '0-back', '2-back' , 'Instruction' } , { 'Intercept' }                         , {'Mean', ''} },
        {{ 'Intercept' }                                      , { 'Intercept' , '0-back', '2-back', 'Instruction'  }   , {'Var', 'Intercept'}},
        {{ 'Intercept'                            }           , { 'Intercept'                       }   , {'Intercept', ''}   },
        };
elseif (whsim==27)
    % HCP DATA BASIC (NO TD, NO PW) with FIXATION IN VARIANCE
    whs     = 2; % which data set? 8= Ohio State Well Being (edge voxels removed); 0 = HCP WM
    K       = 10; % Number of folds in CV; K=0 do not CV
    seed    = 1;
    Tremove = 10; % How many time steps to remove from data?
    doconstrained = false; % 1 = Use fmincon; 0 = fminunc
    
    prewhiten = false;
    var_log_transform = false;
    TukN = NaN;
    multivariate = false;
    addmotion = true; 
    
    roifile = hcproifile;
    designfile = hcpdesignfile;
    combdesignfile = hcpcombdesignfile;
    output_directory = fullfile( results_directory, 'batch_analyses', 'single_jobs'); 
    
    runmodels  = {
        {{ 'Intercept' , '0-back', '2-back' , 'Instruction' } , { 'Intercept' , '0-back', '2-back', 'Instruction' }    , {'Var+Mean', 'Mean'} },
        {{ 'Intercept' , '0-back', '2-back' , 'Instruction' } , { 'Intercept' }                         , {'Mean', ''} },
        {{ 'Intercept' }                                      , { 'Intercept' , '0-back', '2-back', 'Instruction'  }   , {'Var', 'Intercept'}},
        {{ 'Intercept'                            }           , { 'Intercept'                       }   , {'Intercept', ''}   },
        };
end

%% Change Output Directory to Test if Applicable 
if dotest
   output_directory = [output_directory '_test'];  
end

%% Change doubles to ints for printing

whs = int64(whs);
seed = int64(seed);
Tremove = int64(Tremove);

%% Assign Options to Structure

opts = struct();
opts.whsim = whsim;
opts.whs = whs;
opts.K = K;
opts.seed = seed;
opts.Tremove = Tremove;
opts.doconstrained = doconstrained;
opts.prewhiten = prewhiten;
opts.var_log_transform = var_log_transform;
opts.TukN = TukN;
opts.multivariate = multivariate;
opts.roifile = roifile;
opts.designfile = designfile;
opts.combdesignfile = combdesignfile;
opts.runmodels = runmodels;
opts.output_directory = output_directory;
opts.addmotion = addmotion;     

%% Print Options
LOG.info( 'INFO', sprintf('Running Simulation %d', whsim));

fields = fieldnames(opts);

for f = 1:length(fields)
    val = opts.(fields{f});
    if isa(val, 'char')
        LOG.info('INFO', sprintf('\t%0.20s: %s', fields{f}, val));
    elseif isa(val, 'double')
        LOG.info('INFO', sprintf('\t%0.20s: %f', fields{f}, val));
    elseif or(isa(val, 'logical'), isa(val, 'integer'))
        LOG.info('INFO', sprintf('\t%0.20s: %d', fields{f}, val));
    end
    
end

LOG.info('INFO', sprintf('\trunmodels:'))
for i = 1:length(runmodels)
    LOG.info('INFO', sprintf('\t\tmodel: %s, depends on: %s', opts.runmodels{i}{3}{1}, opts.runmodels{i}{3}{2})); 
end

if and(isnan(TukN), prewhiten)
   LOG.error('ERROR', 'Change Tukey Taper window or don''t prewhiten'); 
end

