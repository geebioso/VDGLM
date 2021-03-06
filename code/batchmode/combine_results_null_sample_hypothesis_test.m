function[] = combine_results_null_sample_hypothesis_test( whsim, isHPC , dotest, logging, set_up_directory_structure )

% This function will read through the directory containing results from
% single subjects and aggregate them into a combined result. Due to the
% size of the data, models will be stored separately (unlike the volumetric
% results). The function load_results.m will take this into account and
% when loaded, the results will be in the same format as the volumetric
% results. 

% INPUT: 
%   numeric whsim: which simulation 
%   bool isHPC: are we doing this on the HPC? 
%   bool dotest: are we in test mode? 

% 1-17-2018: Re-writing code to reflect that we only process one subject at
% a time. At first I wanted it to process batches of subjects, but I never
% run things that way. 

% create the logger reference:
LOG = log4m.getLogger('test_log.txt');
LOG.setCommandWindowLevel(LOG.(logging));
LOG.setLogLevel(LOG.OFF);

%% Set up input and output directories 
[results_directory, images_directory, ROI2NIfTI_directory] = set_results_directory( isHPC, set_up_directory_structure );

input_directory = fullfile( results_directory, 'batch_analyses', 'null_single_jobs'); 
output_directory = fullfile( results_directory, 'batch_analyses', 'null_combined'); 
if dotest
   input_directory = [input_directory '_test'];  
   output_directory = [output_directory '_test'];  
end

LOG.info('INFO',sprintf('Input directory: %s', input_directory)); 
LOG.info('INFO',sprintf('Output directory: %s\n', output_directory)); 

%% Get list of all single job names 
files = dir(input_directory);
NF = length(files);
filenames = cell(NF, 1);

% get rid of hidden files 
for f = 1:NF
    filenames{f} = files(f).name; 
    if filenames{f}(1) == '.'
      filenames{f} = [];  
   end
end

ii = cellfun( @(x) ~isempty(x), filenames, 'UniformOutput', true);
filenames = filenames(ii); 

% get only files that are from the appropriate simulation 
ii = cellfun( @(x)  contains( x, sprintf('whs%d', whsim)) , filenames, 'UniformOutput', true); 
filenames = filenames( ii ); 

NF = length(filenames); 

%% Get list of all subjects
% The list will be used to make sure that aggregation of jobs into a single 
% result structure preserves subject identity 

bad_files = cell(1,2);  % list of subjects that are processed twice (was an issue at one point) 
all_subjs = cell(NF, 1); 
sub_nums = zeros(NF, 1); 
for f = 1:NF

    filenow = fullfile(input_directory, filenames{f});
    
    % parse subject number
    parts = split(filenames{f}, '_');
    sub_num = str2num(parts{5});
    load(filenow, 'subjs');
    
    % keep track of which subjects have been processed
   
    
        all_subjs{f} = subjs{1};
    
    sub_nums(f) = sub_num; 
end
     
[~, ii] = sort(sub_nums); 
sub_nums = sub_nums(ii); 
all_subjs = all_subjs(ii); 
filenames = filenames(ii); 

if LOG.commandWindowLevel < 3
    TAB = table(all_subjs, filenames, sub_nums)
end

bad_files = bad_files(2:end, :); 
NS = length(all_subjs);

%% Load one subject to get model storage structure
load(filenow);
R = size(models{1}.allparams, 2);
M = length(models);
Nsamp = size(models{1}.allparams, 3); 

%% Create master storage
allbicm = zeros(NS, R, Nsamp, M);
allllsm = zeros(NS, R, Nsamp, M);
allbadchol = zeros(NS, R, M); 

all_models = models;
for m = 1:M
    P = size(all_models{m}.allparams, 4);
    all_models{m}.allparams = zeros(NS, R, Nsamp, P);
    all_models{m}.motionparams = cell(NS, R); 
end

%% Iterate over files and add to master storage
for f = 1:NF
    
    filenow = fullfile(input_directory, filenames{f});
    sub_num_now = sub_nums(f); 
    fprintf('adding %s\n', filenow); 
    
    % load results from the current file and store
    results = load(filenow);
    bicm = results.allbicm;
    llsm = results.allllsm;
    bmBIC = results.bestmodelBIC;
    bmCV = results.bestmodelCV;
    subjsnow = results.subjs;
    subjnow = results.subjs{1};
    modelsnow = results.models;
    badchol = results.allbadchol; 
    
       % For each subject insert the results into master storage using the unique
    % subject identifier
    if length(subjsnow) > 1
       error('combine_results.m only works with batch sizes of 1');  
    end
    
    all_idx = find(ismember(all_subjs, subjnow));  % index into master list of subject
    now_idx = 1; % only works if 
    
    if LOG.commandWindowLevel < 3
        LOG.debug( 'DEBUG', sprintf('check: %s=%s, filename=%s, sub_num=%d', all_subjs{all_idx}, subjnow, filenames{f}, sub_num_now)); 
    elseif LOG.commandWindowLevel >= 3
        LOG.info('INFO', sprintf('%s', filenames{f})); 
    end
    
    % at each subject form the current file to master storage using unique
    % subject identifier
    for s = 1:length(subjsnow)
        
            subjnow = subjsnow{s};
            for m = 1:M
                all_models{m}.allparams( all_idx, :, :, :) = modelsnow{m}.allparams( now_idx, :, :, : );
                
                allbicm(all_idx,:,:,m) = bicm(now_idx,:,:,m);
                allllsm(all_idx,:,:,m) = llsm(now_idx,:,:,m);
                allbadchol(all_idx,:,m) = badchol(now_idx,:,m);
            end
    end
    
  
end

% compute best BIC and best CV 
[~, bestmodelBIC] = min(allbicm, [], 4);
[~, bestmodelCV ] = max(allllsm, [], 4);

% Save it all! Models saved separately 
fprintf('\n'); 
for m = 1:M
    filename = fullfile( output_directory, sprintf('whs%d_allmodels_%d', whsim, m));
    all_models_now = all_models{m};
    fprintf('Saving all_models %d\n', m); 
    save( filename, 'all_models_now' , '-v7.3');
end

fprintf('Saving allbicm\n');
save( fullfile( output_directory, sprintf('whs%d_allbicm.mat', whsim)), 'allbicm');
fprintf('Saving allllsm\n');
save( fullfile( output_directory, sprintf('whs%d_allllsm.mat', whsim)), 'allllsm');
fprintf('Saving bestmodelBIC\n');
save( fullfile( output_directory, sprintf('whs%d_bestmodelBIC.mat', whsim)), 'bestmodelBIC');
fprintf('Saving bestmodelCV\n');
save( fullfile( output_directory, sprintf('whs%d_bestmodelCV.mat', whsim)), 'bestmodelCV');
fprintf('Saving all_subjs\n');
save( fullfile( output_directory, sprintf('whs%d_all_subjs.mat', whsim)), 'all_subjs');
fprintf('Saving sub_nums\n');
save( fullfile( output_directory, sprintf('whs%d_sub_nums.mat', whsim)), 'sub_nums');



