% function[] = combine_results_null_sample_hypothesis_test( isHPC , dotest )


% This function will read through the directory containing results from
% single subjects and aggregate them into a combined result. Due to the
% size of the data, models will be stored separately (unlike the volumetric
% results). The function load_results.m will take this into account and
% when loaded, the results will be in the same format as the volumetric
% results. 

% INPUT: 
%   bool isHPC: are we doing this on the HPC? 
%   bool dotest: are we in test mode? 

%% Set I/O Directories 

if isHPC
    results_directory = '/pub/ggaut/VDGLM/Results';
else
    results_directory = fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'Projects',...
        'varianceGLM', 'Results');     
end

input_directory = fullfile( results_directory, 'batch_analyses', 'null_single_jobs'); 
output_directory = fullfile( results_directory, 'batch_analyses', 'null_combined'); 
if dotest
   input_directory = [input_directory '_test'];  
   output_directory = [output_directory '_test'];  
end

fprintf('Input directory: %s\n', input_directory); 
fprintf('Output directory: %s\n\n', output_directory); 


%% Get list of all single job result filenames 
files = dir(input_directory);
NF = length(files);
filenames = cell(NF, 1);

include_transitions = 0; 

for i = 1:NF
    filenames{i} = files(i).name;
end

if strcmp(filenames{1}, '.')
    filenames = filenames(3:end);
    NF = NF - 2;
end

if strcmp(filenames{1}, '.DS_Store')
    filenames(1) = [];
    NF = NF - 1; 
end

%% Get list of all subjects
% The list will be used to make sure that aggregation of jobs into a single 
% result structure preserves subject identity 
bad_files = cell(1,2); 
for f = 1:NF
    filenow = fullfile(input_directory, filenames{f});
    load(filenow, 'subjs');
    
    % keep track of which subjects have been processed 
    if f == 1
        all_subjs = subjs';
    else
        sub_already_in_list = ismember(all_subjs, subjs); 
        if ~any(sub_already_in_list)
           all_subjs = [all_subjs; subjs'];
        else % find bad filename 
           bf_temp = cell(1,2); 
           bf_temp{1} = filenow; 
           bf_temp{2} = fullfile( input_directory, filenames{sub_already_in_list}); 
           bad_files = [bad_files; bf_temp];
        end
     
    end
end
bad_files = bad_files(2:end, :); 
NS = length(all_subjs);

%% Load one batch to get model storage structure
load(filenow);
R = size(models{1}.allparams, 2);
M = length(models);
Nsamp = size(models{1}.allparams, 3); 

%% Create master storage

allbicm = zeros(NS, R, Nsamp, M);
allllsm = zeros(NS, R, Nsamp, M);

all_models = models;
for m = 1:M
    P = size(models{m}.allparams, 4);
    all_models{m}.allparams = zeros(NS, R, Nsamp, P);
end

%% Iterate over files and add to master storage
for f = 1:NF
    
    filenow = fullfile(input_directory, filenames{f});
    fprintf('adding %s\n', filenow); 
    
    % load results from the current file and store
    results = load(filenow);
    bicm = results.allbicm;
    llsm = results.allllsm;
    bmBIC = results.bestmodelBIC;
    bmCV = results.bestmodelCV;
    subjsnow = results.subjs;
    modelsnow = results.models;
    
    % at each subject form the current file to master storage using unique
    % subject identifier
    for s = 1:length(subjsnow)
        
            subjnow = subjsnow{s};
            all_idx = find(ismember(all_subjs, subjnow));  % index into master list of subject
            now_idx = find(ismember(subjsnow, subjnow));   % index into the current list of subjects (from file we are adding to the master list)
            
            for m = 1:M
                all_models{m}.allparams( all_idx, :, :, :) = modelsnow{m}.allparams( now_idx, :, :, : );
                
                allbicm(all_idx,:,:,m) = bicm(now_idx,:,:,m);
                allllsm(all_idx,:,:,m) = llsm(now_idx,:,:,m);
            end
    end
end

% compute best BIC and best CV 
[~, bestmodelBIC] = min(allbicm, [], 4);
[~, bestmodelCV ] = max(allllsm, [], 4);

% Save it all! Models saved separately 
fprintf('\n'); 
for m = 1:M
    filename = fullfile( output_directory, sprintf('allmodels_%d', m));
    all_models_now = all_models{m};
    fprintf('Saving all_models %d\n', m); 
    save( filename, 'all_models_now' , '-v7.3');
end

fprintf('Saving allbicm\n');
save( fullfile( output_directory, 'allbicm.mat'), 'allbicm');
fprintf('Saving allllsm\n');
save( fullfile( output_directory, 'allllsm.mat'), 'allllsm');
fprintf('Saving bestmodelBIC\n');
save( fullfile( output_directory, 'bestmodelBIC.mat'), 'bestmodelBIC');
fprintf('Saving bestmodelCV\n');
save( fullfile( output_directory, 'bestmodelCV.mat'), 'bestmodelCV');
fprintf('Saving all_subjs\n');
save( fullfile( output_directory, 'all_subjs.mat'), 'all_subjs');



