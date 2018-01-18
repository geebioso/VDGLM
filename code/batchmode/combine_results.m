function[] = combine_results( isHPC , dotest, logging)

% This function will read through the directory containing results from
% single subjects and aggregate them into a combined result. Due to the
% size of the data, models will be stored separately (unlike the volumetric
% results). The function load_results.m will take this into account and
% when loaded, the results will be in the same format as the volumetric
% results. 

% INPUT: 
%   bool isHPC: are we doing this on the HPC? 
%   bool dotest: are we in test mode? 

% 1-17-2018: Re-writing code to reflect that we only process one subject at
% a time. At first I wanted it to process batches of subjects, but I never
% run things that way. 

% create the logger reference:
LOG = log4m.getLogger('');
LOG.setCommandWindowLevel(LOG.(logging));
LOG.setLogLevel(LOG.OFF);


%% Set up input and output directories 
[results_directory] = set_results_directory( isHPC ) 

input_directory = fullfile( results_directory, 'batch_analyses', 'single_jobs'); 
output_directory = fullfile( results_directory, 'batch_analyses', 'combined'); 
if dotest
   input_directory = [input_directory '_test'];  
   output_directory = [output_directory '_test'];  
end

fprintf('Input directory: %s\n', input_directory); 
fprintf('Output directory: %s\n\n', output_directory); 

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
    sub_num = str2num(parts{4});
    load(filenow, 'subjs');
    
    % keep track of which subjects have been processed
   
    
        all_subjs{f} = subjs{1};
    
    sub_nums(f) = sub_num; 
end

%         sub_already_in_list = ismember(all_subjs, subjs); 
%         if ~any(sub_already_in_list)
%             all_subjs = [all_subjs; subjs];
%             sub_nums = [sub_nums; sub_num]; 
%         else % find bad filename 
%            bf_temp = cell(1,2); 
%            bf_temp{1} = filenow; 
%            bf_temp{2} = fullfile( input_directory, filenames{sub_already_in_list}); 
%            bad_files = [bad_files; bf_temp];
%         end
     
[~, ii] = sort(sub_nums); 
sub_nums = sub_nums(ii); 
all_subjs = all_subjs(ii); 
filenames = filenames(ii); 

if LOG.commandWindowLevel < 3
    TAB = table(all_subjs, filenames, sub_nums)
end

bad_files = bad_files(2:end, :); 
NS = length(all_subjs);

%% Load one batch to get model storage structure
load(filenow);
R = size(models{1}.allpredsm, 3);
T = size(models{1}.allpredsm, 1);
M = length(models);

%% Create master storage
allbicm = zeros(NS, R, M);
allllsm = zeros(NS, R, M);
bad_subjs = cell(1); 

all_models = models;
for m = 1:M
    P = size(all_models{m}.allparams, 2);
    all_models{m}.allparams = zeros(NS, P, R);
    all_models{m}.motionparams = cell(NS, R); 
    all_models{m}.allpredsm = zeros(T, NS, R);
    all_models{m}.allpredsv = zeros(T, NS, R);
    
end

%% Iterate over files and add to master storage
LOG.info('INFO', sprintf('Adding from directory %s', input_directory))
for f = 1:NF
    
    filenow = fullfile(input_directory, filenames{f});
    sub_num_now = sub_nums(f); 
    
    % load results from the current file and store
    results = load(filenow);
    bicm = results.allbicm;
    llsm = results.allllsm;
    bmBIC = results.bestmodelBIC;
    bmCV = results.bestmodelCV;
    subjsnow = results.subjs;  
    subjnow = results.subjs{1};
    modelsnow = results.models;
    
    % For each subject insert the results into master storage using the unique
    % subject identifier
    if length(subjsnow) > 1
       error('combine_results.m only works with batch sizes of 1');  
    end
    
    all_idx = find(ismember(all_subjs, subjnow));  % index into master list of subject
    now_idx = 1; % only works if 
    
    if sub_num_now == 16
       x = 1;  
    end
    
    if LOG.commandWindowLevel < 3
        LOG.debug( 'DEBUG', sprintf('check: %s=%s, filename=%s, sub_num=%d', all_subjs{all_idx}, subjnow, filenames{f}, sub_num_now)); 
    elseif LOG.commendWindow >= 3
        LOG.info('INFO', sprintf('%s', filenames{f})); 
    end
    for m = 1:M
        all_models{m}.allparams( all_idx, :, :) = modelsnow{m}.allparams( now_idx, :, : );
        all_models{m}.motionparams( all_idx, : ) = modelsnow{m}.motionparams( now_idx, : );
        all_models{m}.allpredsm( :, all_idx, :) = modelsnow{m}.allpredsm( :, now_idx, : );
        all_models{m}.allpredsv( :, all_idx, :) = modelsnow{m}.allpredsv( :, now_idx, : );
        
        allbicm(all_idx,:,m) = bicm(now_idx,:,m);
        allllsm(all_idx,:,m) = llsm(now_idx,:,m);
    end
    
end

%% compute best BIC and best CV 
[~, bestmodelBIC] = min(allbicm, [], 3);
[~, bestmodelCV ] = max(allllsm, [], 3); 


%% Save it all! Models saved separately 
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
fprintf('Saving sub_nums\n');
save( fullfile( output_directory, 'sub_nums.mat'), 'sub_nums');



