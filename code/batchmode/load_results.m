function [models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, sub_nums] = ...
    load_results(results_directory, whsim, dotest, LOG, whmode)

% INPUT: 
%   results_directory: path to where we store results 
%   whsim: which simulation to run 
%   dotest: test mode? 
%   LOG: log4 logger 
%   whmode: are we fitting the real data or doing null distribution
%       sampling ('analyze' or 'null') 

LOG.info('INFO', sprintf('Loading Results, whsim = %d', whsim));
if whsim > 25
    
    switch whmode
        case 'analyze'
            input_directory = fullfile(results_directory, 'batch_analyses', 'combined');
        case 'null'
            input_directory = fullfile(results_directory, 'batch_analyses', 'null_combined');
    end
    if dotest
        input_directory = [input_directory '_test'];
    end
    
    % get numbers of all_models files
    files = dir(input_directory);
    NF = length(files);
    filenames = cell(NF, 1);
    for i = 1:NF
        filenames{i} = files(i).name;
    end
    
    % get simulation file names
    % get only files that are from the appropriate simulation
    ii = cellfun( @(x)  contains( x, sprintf('whs%d', whsim)) , filenames, 'UniformOutput', true);
    filenames = filenames( ii );
    
    % get model files
    ii = cellfun( @(x) regexp( x , 'allmodels' ), filenames, 'UniformOutput', 0);
    ii = cellfun( @(x) ~isempty(x), ii);
    filenames = filenames(ii);
    NF = length(filenames);
    
    % parse model number from each filename
    mods = zeros(NF, 1);
    for f = 1:NF
        temp = split(filenames{f}, '_');
        temp = split(temp{3}, '.');
        mods(f) = str2num(temp{1});
    end
    
    M = length(mods);
    models = cell(M,1);
    
    for i = 1:M
        m = mods(i);
        filename = fullfile( input_directory, sprintf( 'whs%d_allmodels_%d.mat', whsim, m));
        LOG.info('INFO', sprintf('\t%s', filename));
        load(filename);
        models{i} = all_models_now;
        
        models{i}.allparams = single(models{i}.allparams);
        
        if isfield(models{i}, 'allpredsm')
            models{i}.allpredsm = single(models{i}.allpredsm);
            models{i}.allpredsv = single(models{i}.allpredsv);
        end
    end
    
    filenames = {'allbicm.mat'; ...
        'allllsm.mat';...
        'bestmodelBIC.mat'; ...
        'bestmodelCV.mat'; ...
        'all_subjs.mat'; ...
        'sub_nums.mat'
        };
    
    for f = 1:length(filenames)
        filename = fullfile( input_directory, sprintf('whs%d_%s', whsim, filenames{f}));
        load(filename);
        LOG.info('INFO', sprintf('\t%s', filename));
    end
    
else
    
    input_directory = fullfile( results_directory, 'single_analyses');
    if dotest
        filenm = fullfile(input_directory, sprintf( 'paramswm_cc_whs%d_test' ,whsim));
    else
        filenm = fullfile(input_directory, sprintf( 'paramswm_cc_whs%d' ,whsim));
    end
    load(filenm);
    LOG.info('INFO', sprintf('\t%s', filename));
    all_subjs = [];
    
end

end