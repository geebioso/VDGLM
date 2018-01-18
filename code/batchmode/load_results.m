function [models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, bad_subjs, sub_nums] = ...
    load_results(results_directory, whsim, dotest, LOG)

LOG.info('INFO', sprintf('Loading Results, whsim = %d', whsim));
if ismember(whsim, [26, 27])
    
    input_directory = fullfile(results_directory, 'batch_analyses', 'combined');
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
    
    % get model files
    ii = cellfun( @(x) regexp( x , 'allmodels' ), filenames, 'UniformOutput', 0);
    ii = cellfun( @(x) ~isempty(x), ii);
    filenames = filenames(ii);
    
    % parse model number from each filename
    mods = cellfun( @(x) str2num( x(regexp( x, '\d*'))), filenames, 'UniformOutput', 1);
    
    M = length(mods);
    models = cell(M,1);
    
    for i = 1:M
        m = mods(i);
        filename = fullfile( input_directory, sprintf( 'allmodels_%d.mat', m));
        load(filename);
        LOG.info('INFO', sprintf('\t%s', filename));
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
        'bad_subjs.mat'; ... 
        'sub_nums.mat'
        };
    
    for f = 1:length(filenames)
        filename = fullfile( input_directory, filenames{f});
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
    bad_subjs = [];
    all_subjs = [];
    
end

end