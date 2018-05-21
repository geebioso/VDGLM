
whsim = 26;
simfield = sprintf('sim%d', whsim);
dotest = 0;
isHPC = 0;
whmodel = 1;
modelfield = sprintf('model%d', whmodel);

do_compute_tstats=1;

[results_directory] = set_results_directory( isHPC );

addpath('/Users/Garren/Dropbox/FMRI/BrainVisualization/NIfTI/NIfTI');
addpath('/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/ROI2NIfTI');
addpath('/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/ROI2NIfTI/dicm2nii');

%% Load Data

% [models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, sub_nums] = ...
%     load_results(results_directory, whsim, dotest, LOG, 'analyze');


[NS, P, R] = size( models{whmodel}.allparams);

% Bonferroni Correction 
significance_level = 0.05/R;

%% Define Contrasts
contrasts = struct();
% contrasts.sim7.model1 = {[ 0 1 0 -1 0 0 0], [ 0 0 0 0 0 0 1], [ 0 0 0 0 0 1 0]};
% contrasts.sim7.model2 = {[ 0 1 0 -1 0 ]};
%
% contrasts.sim1.model1 = {[ 0 1 -1 0 0 0], [ 0 0 0 0 0 1]};
% contrasts.sim1.model2 = {[ 0 1 -1 0 0]};

contrasts.sim7.model1 = {
    [ -1 1 0 0 0 0 0], ...
    [ 0 0 0 0 -1 1 0], ...
    [ 0 1 0 0 0 0 0], ...
    [ 0 0 0 0 0 1 0], ...
    [ 1 0 0 0 0 0 0], ...
    [ 0 0 0 0 1 0 0], ...
    [ 0 0 1 0 0 0 0], ...
    [ 0 0 0 0 0 0 1]
    };

contrasts.sim7.names = {
    '2bck_minus_0bck_mean', ...
    '2bck_minus_0bck_var',...
    '2bck_minus_baseline_mean',...
    '2bck_minus_baseline_var', ...
    '0bck_minus_baseline_mean', ...
    '0bck_minus_baseline_var', ...
    'Instruction_minus_baseline_mean',...
    'Instruction_minus_baseline_var'
    };

contrasts.sim1.model1 = {[ 0 1 -1 0 0 0], [ 0 0 0 0 0 1]};
contrasts.sim1.names= {'2bck_minus_0bck_mean', '2bck_minus_baseline_variance'};

% same contrasts for pre-whitened and CIfTI analyses of the HCP data
contrasts.sim26 = contrasts.sim7;
contrasts.sim27 = contrasts.sim7;
contrasts.sim6 = contrasts.sim7;




%% Compute t-tests
allparams = models{whmodel}.allparams;
NC = length( contrasts.(simfield).names);

if do_compute_tstats
    TSTAT  = struct();
    
    TSTAT.(simfield).contrast_names = cell(NC,1);
    TSTAT.(simfield).pvals = cell(NC,1);
    TSTAT.(simfield).tstats = cell(NC,1);
    
    for c = 1:NC
        
        C = contrasts.(simfield).(modelfield){c};
        contrast_str = contrasts.(simfield).names{c};
        contrast_str = strrep(contrast_str,' ', '');
        
        % compute individual contrasts
        C2 = zeros(1,1,P);
        C2(:) = C;
        copes = sum( bsxfun( @times, permute(allparams, [1 3 2]), C2), 3); % checked that this is the same result as using a for loop
        
        temp_pvals = zeros(1, R);
        temp_tstats = zeros(1, R);
        
        for r = 1:R
                [hyp,pval,ci,stats] = ttest(copes(:,r) );
                temp_tstats(r) = stats.tstat;
                temp_pvals(r) = pval;
        end
        
        df = stats.df; 
        TSTAT.(simfield).pvals{c} = temp_pvals;
        TSTAT.(simfield).tstats{c} = temp_tstats;
        TSTAT.(simfield).contrast_names{c} = contrasts.(simfield).names{c}; 
        
    end
    
    %% Compute t-quantile 
    t_threshold = tinv( 1 - significance_level, df);
    
    filename = fullfile( results_directory, 'single_analyses', ...
        sprintf('balanced_tstats_whs%d_whmodel%d', whsim, whmodel));
    
    save(filename, 'TSTAT'); 
    
    
    %% Output NIfTI images 
    
       
    filename = fullfile(results_directory, '../ROI2NIfTI', 'files', ...
        sprintf('balanced_tstats_whs%d_whmodel%d', whsim, whmodel));
    
    temp = zeros(R, NC);
    for c = 1:NC
        temp(:,c) = TSTAT.(simfield).tstats{c};
    end
    
    ROI2dscalar_nii_v2(temp, filename, TSTAT.(simfield).contrast_names, 'pct');
    
end

