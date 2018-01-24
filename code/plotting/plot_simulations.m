function [sims] = plot_simulations(whsims, isHPC, dotest)

% This function plots:
%   1) Parameter histograms
%   2) Brain visualization of model preference
%   3) Predicted vs actual
%   4) Area plots of model preference

% You can turn these plots on and off in the options section. 

% NOTE: Brain visualization plotting is not longer supported in this
% function

% INPUT:
%   array[numeric] whsims: which simulations do we want to plot?
%   bool isHPC: are we running on the HPC?
%   bool dotest: are we in test mode?

% OUTPUT:
%   struct sims: a structure containing the model run, proportion of data points
%       that prefer each model, design, and model descriptions for each
%       simulation

% The files needed to run this function:
%     '/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/code/batchmode/set_analysis_options_v2.m'
%     '/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/code/utils/tight_subplot.m'
%     '/Users/Garren/Dropbox/FMRI/PRojects/varianceGLM/code/plotting/plot_best_models.m'
%     '/Users/Garren/Dropbox/FMRI/PRojects/varianceGLM/code/plotting/plot_predicted_vs_actual.m'
%     '/Users/Garren/Dropbox/FMRI/PRojects/varianceGLM/code/plotting/plot_predicted_vs_actual_varmean_only.m'
%     '/Users/Garren/Dropbox/FMRI/PRojects/varianceGLM/code/plotting/plot_spatial_visualization'


addpath( fullfile(getenv('HOME'), 'Dropbox', 'MATLAButils'));

% set up logger
LOG = log4m.getLogger('crap.txt');
LOG.setCommandWindowLevel(LOG.INFO);
LOG.setLogLevel(LOG.OFF);

%% Options

plot_params = false;         % plot histograms of the parameters?
plot_brain = false;         % do brain plots?
plot_preds = true;         % do predicted vs actual plots
plot_best = true;           % do best model plots

var_method = 'meanpred';      % choose fom 'sample', 'meanpred'
nsamp = 100;                % numer of samples to use if var_method is 'sample'
models_to_plot = [1:4];     % which models to plot

[results_directory] = set_results_directory( isHPC ); 


%%
% list of simulations to plot figure for. Only the last will be currently displayed
Nsims =  length(whsims);

% Brain Net Viewer Options
opt_file = 'roi_draw_opts.mat';
surf_file = fullfile('..','..', 'BrainVisualization','BrainNetViewer','Data','SurfTemplate','BrainMesh_ICBM152.nv');
plot_mode = 'all';
metric = 'CV';

sims = struct();
for i = 1:Nsims
    whsim = whsims(i);
    
    % plot for a single simulation
    [propTab, models, design, descriptions, rst] = plot_one_sim(results_directory,...
        whsim, plot_params, opt_file, surf_file, plot_mode, metric,...
        plot_brain, isHPC, var_method, dotest, models_to_plot, LOG, plot_preds, plot_best);
    
    % store simulation output for group comparisons
    % have to add 'sim' to the fieldname because numbers aren't allowed to
    % begin fields of a structure
    sims.(['sim' num2str(whsim)]) = struct();
    sims.(['sim' num2str(whsim)]).propTab= propTab;
    sims.(['sim' num2str(whsim)]).models = models;
    sims.(['sim' num2str(whsim)]).design = design;
    sims.(['sim' num2str(whsim)]).descriptions = descriptions;
    sims.(['sim' num2str(whsim)]).rst = rst;
end

end


%%
function [propTab, models, design, descriptions, rst] = plot_one_sim(results_directory, ...
    whsim, plot_params, opt_file, surf_file, plot_mode, metric, ...
    plot_brain, isHPC, var_method, dotest, models_to_plot, LOG, plot_preds, plot_best)
% Function to plot the results from a single simulation

% INPUT:
%   char results_directory:
%   numeric whsim: which simulation results to plot
%   bool plot_params: do we want to plot histograms of parameter estimates?
%       Default plots are preferred models and model predictions versus
%       actual data.
%   char opt_file: options file for brain plots
%   char surf_file: surface file for brain plots
%   char plot_mode: something for brain plots
%   char metric: 'BIC' or 'CV'
%   bool plot_brain: do brain plots?
%   bool isHPC: are we running on HPC?
%   char var_method: for plotting model predictions, how to we want to compute
%       the variance (`sample` or `meanpred`)
%   bool dotest: testing?
%   array[numeric] models_to_plot: which models to plot
%   log4m LOG: logger

% OUTPUT:
%   tab propTab: lists the proportion of data points that prefer each
%       model
%   cellarray models: contains information about models run. E.g,
%         models{1}
%
%         ans =
%
%           struct with fields:
%
%                 code: {{1×4 cell}  {1×2 cell}  'Var+Mean+Fix+Ins'}
%                  meaneffects: {'Intercept'  '2-back'  'Instruction'  'Fixation'}
%                   vareffects: {'Intercept'  '2-back'}
%                  description: 'Var+Mean+Fix+Ins'
%                     meancols: [1 3 4 5]
%                      varcols: [1 3]
%                 designlabels: {'Intercept'  '2-back'  'Instruction'  'Fixation'}
%                        inits: [0 0.1000 0.1000 0.1000 1 0.1000]
%                  paramlabels: {6×1 cell}
%                        whsim: 2
%                  ismeanparam: [1 1 1 1 0 0]
%                   isvarparam: [0 0 0 0 1 1]
%             isinterceptparam: [1 0 0 0 1 0]
%                    allparams: [142×6×299 double]
%                    allpredsm: [395×142×299 single]
%                    allpredsv: [395×142×299 single]
%   design: the overall design matrix for that simulation (combination of
%       all columns used for all models)
%   descriptions: plain text description of each model run

%% Set Simulation
[ opts, dotest] = set_analysis_options_v2(whsim, isHPC, dotest, LOG);

%% Load the ROI timecourse data and Design
[ dat ] = load_data_and_design( opts, dotest, LOG );

%% Unpack Options and Data
% Options
whsim = opts.whsim;
whs = opts.whs;
K = opts.K;
seed = opts.seed;
Tremove = opts.Tremove;
doconstrained = opts.doconstrained;
prewhiten = opts.prewhiten;
var_log_transform = opts.var_log_transform;
TukN = opts.TukN;
multivariate = opts.multivariate;
roifile = opts.roifile;
designfile = opts.designfile;
combdesignfile = opts.combdesignfile;
runmodels = opts.runmodels;

% Data
tcn = dat.tcn;
design = dat.design;
designlabels = dat.designlabels;
T = dat.T;
NS = dat.NS;
R = dat.R;

%% Load Results

[models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, sub_nums] = ...
    load_results(results_directory, whsim, dotest, LOG);

NS = size( models{1}.allpredsm, 2 );
R  = size( models{1}.allpredsm, 3 ); 
M = length(models);

% get rid of subjects with parameters set to infinite values 
inf_params = false( NS, R ); 
for m = 1:M 
   inf_params = or( inf_params, squeeze( any( isinf( models{m}.allparams ), 2 ))); 
end

inf_params = sum( inf_params, 2 ) > 0; 
for m = 1:M 
   models{m}.allparams = models{m}.allparams( ~inf_params, :, : ); 
   models{m}.allpredsm = models{m}.allpredsm( :, ~inf_params, : ); 
   models{m}.allpredsv = models{m}.allpredsv( :, ~inf_params, : ); 
end 

allbicm( inf_params, :, : ) = []; 
allllsm( inf_params, :, : ) = []; 
bestmodelBIC( inf_params, : ) = []; 
bestmodelCV( inf_params, : ) = []; 
all_subjs = all_subjs( ~inf_params ); 
NS = sum( ~inf_params ); 


%%
% -----------------------------------------------
%      Plot Model Results
%  -----------------------------------------------
% NOTE: FOR plotting the HCP data, we are using the design matrix from the
% 1st subject. The design matrices are very similar, so for visualization,
% this should be ok.

% plotting options
lw = 1.5;

%% Convert model names in structure to cell array
descriptions = cell( M  , 1 );
for m=1:M
    descriptions{ m } = models{ m }.description{1};
    descriptions{ m } = strrep(descriptions{ m }, '+Fix+Ins', ''); % remove some of the HCP labeling for the paper
    descriptions{ m } = strrep(descriptions{ m }, '(prewhite)', ''); % remove '(pre-white)' for the paper
end


%%   1) Parameter histograms
if plot_params
    % plot parameters for each model m
    for m=1:M
        paramsnow   = models{ m }.allparams;
        paramlabels = models{ m }.paramlabels;
        
        % Reshape the parameters so it is in the same shape as the bestmodelBIC array
        P = size( paramsnow , 2 );
        params = permute( paramsnow , [ 1 3 2 ] );
        params = reshape( params , [ NS * R  P ] );
        
        % get data points that prefer model m
        wh = find( bestmodelBIC(:) == m );
        params = params( wh , : );
        
        h = figure( 70 + m ); clf;
        h.Position = [137 5 560 700];
        set( gcf , 'Name' , sprintf( 'Histogram %s' , models{ m }.description{1} ));
        for p=1:P
            ps = params( : , p );
            subplot( P , 1 , p );
            histogram( ps , 100 );
            str =  paramlabels{ p };
            title( str );
        end
    end
end

%% 2) Brain visualization of model preference
if plot_brain
    if ispc % Brain Net Viewer is screwing up on my MAC
        plot_spatial_visualization( whsim, roifile, input_directory,...
            plot_mode, metric, multivariate, whs, opt_file, surf_file)
    end
end

switch multivariate
    %% Mulativariate
    case 0
        if plot_preds
            %% 3) Predicted vs actual
            % all models CV
            fig = 299;
            name = 'Model Predictions Overall + Data (CV model Selection)';
            filenm = fullfile( results_directory, '..', 'images', sprintf('model_predictions_overall_CV_whsim%d', whsim));
            plot_predicted_vs_actual( tcn, fig, name, models, bestmodelCV, ...
                models_to_plot, design, var_method, filenm, sub_nums)
            
            
            
            % all models BIC
            fig = 298;
            name = 'Model Predictions Overall + Data (BIC model Selection)';
            filenm = fullfile( results_directory, '..', 'images', sprintf('model_predictions_overall_BIC_whsim%d', whsim));
            plot_predicted_vs_actual( tcn, fig, name, models, bestmodelBIC, ...
                models_to_plot, design, var_method, filenm, sub_nums)
            
            
            
            % just the Var+Mean and Var models CV
            fig = 300;
            groups = {[1,3]};
            name = 'Model Predictions Overall + Data (CV model Selection)';
            filenm = fullfile( results_directory, '..', 'images', sprintf('model_predictions_overall_CV_whsim%d_varmeanonly', whsim));
            plot_predicted_vs_actual_varmean_only( tcn, fig, name, models, bestmodelCV, design, filenm)
        end
end

%% 4) Area plots of model preference
if plot_best
    fig = 2;
    name =  'Best Models per ROI (BIC)' ;
    plotby = 'roi';
    filename = fullfile( results_directory, '..', 'images', sprintf('best_models_per_roi_BIC_whsim%d', whsim));
    plot_best_models( fig, models_to_plot, bestmodelBIC, name, descriptions, R, NS, plotby, filename)
    
    fig = 12;
    name = 'Best Models per ROI (CV)' ;
    filename = fullfile( results_directory, '..', 'images', sprintf('best_models_per_roi_CV_whsim%d', whsim));
    plot_best_models( fig, models_to_plot, bestmodelCV, name, descriptions, R, NS, plotby, filename)
    
    fig = 3;
    name = 'Best Models per Subject (BIC)';
    plotby = 'subj';
    filename = fullfile( results_directory, '..', 'images', sprintf('best_models_per_subject_BIC_whsim%d', whsim));
    plot_best_models( fig, models_to_plot, bestmodelBIC, name, descriptions, R, NS, plotby, filename)
    
    fig = 4;
    name = 'Best Models per Subject (CV)';
    filename = fullfile( results_directory, '..', 'images', sprintf('best_models_per_subject_CV_whsim%d', whsim));
    plot_best_models( fig, models_to_plot, bestmodelCV, name, descriptions, R, NS, plotby, filename)
end

%% Get percent of regions that preferred a model with a variance term

% get the number of data points that prefer each model
% -1 corresponds to time series of all zeros
biccount = histc(bestmodelBIC(:), [-1 1:length(descriptions)]);
cvcount  = histc(bestmodelCV(:),  [-1 1:length(descriptions)]);

% get proportions
bicprop = biccount/sum(biccount);
cvprop  = cvcount/sum(cvcount);

modellabels = ['NaN'; descriptions];

% create proportion table
propTab = table( bicprop, cvprop, modellabels, 'VariableNames', {'BIC', 'CV','Models'});

% get index into any model labels with variance in the title
ii = cell2mat(cellfun( @(x) ~isempty(strfind( x , 'Var')), modellabels, 'UniformOutput', 0));

% output proportion of regions that favor a model with 'Var' in the title
sum( propTab.BIC(ii))
sum(propTab.CV(ii))

% proportion of models favoring intercept
ii = strcmp( 'Intercept', modellabels);

propTab.BIC(ii)
propTab.CV(ii)

% Add rst if doesn't exist
if ~exist('rst', 'var')
    rst = NaN;
end
end
