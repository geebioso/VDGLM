% function [] = plot_null_best(whsim, isHPC, dotest)

% This function plots:
%   1) Area plots of model preference

% You can turn these plots on and off in the options section.

% NOTE: Brain visualization plotting is not longer supported in this
% function

% INPUT:
%   array[numeric] whsims: which simulations do we want to plot?
%   bool isHPC: are we running on the HPC?
%   bool dotest: are we in test mode?
%   char results_directory:
%   numeric whsim: which simulation results to plot
%   bool plot_params: do we want to plot histograms of parameter estimates?
%       Default plots are preferred models and model predictions versus
%       actual data.
%   char metric: 'BIC' or 'CV'
%   bool plot_brain: do brain plots?
%   bool isHPC: are we running on HPC?
%   char var_method: for plotting model predictions, how to we want to compute
%       the variance (`sample` or `meanpred`)
%   bool dotest: testing?
%   array[numeric] models_to_plot: which models to plot
%   log4m LOG: logger

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

models_to_plot = 1:2;       % which models to plots
[results_directory] = set_results_directory( isHPC );

%%
% list of simulations to plot figure for. Only the last will be currently displayed
% Brain Net Viewer Options
plot_mode = 'all';
metric = 'CV';

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

load(fullfile( results_directory, 'batch_analyses', 'combined', sprintf('whs%d_allllsm.mat', whsim))); 
bestmodelCV  = NaN( NS , R );
for s=1:NS
    for j=1:R
            [ a , whmax ] = max( allllsm( s,j,1:2) );
            if isnan(a)
                bestmodelCV( s , j ) = -1; % set NaN to modelid -1
            else
                bestmodelCV( s , j ) = whmax;
            end
    end
end

x = load(fullfile( results_directory, 'batch_analyses', 'null_combined', sprintf('whs%d_bestmodelCV.mat', whsim))); 
bestmodelCV_null = x.bestmodelCV; 

NS = size( models{1}.allparams, 1 );
R  = size( models{1}.allparams, 2 );
M = length(models);

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

%% 4) Area plots of model preference

fig = 4;
name = 'Best Models per Subject (CV)';
filename = fullfile( results_directory, '..', 'images', sprintf('best_models_per_subject_CV_whsim%d_null_sample', whsim));
%% Create Area Plot 

figure( fig ); clf;
lw = 1.2; 
py = zeros( NS , length(models_to_plot) );
py_null = zeros( NS , length(models_to_plot) );
for j=1:NS
    ps = hist( bestmodelCV(j,:) , models_to_plot );
    py( j , : ) = ps * 100 / sum( ps );
    
    ps_null = hist( bestmodelCV_null(j,:) , models_to_plot );
    py_null( j , : ) = ps_null * 100 / sum( ps_null );
end
% Rearranging ROIs according to the areas best fit by first model
[ ~ , index ] = sort( py(:,1) , 1 , 'descend' );
plot( py( index , 1 ) , 'LineWidth', lw); hold on; 
plot( py_null(index, 1), 'LineWidth', lw); hold on; 
% plot( py( index , 2 ) ); hold on; 
% plot( py_null(index, 2)); hold on; 

legend( descriptions(models_to_plot) );
xlabel( 'Subject' );
ylabel( '% ROIs' );
ylim( [ 0 100 ] );
xlim( [ 1 NS ] );

legend( 'Location', 'best', {'Real data', 'Null samples'}); % , 'Real M', 'Null M'}); 

hold on ; 
title( sprintf('Percent of ROIs that prefer the VDGLM model'));
set( gcf , 'Name' , name );
set( gca, 'FontSize', 13);

print( filename, '-depsc');

%% Load Null Results 
[models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, sub_nums] = ...
    load_results(results_directory, whsim, dotest, LOG, 'null');


