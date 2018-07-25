% function [] = plot_design(whsim)

whsim  = 26;
isHPC  = 0;
dotest = 0;
set_up_directory_structure = 1;

%To create the logger reference:
LOG = log4m.getLogger('crap.txt');
LOG.setCommandWindowLevel(LOG.('INFO'));
LOG.setLogLevel(LOG.OFF);

%% Get directories 

[results_directory, images_directory, ROI2NIfTI_directory] = set_results_directory( isHPC, set_up_directory_structure );

%% Set Simulation
[ opts, dotest] = set_analysis_options(whsim, isHPC, dotest, set_up_directory_structure, LOG);

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
output_directory = opts.output_directory;
addmotion = opts.addmotion;

% Data
tcn = dat.tcn;
design = dat.design;
designlabels = dat.designlabels;
T = dat.T;
R = dat.R;
motionX = dat.motionX;
scrubX = dat.scrubX;

LOG.info('INFO', sprintf('\tdotest = %d', dotest));
LOG.debug('DEBUG', sprintf('\tR = %d', R));
LOG.debug('DEBUG', sprintf('\tT = %d', T));

fprintf('Plotting Design for Analysis %d\n', whsim);

%% Plot Design

cols_to_plot = {'0-back', '2-back', 'Instruction', 'Fixation'}; 
ii = ismember( designlabels, cols_to_plot); 

designlabels = designlabels(ii); 
design = design(:,ii); 

T = size(design, 1);

figure( 78 ); clf;
for i=1:size(design,2)
    plot( design(:,i) *0.7 + i + 1 ); hold on;
end
set( gca , 'YTick' , (2:(size(design,2)+1))  );
set( gca , 'YTickLabel' ,  designlabels);
set( gca , 'TickLabelInterpreter' , 'none' );
xlim( [ 0 T ] );
%xlabel('time'); 
set( gca, 'XTick', []); 

title( 'HCP WM Experimental Design' ); 

filename = fullfile(images_directory, 'HCP_WM_design');
print(filename, '-depsc'); 
