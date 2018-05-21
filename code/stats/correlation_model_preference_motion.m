

whsim = 26; 
dotest = 0; 
isHPC = 0; 

LOG = log4m.getLogger('crap.txt');
LOG.setCommandWindowLevel(LOG.INFO);
LOG.setLogLevel(LOG.OFF);

[results_directory] = set_results_directory( isHPC );

%% Set Simulation
[ opts, dotest] = set_analysis_options_v2(whsim, isHPC, dotest, LOG);

%% Load the ROI timecourse data and Design
[ dat ] = load_data_and_design( opts, dotest, LOG );

%% Unpack Options and Data

% Data
tcn = dat.tcn;
design = dat.design;
designlabels = dat.designlabels;
T = dat.T;
R = dat.R;
NS = size(tcn, 2); 
motionX = dat.motionX;
scrubX = dat.scrubX;

%% Load Results 
load(fullfile( results_directory, 'batch_analyses', 'combined', sprintf('whs%d_allllsm.mat', whsim))); % allllsm: OOSLL values 
load(fullfile( results_directory, 'batch_analyses', 'combined', sprintf('whs%d_all_subjs.mat', whsim))); % all_subjs: list of subjects 
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
pct_ROI = mean( bestmodelCV == 1, 2);

%% Compute correlation between ROI preference and motion
meanFD = zeros(NS,1); 

for s = 1:NS

    b = motionX{s}; 

    FD = [zeros(1,6); diff(b)];
    FD(:,1:3) = FD(:,1:3) * 50; % 50 mm is radius of brain
    FD = sum(abs(FD), 2); % Power et al method

    meanFD(s) = mean(FD); 
end

%% 
h = figure(1); clf; 
h.Position =  [34 285 584 420]; 
mdl = fitlm( meanFD, pct_ROI); 
plot( mdl ); 

xlabel( 'Mean FD'); 
ylabel( '%% ROI prefers V+M'); 

title( 'Head motion vs V+M preference'); 

filename = fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', 'varianceGLM', ...
    'images', sprintf('head_motion_vs_VM_preference_whs%d', whsim)); 

print( filename, '-dpng'); 

%% See if GMIN correlates with model performance 

gmin = load(fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'restingstatedata', 'WM_LR_gmin.mat')); 
gmin_subjs = gmin.subjs; 
gmin = gmin.gmin; 

% index from gmin to allsubjs
ii = ismember( all_subjs, gmin_subjs); 
roi_pref = mean( bestmodelCV, 2); 
roi_pref = roi_pref(ii); 
    
if isequal( all_subjs(ii), gmin_subjs ) 

    mdl = fitlm( gmin, roi_pref ); 
    plot( mdl ); 
end







