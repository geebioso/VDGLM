

addpath(fullfile(getenv('HOME'), 'Dropbox', 'MATLAButils'));
addpath('..'); 

input_directory = 'Perm_Test_Results_all'; 
whsim = 26; 
dotest = 0; 

LOG = log4m.getLogger('crap.txt');
LOG.setCommandWindowLevel(LOG.INFO);
LOG.setLogLevel(LOG.OFF);

%% Load Data 
[models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, bad_subjs] = load_results(input_directory, whsim, dotest, LOG); 

[NS, R, Nsamp, Pv] = size(models{1}); 
Pm = size(models{2}, 4); 

models_to_plot = 1:2; 
descriptions = {'Var+Mean', 'Mean'}; 

%% Compute Error Bars at 95% Level 
bic = struct(); 
bic.sub_vm  = squeeze( mean( bestmodelBIC == 1, 1) ); 
bic.sublow  = quantile( bic.sub_vm, 0.025, 2); 
bic.subhigh = quantile( bic.sub_vm, 0.975, 2); 

bic.roi_vm  = squeeze( mean( bestmodelBIC == 1, 2)); 
bic.roilow  = quantile( bic.roi_vm, 0.025, 2); 
bic.roihigh = quantile( bic.roi_vm, 0.975, 2); 

cv = struct(); 
cv.sub_vm  = squeeze( mean( bestmodelCV == 1, 1) ); 
cv.sublow  = quantile( cv.sub_vm, 0.025, 2); 
cv.subhigh = quantile( cv.sub_vm, 0.975, 2); 

cv.roi_vm  = squeeze( mean( bestmodelCV == 1, 2)); 
cv.roilow  = quantile( cv.roi_vm, 0.025, 2); 
cv.roihigh = quantile( cv.roi_vm, 0.975, 2); 

%% Load data 

% input_directory = 'Results_all'; 
% whsim = 26; 
% 
% [models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, bad_subjs] = load_results(input_directory, whsim, dotest, LOG); 

%% Determine Best Model for Subject + ROI (using BIC)
% Area plots across subjects
figure( 2 ); clf;
py = zeros( R , length(models_to_plot) );
for j=1:R
    ps = hist( bestmodelBIC(:,j) , models_to_plot );
    py( j , : ) = ps * 100 / sum( ps );
end
% Rearranging ROIs according to the areas best fit by first model
[ ~ , index ] = sort( py(:,1) , 1 , 'descend' );
area( py( index , : ) );
legend( descriptions );
xlabel( 'ROI' );
ylabel( '% Subjects' );
ylim( [ 0 100 ] );
xlim( [ 1 R+0.5 ] );
title( sprintf('Best Model for each ROI and Subject (using BIC, pw=%d)',prewhiten));
set( gcf , 'Name' , 'Best Models per ROI (BIC)' );
set( gca, 'FontSize', 13);

%% Determine Best Model for Subject + ROI (CV)
% Area plots across subjects
figure( 12 ); clf;
py = zeros( R , length(models_to_plot) );
for j=1:R
    ps = hist( bestmodelCV(:,j) , models_to_plot );
    py( j , : ) = ps * 100 / sum( ps );
end
% Rearranging ROIs according to the areas best fit by first model
[ ~ , index ] = sort( py(:,1) , 1 , 'descend' );
area( py( index , : ) );
legend( descriptions );
xlabel( 'ROI' );
ylabel( '% Subjects' );
ylim( [ 0 100 ] );
xlim( [ 1 R+0.5 ] );
title( sprintf('Best Model for each ROI and Subject (using CV, pw=%d)', prewhiten));
set( gcf , 'Name' , 'Best Models per ROI (CV)' );
set( gca, 'FontSize', 13);

%% Determine Best Model for Subject + ROI (using BIC)
% Area plots across subjects
figure( 3 ); clf;
py = zeros( NS , length(models_to_plot) );
for j=1:NS
    ps = hist( bestmodelBIC(j,:) , models_to_plot );
    py( j , : ) = ps * 100 / sum( ps );
end
% Rearranging ROIs according to the areas best fit by first model
[ ~ , index ] = sort( py(:,1) , 1 , 'descend' );
area( py( index , : ) );
legend( descriptions );
xlabel( 'Subject' );
ylabel( '% ROIs' );
ylim( [ 0 100 ] );
xlim( [ 1 NS ] );
%xlim( [ 1 NS+0.5 ] );
% title( sprintf('Best Model for each ROI and Subject (using BIC)\nhp=%d, pw=%d', ...
%     dohighpassfilter, prewhiten));
set( gcf , 'Name' , 'Best Models per Subject (BIC)' );
set( gca, 'FontSize', 13);

filename = fullfile('images', sprintf('perm_test_best_models_per_subject_BIC_whsim%d', whsim));
print( filename, '-depsc');

%% Determine Best Model for Subject + ROI (using CV)
% Area plots across subjects
figure( 4 ); clf;
py = zeros( NS , length(models_to_plot) );
for j=1:NS
    ps = hist( bestmodelCV(j,:) , models_to_plot );
    py( j , : ) = ps * 100 / sum( ps );
end
% Rearranging ROIs according to the areas best fit by first model
[ ~ , index ] = sort( py(:,1) , 1 , 'descend' );
area( py( index , : ) );
legend( descriptions );
xlabel( 'Subject' );
ylabel( '% ROIs' );
ylim( [ 0 100 ] );
xlim( [ 1 NS ] );
%xlim( [ 1 NS+0.5 ] );
% title( sprintf('Best Model for each ROI and Subject (using CV)\nhp=%d, pw=%d', ...
%     dohighpassfilter, prewhiten ));
set( gcf , 'Name' , 'Best Models per Subject (CV)' );
set( gca, 'FontSize', 13);

filename = fullfile('images', sprintf('perm_test_best_models_per_subject_CV_whsim%d', whsim));
print( filename, '-depsc');