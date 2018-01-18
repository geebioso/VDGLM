function [] = plot_spatial_visualization( whsim, roifile, input_directory,...
    plot_mode, metric, multivariate, whs, opt_file, surf_file)

% Var+Mean: red
% Mean: Orange
% Var: 
% Transitions: Light Green  
% Intercept: Green

% this function will visualize the preferred model for each ROI

% INPUT:
%   whsims: which simulations do we want to plot?
%   roifile: the nii file that gives us ROI locations
%   input_directory: where are model results stored? 
%   bv_directory: where is brain viewer?
%   mode: which mode do we want 'VDGLMvGLM' or 'all'?
%       - VDGLMvGLM: plots the ROIs for which comparison metric is higher
%       for VDGLM versus GLM
%       - all: plots the models with the highest model comparison score
%   metric: choice of metric 'BIC' or 'CV'

bv_directory = fullfile('..','..','BrainNetViewer_20170403');
nifti_directory = fullfile('..','..','BrainVisualization','NIfTI','NIfTI');

addpath(bv_directory); 
addpath(nifti_directory); 

%% Get Preferred Regions 
roifile = fullfile('..', '..', 'BrainVisualization','ROI_299_3mm.nii.gz'); 
rois = load_nii(roifile); 

filenm = fullfile(input_directory, sprintf('paramswm_cc_whs%d.mat', whsim)); 
load(filenm); 
M = length(models); 

% Load 269 index into 299 
switch whs
    case 8
        old_roi = load( fullfile( '..', '..', '..', 'anatomical', 'ICAresults', 'tc_csf_wm_motion_out_globalMask5000_203subj.mat'), 'roi');
        old_roi = old_roi.roi;
        un_old_roi = unique(old_roi);
        un_old_roi( un_old_roi == 0 ) = []; % we don't care about 0 for our ROIs

        unique_rois = un_old_roi; 
    case 0
        unique_rois = unique(rois.img(:)); 
        unique_rois( unique_rois == 0 ) = []; % don't care about 0 for our ROIs 
end 

R = length(unique_rois); 

%% Get Preferred Models 
% take the mode for univariate analysis 
switch multivariate 
    case 0
        if strcmp(metric, 'BIC')
            roi_best = mode(bestmodelBIC);
        elseif strcmp(metric, 'CV')
            roi_best = mode(bestmodelCV);
        end
    case 1
        if strcmp(metric, 'BIC')
            roi_best = bestmodelBIC;
        elseif strcmp(metric, 'CV')
            roi_best = bestmodelCV;
        end
end 

%% Compute the Preferred Models Image 
preferred = rois; 
preferred.img( ~ismember(preferred.img(:), unique_rois) ) = 0; 
for r = 1:R
    nmatch = sum( preferred.img(:) == unique_rois(r)); 
    fprintf('ROI %d: nmatch=%d\n', unique_rois(r), nmatch); 
    preferred.img( preferred.img == unique_rois(r) ) = roi_best(r); 
end


% view_nii(preferred); 
save_file = fullfile( input_directory, sprintf('preferred_rois_whs%d.nii', whsim)); 
save_nii( preferred, save_file); 

%% Brain Net: Surface, volume, options, output
save_filenm = fullfile('images', sprintf('preferred_rois_whs%d.eps', whsim) );
colormap(jet(M)); 
h = BrainNet_MapCfg( ...
    surf_file,...
    save_file,...
    opt_file, ...
    save_filenm);

axes(h.Children(end))
hb = colorbar('southoutside'); 


hb.Position(2) = hb.Position(2) - 0.01; 
hb.Position(4) = hb.Position(4) - 0.02; 

%caxis([ 1 M ])
hb.Ticks = 1:M; 
hb.TickLabels = {'V+M', 'M', 'V', 'T', 'I'}; 

CD = h.Children(3).Children(4);

end


