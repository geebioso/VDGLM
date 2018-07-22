

% this script computes the color bounds for our contrasts with the
% following workbench parameter settings: 

posmin=0.04;
posmax=0.96;
negmin=0.02;
negmax=0.98;

% These values can either refer to fixed values (across contrasts) or
% percentage thresholds. 

dofixed = 1; % should we compute bounds for fixed bounds (across contrasts), or to have a separate bound for each contrast (0/1)

addpath('/Users/Garren/Dropbox/FMRI/BrainVisualization/NIfTI/NIfTI');
addpath('/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/ROI2NIfTI/dicm2nii');

% Add simulation information 
whsim = 26; 
whmodel = 2; 
isHPC = 0; 
dotstats = 0; 

[results_directory] = set_results_directory( isHPC, set_up_directory_structure);

%% Load data 

load( fullfile( results_directory, '..', 'code', 'stats',...
            sprintf('cohensd_whs%d_whmodel%d.mat', whsim, whmodel))); 
        
NC = length(D.sim26.cohensd); 

%% Compute quantiles 
p = 0; 

ismean = cellfun( @(x) isempty( strfind( x, 'var') ), D.sim26.contrast_names); 
isvar = ~ismean; 


all_mean_cohensd = cell2mat(D.sim26.cohensd(ismean));
all_var_cohensd = cell2mat(D.sim26.cohensd(isvar));

varmin = quantile( all_var_cohensd(:), 0); 
varmax = quantile( all_var_cohensd(:), 1);

meanmin = quantile( all_mean_cohensd(:), p); 
meanmax = quantile( all_mean_cohensd(:), 1- p); 

f = fopen( fullfile( '..','plotting', 'color_bounds.txt'), 'w'); 

fprintf( f, '%2.2f\n', meanmin ); 
fprintf( f, '%2.2f\n', meanmax ); 
fprintf( f, '%2.2f\n', varmin ); 
fprintf( f, '%2.2f', varmax ); 

fclose(f); 

%% Get number of voxels per ROI 
filepath = fullfile( getenv( 'HOME'), 'Dropbox', 'FMRI', 'Projects', 'varianceGLM', 'ROI2NIfTI'); 
roi = nii_tool('img', [filepath '/Gordon333.32k_fs_LR.dlabel.nii']);
roi = squeeze(roi);
roi = roi(1:59412);

roi_voxel_count = zeros(333,1); 
for i = 1:333
    roi_voxel_count(i)= sum( roi==i );
    
end


% quantile( cohensd_expanded, 0.98)

%% Compute Cohen's d MATLAB colors by linear interpolation

map = fsl_colormap(); 

close all; 

D.sim26.bounds = zeros(NC, 2); 
for c = 1:NC
    
    contrast_now = D.sim26.contrast_names{c}; 
    cohensd_now = D.sim26.cohensd{c}; 
    
    % Expand Cohen's D (workbench uses voxels, not ROIs to compute
    % quantiles) 
    cohensd_expanded = cell(333,1); 
    for i = 1:333
        cohensd_expanded{i} = cohensd_now(i)*ones(roi_voxel_count(i), 1); 
    end
    cohensd_expanded = cell2mat( cohensd_expanded ) ; 

    posmax = quantile( cohensd_expanded( cohensd_expanded > 0 ), 0.96);
    negmin = quantile( cohensd_expanded( cohensd_expanded < 0 ), 0.02);
    
    ii = and( cohensd_now > negmin, cohensd_now < posmax ); 
    cohensd_now = cohensd_now(ii); 
    cohensd_now(1) = posmax; 
    cohensd_now(2) = negmin; 
    
    D.sim26.bounds(c,:) = [negmin, posmax]; 
    
    figure(c); clf; imagesc(cohensd_now); h = colorbar; colormap(map); title( strrep( contrast_now, '_', ' ') );
    delta = (posmax - negmin);
    h.Ticks = [negmin + 0.03*delta delta/2 + negmin posmax - 0.03*delta]; 
    h.TickLabels = {num2str(round(negmin, 2)), '±0.2',  num2str(round(posmax, 2))}; 
    h.TickLength = 0; 
    h.FontSize = 20;
    
    
    filename = fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', ...
        'varianceGLM', 'images', sprintf( 'colorbar_cohensd_mdl%d_%s', whmodel, contrast_now));
    print(filename, '-depsc');


    if strcmp(contrast_now, '0bck_minus_baseline_mean')
        x = 1;
    end
end


%% Compute Tstat MATLAB colors by linear interpolation

if dotstats
filename = fullfile( results_directory, 'single_analyses', ...
        sprintf('balanced_tstats_whs%d_whmodel%d', whsim, whmodel));
    
load(filename);  

map = fsl_colormap(); 
close all; 

TSTAT.sim26.bounds = zeros(NC, 2); 
for c = 1:NC
    
    contrast_now = TSTAT.sim26.contrast_names{c}; 
    tstat_now = TSTAT.sim26.tstats{c}; 
    
    % Expand Cohen's D (workbench uses voxels, not ROIs to compute
    % quantiles) 
    tstat_expanded = cell(333,1); 
    for i = 1:333
        tstat_expanded{i} = tstat_now(i)*ones(roi_voxel_count(i), 1); 
    end
    tstat_expanded = cell2mat( tstat_expanded ) ; 

    posmax = quantile( tstat_expanded( tstat_expanded > 0 ), 0.96);
    negmin = quantile( tstat_expanded( tstat_expanded < 0 ), 0.02);
    
    ii = and( tstat_now > negmin, tstat_now < posmax ); 
    tstat_now = tstat_now(ii); 
    tstat_now(1) = posmax; 
    tstat_now(2) = negmin; 
    
    D.sim26.bounds(c,:) = [negmin, posmax]; 
    
    figure(c+10); clf; imagesc(tstat_now); h = colorbar; colormap(map); title( strrep( contrast_now, '_', ' ') );
    delta = (posmax - negmin);
    h.Ticks = [negmin + 0.03*delta delta/2 + negmin posmax - 0.03*delta]; 
    h.TickLabels = {num2str(round(negmin, 2)), '±0.2',  num2str(round(posmax, 2))}; 
    h.TickLength = 0; 
    h.FontSize = 20;
    
    
    filename = fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', ...
        'varianceGLM', 'images', sprintf( 'colorbar_tstat_%s', contrast_now));
    print(filename, '-depsc', '-opengl', '-r0');


    if strcmp(contrast_now, '0bck_minus_baseline_mean')
        x = 1;
    end
end

%% Compute HCP colorbars by linear interpolation 

addpath( fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'cifti-matlab')); 
addpath( fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'cifti-matlab', '@xmltree')); 

filename = fullfile( results_directory, '..', 'ROI2NIfTI', ...
    'HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMSulc.dscalar.nii');
    
tmpfile = fullfile( results_directory, '..', 'ROI2NIfTI', ...
    'HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMSulc.gii');

cifti = gifti(tmpfile);

% get bounds 

% Define maps 
cdata_now = cifti.cdata(:, 9:11); 
posmax = max(max(cdata_now(:))); 
negmin = min(min(cdata_now(:))); 

filenm =  fullfile('..', 'plotting', sprintf('bounds_hcp.txt', whsim, whmodel));
f = fopen(filenm, 'w');
fprintf(f, '%2.2f\n%2.2f\n', negmin, posmax);
fclose(f);

% print fixed colorbar 
figure(1); clf; imagesc(cdata_now(:)); h = colorbar; colormap(map); title( strrep( contrast_now, '_', ' ') );
delta = (posmax - negmin);
h.Ticks = [negmin + 0.03*delta delta/2 + negmin posmax - 0.03*delta];
h.TickLabels = {num2str(round(negmin, 2)), '±0.2',  num2str(round(posmax, 2))};
h.TickLength = 0;
h.FontSize = 20;

 filename = fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', ...
        'varianceGLM', 'images', sprintf( 'colorbar_hcp_fixed_mean', contrast_now));
    print(filename, '-dpng'); 
    
    
    

cont = struct(); 
cont.idx = [9,10,11]; 
cont.names = {'2-back_minus_Fixation', '0-back_minus_Fixation', '2-back_minus_0-back'}; 

% 

map = fsl_colormap(); 

for i = 1:3
    
    idx = cont.idx(i); 
    contrast_now = cont.names{i}; 
    cdata_now = cifti.cdata(:,idx); 
    
    
    posmax = quantile( cdata_now( cdata_now > 0 ), 0.96);
    negmin = quantile( cdata_now( cdata_now < 0 ), 0.02);
    
    ii = and( cdata_now > negmin, cdata_now < posmax ); 
    cdata_now = cdata_now(ii); 
    cdata_now(1) = posmax; 
    cdata_now(2) = negmin; 
    
    
    figure(i); clf; imagesc(cdata_now); h = colorbar; colormap(map); title( strrep( contrast_now, '_', ' ') );
    delta = (posmax - negmin);
    h.Ticks = [negmin + 0.03*delta delta/2 + negmin posmax - 0.03*delta]; 
    h.TickLabels = {num2str(round(negmin, 2)), '±0.2',  num2str(round(posmax, 2))}; 
    h.TickLength = 0; 
    h.FontSize = 20;
    
    
    filename = fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', ...
        'varianceGLM', 'images', sprintf( 'colorbar_hcp_%s', contrast_now));
    print(filename, '-depsc', '-opengl', '-r0');
    
end

end

%% Color bounds for fixed color scale 
if dofixed 
    
    map = fsl_colormap();
    
    close all;
    
    % save fixed min/max bounds for mean and variance contrats
    ismean = cellfun( @(x) contains( x, 'mean'), D.sim26.contrast_names);
    isvar = ~ismean;
    is_inst = cellfun( @(x) contains( x, 'Instruction'), D.sim26.contrast_names);
    
    filenm =  fullfile('..', 'plotting', sprintf('bounds_whs%d_mdl%d_mean.txt', whsim, whmodel)); 
    mean_bounds = readtable( filenm, 'ReadVariableNames', 0);
    mean_bounds = mean_bounds.Var1;
    posmax = mean_bounds(2); 
    negmin = mean_bounds(1); 
    
    idx = and(ismean, ~is_inst);
    figure(1); clf; imagesc([D.sim26.cohensd{idx}]); h = colorbar; colormap(map); title( strrep( contrast_now, '_', ' ') );
    delta = (posmax - negmin);
    h.Ticks = [negmin + 0.03*delta delta/2 + negmin posmax - 0.03*delta]; 
    h.TickLabels = {num2str(round(negmin, 2)), '±0.2',  num2str(round(posmax, 2))}; 
    h.TickLength = 0; 
    h.FontSize = 20;
    
    filename = fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', ...
        'varianceGLM', 'images', sprintf( 'colorbar_cohensd_fixed_mdl%d_mean', whmodel));
    print(filename, '-dpng');
    
    filenm =  fullfile('..', 'plotting', sprintf('bounds_whs%d_mdl%d_var.txt', whsim, whmodel)); 
    var_bounds = readtable( filenm, 'ReadVariableNames', 0);
    var_bounds = var_bounds.Var1;
    posmax = var_bounds(2); 
    negmin = var_bounds(1); 
    
    if whmodel == 1
        idx = and(isvar, ~is_inst);
        figure(2); clf; imagesc([D.sim26.cohensd{idx}]); h = colorbar; colormap(map); title( strrep( contrast_now, '_', ' ') );
        delta = (posmax - negmin);
        h.Ticks = [negmin + 0.03*delta delta/2 + negmin posmax - 0.03*delta];
        h.TickLabels = {num2str(round(negmin, 2)), '±0.2',  num2str(round(posmax, 2))};
        h.TickLength = 0;
        h.FontSize = 20;
        
        
        filename = fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', ...
            'varianceGLM', 'images', sprintf( 'colorbar_cohensd_fixed_mdl%d_var', whmodel));
        print(filename, '-dpng');
    end
    
    
    
end 


