

whsim = 26; 
whmodel = 1; 
simfield = sprintf('sim%d', whsim); 

% load D 
load( fullfile( results_directory, '..', 'code', 'stats',...
            sprintf('cohensd_whs%d_whmodel%d.mat', whsim, whmodel)), 'D'); 
        
        
%% Load Gordon Atlas Labels 
filenm = fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', 'varianceGLM', ...
    'ROI2NIfTI', 'GordonParcels', 'Parcels.xlsx'); 
tab = readtable(filenm); 

%% Append Cohens d for each Contrast 
NC = length(D.(simfield).cohensd); 

for c = 1:NC
   tab = [tab, table( D.(simfield).cohensd{c}', 'VariableNames', { ['x' D.(simfield).contrast_names{c}] }) ];  
end

%% Sort table based on each contrasts 

cont_type = 'mean'; 
contrast = ['x2bck_minus_baseline' '_' cont_type]; 

[~, ii] = sort(abs(tab.(contrast)), 'descend'); 

tab2 = table( tab.ParcelID(ii), tab.Hem(ii), tab.SurfaceArea_mm2_(ii), tab.Community(ii), tab.(contrast)(ii)); 