

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

%% Get baseline community percentages 
% get number of community types 
community_names = unique(tab.Community); 
com_num = zeros( length(community_names), 1); 
for c = 1:length(community_names)
    com_num(c) = sum( ismember( tab.Community, community_names{c})); 
end

% create table of percent of community types 
com_tab = table(com_num, community_names, 'VariableNames', {'num', 'name'});
[~, ii] = sort(com_num, 'descend');
com_tab = com_tab(ii,:); 

%% Where Effects Occur 
threshold = 0.5;

% cont_type = 'var'; 
% contrast = ['x2bck_minus_baseline' '_' cont_type]; 

cont_type = 'mean'; 
contrast = ['x2bck_minus_baseline' '_' cont_type]; 

% get all regions with absolute cohen's d above threshold 
ii = abs(tab.(contrast)) >= threshold; 
tab2 = table( tab.ParcelID(ii), tab.Hem(ii), tab.SurfaceArea_mm2_(ii), ...
    tab.Community(ii), tab.(contrast)(ii), ...
    'VariableNames', {'ParcelID', 'Hem', 'SurfaceArea_mm2_', 'Community', contrast}); 
    Nabove_thresh = sum(ii); 
    
% get number of community types above threshold 
community_names = unique(tab2.Community); 
prec = zeros( length(community_names), 1); 
tpr  = zeros( length(community_names), 1); 
for c = 1:length(community_names)
    total_c = com_tab.num( strcmp( com_tab.name , community_names{c})); 
    prec(c) = sum( ismember( tab2.Community, community_names{c}))/total_c; 
    tpr(c)  = sum( ismember( tab2.Community, community_names{c}))/Nabove_thresh;
end

% get number 
com_tab2 = table(community_names, prec, tpr, prec.*tpr, 'VariableNames', {'name', 'prec', 'tpr', 'prec_tpr'});
[~, ii] = sort(com_tab2.prec_tpr, 'descend');
com_tab2 = com_tab2(ii,:); 


