function [] = plot_cv_scores_whsim26(isHPC, dotest)

% This function plots the out of sample log likelihood scores 

%% Set directory 
[results_directory] = set_results_directory( isHPC ); 

input_directory = fullfile( results_directory, 'batch_analyses', 'combined');
if dotest
    input_directory = [input_directory '_test'];
end

fprintf('Input directory: %s\n', input_directory);

%% Load Results 
load(fullfile( input_directory, 'allllsm.mat')); 
load(fullfile( input_directory, 'bestmodelCV.mat')); 

% sort subjects by how many ROIS prefer the Var+Mean model 
pref_vm = bestmodelCV == 1; 
mean_pref_vm = mean(pref_vm,2); 
[pct_pref, ii] = sort(mean_pref_vm);

% get subject with lowest, middle, and highest percent of regions that
% prefer the Var+Mean model 
subs = [ii(1), ii(round(333/2)), ii(end-1)]; 

figure(1); clf; 
for i = 1:length(subs)
    s = subs(i); 
    
    subplot(1,3,i); 
    vm_ll = allllsm(s,:,1); 
    m_ll  = allllsm(s,:,2); 
    
    plot(m_ll, vm_ll, 'x'); 
    refline(1,0); 
    
    xlabel('Mean OOSLL')
    ylabel('Var+Mean OOSLL')
    
    title(sprintf( 'Subject %d Var+Mean vs. Mean OOSLL', s)); 
    
end

filename = fullfile('..','..', 'images', 'vm_m_cv_scores.png'); 
print(filename, '-dpng'); 



