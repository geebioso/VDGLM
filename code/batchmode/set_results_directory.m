function [results_directory] = set_results_directory( isHPC, set_up_directory_structure )

% INPUT: 
%   isHPC: are we running analysis on the HPC (0/1) 
%   set_up_directory_structure: do you want to automatically populate the
%       results directory with the appropriate subdirectories? (0/1) 
%       MATLAB will issue a warning if the directories already exist 

if isHPC
    results_directory = '/pub/ggaut/VDGLM/Results';
else
    results_directory = fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'Projects',...
        'varianceGLM', 'Results');
    
end

if set_up_directory_structure
    
    wd = pwd; 
    
    cd(results_directory) 
   
    mkdir single_analyses
    mkdir batch_analyses
    mkdir batch_analyses/combined
    mkdir batch_analyses/single_jobs
    mkdir batch_analyses/null_single_jobs
    mkdir batch_analyses/null_combined 
    
    mkdir batch_analyses/combined_test
    mkdir batch_analyses/single_jobs_test
    mkdir batch_analyses/null_single_jobs_test
    mkdir batch_analyses/null_combined_test
    
    cd(wd); 
end 

