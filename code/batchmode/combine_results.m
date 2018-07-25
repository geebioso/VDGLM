function [] = combine_results( whsim, isHPC , dotest, logging, set_up_directory_structure, jobtype )

% This is just a wrapper function for combine_results_analyze_data.m and
% combine_results_null_sample_hypothesis_test.m
% LAZY, LAZY Coding, but it works 

% INPUT: 
%   numeric whsim: which simulation 
%   bool isHPC: are we doing this on the HPC? 
%   bool dotest: are we in test mode? 
%   str logging: which logging mode
%   set_up_directory_structure: autopopulate directory structures? 
%   jobtype: 'null' or 'analyze'

%% Combine Results 

if strcmp( jobtype , 'analyze' )
    combine_results_analyze_data( whsim, isHPC , dotest, logging, set_up_directory_structure )
elseif strcmp( jobtype, 'null')
    combine_results_null_sample_hypothesis_test( whsim, isHPC , dotest, logging, set_up_directory_structure )
end