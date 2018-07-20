
% This function will run things in parallel on my desktop. For when
% the queue is terrible on the HPC.

jobtype = 'analyze'; % 'null' or 'analyze': 'null' corresponds to the null hypothesis testing routine, 'analyze` omits the null hypothesis testing 
whsim = 36;          % which simulation to run (see set_analysis_options_v2.m) 
isHPC = 0;           % are we running analysis remotely or on the UCI HPC? 
dotest = 0;          % are we running in test mode? NOTE: test mode only supports 5 subjects and 3 ROIs 
start_sub = 1;       % subject to start on  
end_sub = 875;       % last subject to run 
Nworkers = 8;        % number of parallel cores to use 
logging = 'INFO';    % set logging console information (see log4 for MATLAB contained in utils) 
logfile = 'test_log.txt'; % set log outfile, this is currently deprecated because i set file output to be off, but the options is required by log4
subs = start_sub:end_sub;
Nsamp = 1;           % number of samples for the null hypothesis testing routine 

%% Start parallel loop

% set number of workers to minimum of number of sims, or number of cores

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end

if ~(poolsize == Nworkers)
    delete(gcp('nocreate'))
    
    parpool(Nworkers);
    fprintf('Parallel using %d cores\n', Nworkers);
end

%% Run Jobs
parfor i = 1:length(subs)
    subnow = subs(i);
    
    switch jobtype
        case 'analyze'
            analyzedata_batch(whsim, dotest, isHPC, subnow, subnow, logging, logfile);
        case 'null'
            null_sample_hypothesis_test(whsim, dotest, isHPC, subnow, subnow, Nsamp, logging, logfile)
    end
    
end
