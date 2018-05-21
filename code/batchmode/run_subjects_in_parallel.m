
% This function will run things in parallel on my desktop. For when
% the queue is terrible on the HPC.

jobtype = 'analyze'; % 'null' or 'analyze'
whsim = 36;
isHPC = 0;
dotest = 0;
start_sub = 1;
end_sub = 875;
Nworkers = 8;
logging = 'INFO';
logfile = 'test_log.txt';
subs = start_sub:end_sub;
Nsamp = 1; 

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
            analyzedata_batch_v2(whsim, dotest, isHPC, subnow, subnow, logging, logfile);
        case 'null'
            null_sample_hypothesis_test(whsim, dotest, isHPC, subnow, subnow, Nsamp, logging, logfile)
    end
    
end
