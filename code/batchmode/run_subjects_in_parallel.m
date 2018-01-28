
% This function will run things in parallel on my desktop. For when
% the queue is terrible on the HPC.


whsim = 26;
isHPC = 0;
dotest = 0;
start_sub = 9;
end_sub = 875;
Nworkers = 8;
logging = 'INFO';
logfile = 'test_log.txt';
subs = start_sub:end_sub;

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
    analyzedata_batch_v2(whsim, dotest, isHPC, subnow, subnow, logging, logfile);
    
end
