%% varGLMconvergence.m
function [rhat] = computeRhat( filename, df, verbose )

% df is degrees of freedom of t distribution

if df < 3
    error('degrees of freedom are too small');
end

% filename = fullfile('Simulation', 'inference.mat');
load(filename);
%samples = src.extract('permuted', false);
nchains = size(samples,2);
niter = size(samples(1).lp__,1);

fieldnms = fieldnames(samples);
rhat = struct();
for i = 1:length(fieldnms)
    
    fieldnow = fieldnms{i};
    dims = size(samples(1).(fieldnow));
    
    allsamps = [];
    withinvar = zeros([nchains, dims(2:end)]);
    for c = 1:nchains
        sampnow = samples(c).(fieldnow);
        allsamps = [allsamps; sampnow];
        withinvar(c,:,:,:,:) = var(sampnow,1);
    end
    betweenvar = squeeze(var(allsamps,1));
    
    var_ratio = betweenvar./(squeeze(mean( withinvar,1)));
    df_ratio = df/(df-2);
    rhatnow = (niter - 1)/niter + (nchains+1)/(nchains*niter)*var_ratio*df_ratio;
    
    rhat.(fieldnow) = rhatnow;
    if verbose
        if any( abs(1 - rhatnow) > 0.01)
            fprintf('%20.20s DID NOT CONVERGE\n', fieldnow);
        else
            fprintf('%20.20s CONVERGED\n', fieldnow);
        end
    end
    
    
end
end