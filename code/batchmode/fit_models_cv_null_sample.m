function [ paramsnow, ismotionparam, bicnow, lloutofsample, badchol ] = fit_models_cv_null_sample(...
    models, Xm, ismotionparam, Xv, Y, var_log_transform, doconstrained, TukN, ...
    prewhiten, optim_opts, whtrain_sets, whtest_sets, LOG, i, j, Nsamp)

% This function is similar to the function fit_models_cv.m, but only has to
% fit the VDGLM and so has some tweaks to make it faster.

% This function returns the parameter estimates for a pair of models of the
% form (GLM, VDGLM) such that the estimates of the VDGLM
% model depend on the estimates of the GLM model. The
% reason for setting the code up this way is that the prewhitening procedure
% depends on the solution of the GLM. Thus the entire fitting procedure in
% this function is as follows:
%   - Fit GLM0
%       - prewhiten
%           - Fit GLM1
%           - Fit VDGLM0
%   We fit the GLM models using OLS (which is faster than our optimization
%   procedure).

%   Model dependencies are hard-coded (see set_analysis_options_v2.m) in the
%   following way:
%       Int  <- Var
%       Mean <- Var+Mean
%   Where <- indicates that the left model must run for the right model to run.

% INPUT:
%   cellarray models: storage for results
%   numeric X: full design matrix (all mean and var columns), may have some motions
%       regressor columns
%   numeric motionnow: the motion for this subject/region, will be added to design
%       matrix
%   numeric Y: BOLD response
%   bool var_log_transform: are we fitting the transformed version of the VDGLM?
%   bool doconstrained: are we running constrained optimization?
%   numeric or NaN TukN: size of the Tukey Taper
%   struct optim_opts: optimization options
%   cell whtrain_sets: training fold indices
%   cell whtest_sets: test fold indices
%   log4m LOG: logger
%   int i: subject number
%   int j: region number
%   int Nsamp: number of samples to take

% OUTPUT:
%   cell paramsnow: paramter estimate for all models fit. index into this
%       array matches the index into the variable 'models'
%   cell bicnow: bic for all models fit
%   cell lloutofsample: out of sample log likelihood
%   cell predm: mean predictions
%   cell predv: variance predictions
%   bool badchol: was there a problem with solving the root of the
%       autocovariance matrix during prewhitening?

LOG.debug('DEBUG', sprintf('fitting subject %d, region %d', i, j));
T = size(Y,1);
M = length(models);
K = length(whtrain_sets) - 1; % the last fold is all data

if ~isdeployed
    addpath(fullfile('..', 'optimization')); 
end

%% Estimate the Whitening Matrix 
if prewhiten
    [ ~, ~, ~, ~, W, badchol, df ] = solve_glm( Xm, Y, prewhiten, TukN);
    
    % get inverse of whitening matrix for generating autocorreled noise 
    Winv = inv(W);
    if rcond(Winv) < eps
        badchol = 1;
    else
        badchol = 0;
    end
    
end

% solve the non-prewhitened GLM
[~, ~, B, sigma2, ~, ~, df] = solve_glm( Xm, Y, 0, TukN);
Yhat = Xm*B;

if badchol
    LOG.warn('WARN', sprintf('S is not positive definite subject %d, region %d', i, j));
end

%% Run Sampling Routine 
paramsnow = cell( M , Nsamp );
bicnow = cell( M , Nsamp );
lloutofsample = cell( K+1 , M , Nsamp );

vdglm_idx = 1;
glm_idx = 2;

for n = 1:Nsamp
    
    %% Generate Sample
    % add noise to the non-prewhitened mean trend
    % use the biased estimation of the variance because the autocorrelation
    % estimate is also biased (biased led to better AR parameter recovery)
    noise = mvnrnd(zeros(T, 1), sigma2*(df/T)*eye(T))';
    Ysamp = Yhat + Winv*noise; % removing noise would be Ypre = W * Ysamp or Ypre = Winv \ Ysamp;
    
    %% Model Fitting 
    
    % prewhiten generated time series 
    if prewhiten
        [Y_pre, Xm_pre, ~, ~, ~, badchol, df ] = solve_glm( Xm, Ysamp, prewhiten, TukN);
    end
    
    %% Fit Each Fold
    for f = 1:K+1
        
        trainnow = whtrain_sets{f};
        testnow = whtest_sets{f};
        
        ntrain = sum( trainnow );
        
        % set train and test for variance matrix
        Xv_train = Xv( trainnow, : );
        Xv_test  = Xv( testnow,  : );
        
        if prewhiten
            Xm_train = Xm_pre( trainnow, : );
            Xm_test  = Xm_pre( testnow,  : );
            
            Y_train = Y_pre( trainnow );
            Y_test  = Y_pre( testnow  );
        else
            Xm_train = Xm( trainnow, : );
            Xm_test  = Xm( testnow,  : );
            
            Y_train = Y( trainnow );
            Y_test  = Y( testnow  );
        end
        
        % Define the Likelihoods and Prediction Functions
        if var_log_transform
            ind_fun = @(x) loglik_varmean_matrix_logtransform( x,Xm_train, Xv_train( :, 1 )  , Y_train);
            
            % Define the (negative) out-of-sample likelihood
            ind_funtest = @(x) loglik_varmean_matrix_logtransform( x ,Xm_test, Xv_test( :, 1 )  , Y_test);
            
        else
            ind_fun = @(x) loglik_varmean_matrix_var( x,Xm_train, Xv_train( :, 1 )  , Y_train);
            
            % Define the (negative) out-of-sample likelihood
            ind_funtest = @(x) loglik_varmean_matrix_var( x ,Xm_test, Xv_test( :, 1 )  , Y_test);
            
        end
        
        %% Fit Independent Model (GLM-OLS)
        [~, ~, B, sigma2, ~, ~] = solve_glm( Xm_train, Y_train, 0 ); % do not prewhiten (we will have already done it above)
        
        % store GLM model parameters
        glmparams =  [B; sigma2];
        
        %% Fit Dependent Model (VDGLM-optimization)
        if ~isnan( vdglm_idx )
            
            % Which columns to include for the variance effect
            varcolsnow  = models{ vdglm_idx }.varcols;
            
            %             % Initial parameters for VDGLM optimization
            %             initsnow = [models{ vdglm_idx }.initsmean, ... % mean experiment regressors
            %                 models{ vdglm_idx}.initsmotion*ones(1, sum( ismotionparam{ vdglm_idx} )), ... % motion regressors
            %                 models{ vdglm_idx}.initsvar]'; % variance parameters
            
            % use GLM as an estimate of the mean parameters. Should
            % converge faster.
            initsnow = [ B',  1, models{vdglm_idx}.initsvar(2:end)]';
            
            % Define the Likelihoods and Prediction Functions
            if var_log_transform
                fun = @(x) loglik_varmean_matrix_logtransform( x,Xm_train, Xv_train , Y_train);
                
                % Define the (negative) out-of-sample likelihood
                funtest = @(x) loglik_varmean_matrix_logtransform( x ,Xm_test, Xv_test  , Y_test);
                
                % Define the predicted activation function
                predfun = @(x) preds_varmean_matrix_logtransform( x,Xm_test, Xv_test );
            else
                fun = @(x) loglik_varmean_matrix_var( x,Xm_train, Xv_train , Y_train);
                
                % Define the (negative) out-of-sample likelihood
                funtest = @(x) loglik_varmean_matrix_var( x ,Xm_test, Xv_test , Y_test);
                
                % Define the predicted activation function
                predfun = @(x) preds_varmean_matrix_var( x,Xm_test, Xv_test );
            end
            
            if ~(length(unique(Y)) == 1) % i.e., if the data is from within the brain
                if (doconstrained==0) % run the unconstrained optimizer
                    [vdglmparams,fval,exitflag,output,grad,hessian] = fminunc(fun, initsnow, optim_opts);
                else % run the constrained optimizer
                    x0 = initsnow;
                    A = [];
                    b = [];
                    Aeq = [];
                    beq = [];
                    lb  = [];
                    ub  = [];
                    
                    % Define the non-linear constraint function
                    nonlcon = @(x) varconstraint2( x, Xm_train, Xv_train );
                    
                    % Run the constrained optimizer
                    [vdglmparams,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optim_opts);
                end
            else
                vdglmparams = nan(size(initsnow));
            end
            
        end
        %% Store parameters
        if f == K+1
            
            % compute and store BIC, out-of-sample log-likelihood, predictions
            if ~isnan( vdglm_idx )
                paramsnow{ vdglm_idx, n } = vdglmparams';
                bicnow{ vdglm_idx, n } = log( ntrain )*length(vdglmparams) + 2*fun( vdglmparams );
            end
            
            paramsnow{ glm_idx, n } = glmparams';
            bicnow{ glm_idx, n } = log( ntrain )*length(glmparams) + 2*ind_fun( glmparams );
        else
            
            if ~isnan( vdglm_idx )
                lloutofsample{ f, vdglm_idx, n } = -funtest(vdglmparams);
            end
            lloutofsample{ f, glm_idx, n } = -ind_funtest(glmparams);
        end
        
    end
end
end



