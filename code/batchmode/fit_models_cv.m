function [ paramsnow, ismotionparam, bicnow, lloutofsample, predm, predv, badchol ] = fit_models_cv(...
    models, X, motionnow, scrubnow, Y, var_log_transform, doconstrained, TukN, ...
    prewhiten, optim_opts, whtrain_sets, whtest_sets, LOG, i, j)

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

%% Find model dependencies
model_pairs = [];
for m = 1:M
    
    % get GLM models and find the models that depend on them
    if ~isempty(models{m}.code{3}{2})
        
        vdglm_model = models{m}.code{3}{1};
        glm_model = models{m}.code{3}{2};
        
        % find GLM model
        glm_idx = -1;
        for m2 = 1:M
            
            if strcmp( models{m2}.code{3}{1}, glm_model )
                glm_idx = m2;
                % store that model m must be run before model m2
                model_pair = [m glm_idx];
                
                % check that the design matrices match
                if isequal(models{m}.code{1}, models{glm_idx}.code{1})
                    LOG.debug('DEBUG', sprintf('%s <- %s, mean design matrices match', glm_model, vdglm_model));
                    model_pairs = [model_pairs; model_pair];
                else
                    LOG.error('ERROR', sprintf('%s <- %s: Mean design matrices do not match for GLM and VDGLM models', glm_model, vdglm_model));
                end
            end
        end
        
        % in case GLM model is not found
        if glm_idx == -1
            LOG.error('ERROR', sprintf('GLM model %s is missing for VDGLM model %s', glm_model, vdglm_model));
        end
    elseif and( ~isempty( models{m}.code{3}{1}), ~ismember( m, model_pairs)) % If no VDGLM models, just get GLM models
        glm_idx = m;
        model_pair = [NaN glm_idx];
        model_pairs = [model_pairs; model_pair];
    end
end

%% Scrubbing
for f = 1:K+1
    whtrain_sets{f}(scrubnow) = [];
    whtest_sets{f}(scrubnow)  = [];
end

Y(scrubnow) = [];
X(scrubnow,:) = [];
motionnow(scrubnow, :) = [];

%% Fit models by Model Pairs
paramsnow = cell( M, 1 );
bicnow = cell( M, 1 );
lloutofsample = cell( K+1, M );
predm = cell( M, 1 );
predv = cell( M, 1 );
ismotionparam = cell( M, 1 );

for mp = 1:size( model_pairs , 1 )
    vdglm_idx = model_pairs( mp , 1 );
    glm_idx = model_pairs( mp , 2 );
    
    % Which columns to include for the mean effect (GLM and
    % VDGLM models share mean design matrix) 
    meancolsnow = models{ glm_idx }.meancols;
    
    % Which columns to include for the variance effect
    if ~isnan( vdglm_idx ) 
        varcolsnow  = models{ vdglm_idx }.varcols;
    else
        varcolsnow  = models{ glm_idx }.varcols; 
    end
    
    % set up design matrices
    Xm = [X(:, meancolsnow), motionnow]; % if no motion, motionnow will be empty
    Xv = X(:, varcolsnow );
    
    % indicate which parameters are motion parameters
    ismotionparam{ glm_idx } = [false(1, length(meancolsnow)), true(1, size(motionnow, 2)), false(1)];
    
    %% Prewhiten
    % Prewhiten on all the data
    if prewhiten
        [Y_pre, Xm_pre, B_pre, sigma2_pre, L, badchol ] = solve_glm( Xm, Y, prewhiten, TukN);
    end
    
    %% Fit Each Fold
    for f = 1:K+1
        
        trainnow = whtrain_sets{f};
        testnow = whtest_sets{f};
        
        ntrain = sum( trainnow );
        ntest  = sum( testnow );
        
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
            
            % Define the predicted activation function
            ind_predfun = @(x) preds_varmean_matrix_logtransform( x,Xm_test, Xv_test( :, 1 ) );
        else
            ind_fun = @(x) loglik_varmean_matrix_var( x,Xm_train, Xv_train( :, 1 )  , Y_train);
            
            % Define the (negative) out-of-sample likelihood
            ind_funtest = @(x) loglik_varmean_matrix_var( x ,Xm_test, Xv_test( :, 1 )  , Y_test);
            
            % Define the predicted activation function
            ind_predfun = @(x) preds_varmean_matrix_var( x,Xm_test, Xv_test( :, 1 ) );
        end
        
        %% Fit Independent Model (GLM-OLS)
        [~, ~, B, sigma2, ~, badchol] = solve_glm( Xm_train, Y_train, 0 ); % do not prewhiten (we will have already done it above)
        
        % store GLM model parameters
        indparams =  [B; sigma2];
        
        %% Fit Dependent Model (VDGLM-optimization)
        if ~isnan( vdglm_idx )
            
            % Which columns to include for the variance effect
            varcolsnow  = models{ vdglm_idx }.varcols;
            
            % Initial parameters for VDGLM optimization
            initsnow = [models{ vdglm_idx }.initsmean, ... % mean experiment regressors
                models{ vdglm_idx}.initsmotion*ones(1, size(motionnow,2)), ... % motion regressors
                models{ vdglm_idx}.initsvar]'; % variance parameters
            
            % Indicate motion parameters
            ismotionparam{ vdglm_idx } = [false(1, length(meancolsnow)), true(1, size(motionnow, 2)), false(1, length(varcolsnow))];
            
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
                    [depparams,fval,exitflag,output,grad,hessian] = fminunc(fun, initsnow, optim_opts);
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
                    [depparams,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optim_opts);
                end
            else
                depparams = nan(size(initsnow));
            end
            
        end
        %% Store parameters
        if f == K+1

            % compute and store BIC, out-of-sample log-likelihood, predictions
            if ~isnan( vdglm_idx ) 
                paramsnow{ vdglm_idx } = depparams';
                bicnow{ vdglm_idx } = log( ntrain )*length(depparams) + 2*fun( depparams );
                [ predm{ vdglm_idx } , predv{ vdglm_idx } ] = predfun( depparams );
            end
            
            paramsnow{ glm_idx } = indparams';
            bicnow{ glm_idx } = log( ntrain )*length(indparams) + 2*ind_fun( indparams );
            [ predm{ glm_idx } , predv{ glm_idx } ] = ind_predfun( indparams );
        else
            
            if ~isnan( vdglm_idx ) 
                lloutofsample{ f, vdglm_idx } = -funtest(depparams);
            end
            lloutofsample{ f, glm_idx } = -ind_funtest(indparams);
        end
        
    end
end



