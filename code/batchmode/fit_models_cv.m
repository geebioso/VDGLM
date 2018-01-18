function [ paramsnow, ismotionparam, bicnow, lloutofsample, predm, predv, badchol ] = fit_models_cv(...
    models, X, motionnow, Y, var_log_transform, doconstrained, TukN, ...
    optim_opts, whtrain_sets, whtest_sets, LOG, i, j)

% This function returns the parameter estimates for a pair of models of the
% form (independent, dependent) such that the estimates of the dependent
% model (VDGLM) depend on the estimates of the independent model (GLM). The
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
K = length(whtrain_sets); % the last fold is all data

%% Find model dependencies
model_pairs = [];
for m = 1:M
    
    % get independent models and find the models that depend on them
    if ~isempty(models{m}.code{3}{2})
        
        dep_model = models{m}.code{3}{1};
        ind_model = models{m}.code{3}{2};
        
        % find independent model
        ind_idx = -1;
        for m2 = 1:M
            
            if strcmp( models{m2}.code{3}{1}, ind_model )
                ind_idx = m2;
                % store that model m must be run before model m2
                model_pair = [m ind_idx];
                
                % check that the design matrices match
                if isequal(models{m}.code{1}, models{ind_idx}.code{1})
                    LOG.debug('DEBUG', sprintf('%s <- %s, mean design matrices match', ind_model, dep_model));
                    model_pairs = [model_pairs; model_pair];
                else
                    LOG.error('ERROR', sprintf('%s <- %s: Mean design matrices do not match for dependent models', ind_model, dep_model));
                end
            end
        end
        
        % in case independent model is not found
        if ind_idx == -1
            LOG.error('ERROR', sprintf('independent model %s is missing for dependent model %s', ind_model, dep_model));
        end
    end
end

%% Fit models by Model Pairs

paramsnow = cell( M, 1 );
bicnow = cell( M, 1 );
lloutofsample = cell( K, M );
predm = cell( M, 1 );
predv = cell( M, 1 );
ismotionparam = cell( M, 1 );

for mp = 1:size( model_pairs , 1 )
    dep_idx = model_pairs( mp , 1 );
    ind_idx = model_pairs( mp , 2 );
    
    % Which columns to include for the mean effect
    meancolsnow = models{ dep_idx }.meancols;
    
    % Which columns to include for the variance effect
    varcolsnow  = models{ dep_idx }.varcols;
    
    % set up design matrices
    Xm = [X(:, meancolsnow), motionnow]; % if no motion, motionnow will be empty
    Xv = X(:, varcolsnow );
    
    % Initial parameters for VDGLM optimization
    initsnow = [models{ dep_idx }.initsmean, ... % mean experiment regressors
        models{ dep_idx}.initsmotion*ones(1, size(motionnow,2)), ... % motion regressors
        models{ dep_idx}.initsvar]'; % variance parameters
    
    % indicate which paramters are motion parameters (we won't end up
    % saving them)
    ismotionparam{ ind_idx } = [false(1, length(meancolsnow)), true(1, size(motionnow, 2)), false(1)];
    ismotionparam{ dep_idx } = [false(1, length(meancolsnow)), true(1, size(motionnow, 2)), false(1, length(varcolsnow))];
    
    %% Prewhiten
    
    % %     % Option 1: Prewhiten on each training and test fold
    % %     % Would go inside the fold for loop
    % %     [Y_pre_train, Xm_pre_train, B_pre, sigma2_pre, ~, ~ , ~, badcholtrain ] = solve_glm_prewhiten(...
    % %         Xm( trainnow, : ), Y( trainnow ), TukN );
    % %     [Y_pre_test, Xm_pre_test, ~, ~, ~, ~ , ~, badcholtest ] = solve_glm_prewhiten(...
    % %         Xm( testnow, : ), Y( testnow ), TukN );
    % %
    % %     if badcholtrain
    % %         LOG.info( 'INFO', 'bad training autocovariance matrix');
    % %     end
    % %
    % %     if badcholtest
    % %         LOG.info( 'INFO', 'bad test autocovariance matrix');
    % %     end
    % %
    % %     badchol = or(badcholtrain, badcholtest);
    
    % Option 2: Prewhiten on all the data
    if j == 2
       x = 1;  
    end
    [Y_pre, Xm_pre, B_pre, sigma2_pre, ~, ~ , ~, badchol ] = solve_glm_prewhiten(...
        Xm, Y, TukN );
    
    %% Fit Each Fold 
    for f = 1:K
        
        trainnow = whtrain_sets{f};
        testnow = whtest_sets{f};
        
        ntrain = length( trainnow );
        ntest  = length( testnow );
        
        % set train and test for variance matrix
        Xv_train = Xv( trainnow, : );
        Xv_test  = Xv( testnow, : );
        
        Xm_pre_train = Xm_pre( trainnow, : );
        Xm_pre_test  = Xm_pre( testnow,  : );
        
        Y_pre_train = Y_pre( trainnow );
        Y_pre_test  = Y_pre( testnow  );
        
        %%  Define the Likelihoods and Prediction Functions
        % for independent and dependent models
        if var_log_transform
            fun = @(x) loglik_varmean_matrix_logtransform( x,Xm_pre_train, Xv_train , Y_pre_train);
            ind_fun = @(x) loglik_varmean_matrix_logtransform( x,Xm_pre_train, Xv_train( :, 1 )  , Y_pre_train);
            
            % Define the (negative) out-of-sample likelihood
            funtest = @(x) loglik_varmean_matrix_logtransform( x ,Xm_pre_test, Xv_test  , Y_pre_test);
            ind_funtest = @(x) loglik_varmean_matrix_logtransform( x ,Xm_pre_test, Xv_test( :, 1 )  , Y_pre_test);
            
            % Define the predicted activation function
            predfun = @(x) preds_varmean_matrix_logtransform( x,Xm_pre_test, Xv_test );
            ind_predfun = @(x) preds_varmean_matrix_logtransform( x,Xm_pre_test, Xv_test( :, 1 ) );
        else
            fun = @(x) loglik_varmean_matrix_var( x,Xm_pre_train, Xv_train , Y_pre_train);
            ind_fun = @(x) loglik_varmean_matrix_var( x,Xm_pre_train, Xv_train( :, 1 )  , Y_pre_train);
            
            % Define the (negative) out-of-sample likelihood
            funtest = @(x) loglik_varmean_matrix_var( x ,Xm_pre_test, Xv_test  , Y_pre_test);
            ind_funtest = @(x) loglik_varmean_matrix_var( x ,Xm_pre_test, Xv_test( :, 1 )  , Y_pre_test);
            
            % Define the predicted activation function
            predfun = @(x) preds_varmean_matrix_var( x,Xm_pre_test, Xv_test );
            ind_predfun = @(x) preds_varmean_matrix_var( x,Xm_pre_test, Xv_test( :, 1 ) );
        end
        
        %% Fit Independent Model (GLM-OLS)
        [B_pre, sigma2_pre, ~] = solve_glm(Xm_pre_train, Y_pre_train);
        
        % store independent model parameters
        indparams =  [B_pre; sigma2_pre];
        
        %% Fit Dependent Model (VDGLM-optimization)
        
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
                nonlcon = @(x) varconstraint2( x, Xm_pre_train, Xv_train );
                
                % Run the constrained optimizer
                [depparams,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optim_opts);
            end
        else
            depparams = nan(size(initsnow));
        end
        
        %% Store parameters
        if f == K
            paramsnow{ dep_idx } = depparams';
            paramsnow{ ind_idx } = indparams';
            
            % compute adn store BIC, out-of-sample log-likelihood, predictions
            bicnow{ dep_idx } = log( ntrain )*length(depparams) + 2*fun( depparams );
            bicnow{ ind_idx } = log( ntrain )*length(indparams) + 2*ind_fun( indparams );
            
            % Compute predicted means and stds
            [ predm{ dep_idx } , predv{ dep_idx } ] = predfun( depparams );
            [ predm{ ind_idx } , predv{ ind_idx } ] = ind_predfun( indparams );
        else
            lloutofsample{ f, dep_idx } = -funtest(depparams);
            lloutofsample{ f, ind_idx } = -ind_funtest(indparams);
        end
        
    end
end



