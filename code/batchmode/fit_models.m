function [ paramsnow, ismotionparam, bicnow, lloutofsample, predm, predv, badchol ] = fit_models(...
    models, X, motionnow, Y, var_log_transform, doconstrained, TukN, ...
    optim_opts, trainnow, testnow, LOG)

% This function returns the parameter estimates for a pair of models of the 
% form (independent, dependent) such that the estimates of the dependent
% model (VDGLM) depend on the estimates of the independent model (GLM). The 
% reason for setting the code up this way is that the prewhitening procedure
% depends on the solution of the GLM. Thus the entire fitting procedure in
% this function is as follows: 
%   - Fit GLM0
%       - prewhiten
%       - Fit GLM1
%       - Fit VDGLM0
%   We fit the GLM models using OLS (which is faster than our optimization 
%   procedure). 

%   Model dependencies are hard-coded (see set_analysis_options_v2.m) in the 
%   following way: 
%       Int  <- Var
%       Mean <- Var+Mean 
%   Where <- indicates that the first model must run for the second to run.

% INPUT: 
%   models: storage for results 
%   X: full design matrix (all mean and var columns), may have some motions
%       regressor columns 
%   motionnow: the motion for this subject/region, will be added to design
%       matrix
%   Y: BOLD response 
%   var_log_transform (0/1): are we fitting the transformed version of the VDGLM?
%   doconstrained (0/1): are we running constrained optimization? 
%   TukN: size of the Tukey Taper 
%   optim_opts: optimization options 
%   trainnow: training index 
%   testnow: test index 
%   LOG: logger 

% OUTPUT: 
%   paramsnow: paramter estimate for all models fit. index into this
%       array matches the index into the variable 'models' 
%   bicnow: bic for all models fit 
%   lloutofsample: out of sample log likelihood
%   predm: mean predictions 
%   predv: variance predictions
%   badchol (0/1): was there a problem with solving the root of the
%       autocovariance matrix during prewhitening? 

T = size(Y,1); 
M = length(models);
ntrain = length(trainnow); 
ntest  = length(testnow); 

paramsnow = cell(M,1); 
bicnow = cell(M, 1); 
lloutofsample = cell(M, 1); 
predm = cell(M,1); 
predv = cell(M,1); 
ismotionparam = cell(M,1); 

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
    
    % Initial parameters
    initsnow = [models{ dep_idx }.initsmean, ... % mean experiment regressors
        models{ dep_idx}.initsmotion*ones(1, size(motionnow,2)), ... % motion regressors 
        models{ dep_idx}.initsvar]'; % variance parameters 
    
    % indicate which paramters are motion parameters (we won't end up
    % saving them) 
    ismotionparam{ ind_idx } = [false(1, length(meancolsnow)), true(1, size(motionnow, 2)), false(1)]; 
    ismotionparam{ dep_idx } = [false(1, length(meancolsnow)), true(1, size(motionnow, 2)), false(1, length(varcolsnow))]; 
    
    %% Solve Prewhitened GLM 
    
    % Option 1: Prewhiten on each training and test fold 
    [Y_pre_train, Xm_pre_train, B_pre, sigma2_pre, ~, ~ , ~, badcholtrain ] = solve_glm_prewhiten(...
        Xm( trainnow, : ), Y( trainnow ), TukN ); 
    [Y_pre_test, Xm_pre_test, ~, ~, ~, ~ , ~, badcholtest ] = solve_glm_prewhiten(...
        Xm( testnow, : ), Y( testnow ), TukN ); 
    
    % Option 2: Prewhiten on all the data 
    
    
    indparams =  [B_pre; sigma2_pre]; 
    
    if badcholtrain
       LOG.info( 'INFO', 'bad training autocovariance matrix');  
    end
    
    if badcholtest
       LOG.info( 'INFO', 'bad test autocovariance matrix');  
    end
    
    badchol = or(badcholtrain, badcholtest); 
    
    % set train and test for variance matrix 
    Xv_train = Xv( trainnow, : ); 
    Xv_test  = Xv( testnow, : ); 
    
    %% VDGLM optimization 
    % Define the (negative) likelihood to minimize
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
    
    if ~(length(unique(Y)) == 1)
        % Run the unconstrained optimizer
        if (doconstrained==0)
            [depparams,fval,exitflag,output,grad,hessian] = fminunc(fun, initsnow, optim_opts);
        else
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
    paramsnow{ dep_idx } = depparams'; 
    paramsnow{ ind_idx } = indparams';
    
    % compute adn store BIC, out-of-sample log-likelihood, predictions 
    bicnow{dep_idx} = log( ntrain )*length(depparams) + 2*fun( depparams ); 
    bicnow{ind_idx} = log( ntrain )*length(indparams) + 2*ind_fun( indparams ); 
    
    lloutofsample{dep_idx} = -funtest(depparams); 
    lloutofsample{ind_idx} = -ind_funtest(indparams); 
    
    % Compute predicted means and stds
    [ predm{dep_idx} , predv{dep_idx} ] = predfun( depparams );
    [ predm{ind_idx} , predv{ind_idx} ] = ind_predfun( indparams );

end



