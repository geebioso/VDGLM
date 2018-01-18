
%% Options

whsim = 13;
dohighpassfilter = 0;
prewhiten_mean_models = 0;
prewhiten_var_models = 0;
doarima = 0; % my Econometrics Toolbox trial has expired
filenm = fullfile('Results', sprintf( 'paramswm_cc_whs%d.mat', whsim));
normal = 0;
exponential = 1;
isHPC = 0; 

%% Load a previously prepared parameter set
if ~exist( 'savedwhsim' , 'var' ) || savedwhsim ~= whsim
    fprintf( 'Loading previously saved parameters from: %s\n' , filenm );
    load( filenm );
end

% loads: 'allbicm', 'allllsm'  , 'bestmodelBIC' , 'bestmodelCV' ,'M', 'K' , 'models' , 'savedwhsim'

%% Set Simulation
[ whs, K, seed, Tremove, doconstrained, prewhiten, ...
    var_log_transform, TukN, roifile, designfile, combdesignfile, runmodels,...
    output_directory, multivariate, dotest] = set_analysis_options(whsim, isHPC);

fprintf('Running Simulation %d\n', whsim);
%% Set random seed
rng( seed );

%% Load the ROI timecourse data and Design if not already in memory
[tcn, XS, design, designlabels, T, NS, R ] = load_data_and_design( roifile, ...
    designfile, combdesignfile, multivariate, whs, Tremove, dotest);


%% Simulate Data From Our Model, Prewhiten, and Recover Parameters
if normal
    auto_estimator = 'Tuk'; % TukPAVA, PAVA, Tuk, auto
    nsims = 100;
    
    doplots = 0;
    
    beta =  [ 0; 1 ];
    sigma = [ 1; 1 ] ;
    LB = length(beta);
    initsnow = [beta;sigma];
    X = XS{1}(:, ismember(designlabels', {'2-back', 'Intercept'})); % the 2 back indicator
    % X(:,2) = X(:,2) > 0;
    in2back = X > 0;
    
    % modeltypes = {'glmauto', 'varglm', 'glm'}; % varglm or glmauto
    
    modeltypes = {'varglm'};
    
    if doconstrained == 0
        max_iter = 1000;
        optim_opts = optimoptions('fminunc','Algorithm','trust-region',...
            'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
            'FiniteDifferenceType', 'central', ...
            'Diagnostics', 'off', 'MaxIterations', max_iter); %, 'Display', 'iter-detailed');
    else
        max_iter = 1000;
        % we set the hessian function later 
        optim_opts = optimoptions('fmincon','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
        'Diagnostics', 'off', 'MaxIterations', max_iter,  ...
        'SpecifyConstraintGradient', true);
    end
    
    
    vparams_all = zeros( nsims, length( [beta; sigma]));
    mparams_all = zeros( nsims, length( [beta; 1]));
    predsm     = zeros(nsims, T );
    predsv     = zeros(nsims, T );
    YALL = zeros(nsims, T);
    for m = 1:length(modeltypes)
        modeltype = modeltypes{m};
        for s= 1:nsims
            
            Y= normrnd(X*beta, sqrt(X*sigma));
            YALL(s,:) = Y;
            
            % high pass filter
            addpath( fullfile('/Users/Garren/spm12/external/fieldtrip/external/signal'));
            
            % var model
            fun = @(x) loglik_varmean_matrix_var( x, X, X , Y);
            predfun = @(x) preds_varmean_matrix_var( x, X , X);
            hessfun = @(x, lambda) hessianfcn(x, X, X, Y, lambda); 
            % optim_opts.HessianFcn = hessfun; 
            
            if ~doconstrained
                [paramsnow,fval,exitflag,output,grad,hessian] = fminunc(fun, initsnow, optim_opts);
            else
                x0 = initsnow;
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb  = [];
                ub  = [];
                
                % Define the non-linear constraint function
                nonlcon = @(x) varconstraint2( x, X, X);
                
                % Run the constrained optimizer
                [paramsnow,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optim_opts);
            end
            
            vparams_all( s, :) = paramsnow;
            predsv(s,:) = predfun(paramsnow);
            hessfun = @(x, lambda) hessianfcn(x, X, ones(T,1), Y, lambda); 
            % optim_opts.HessianFcn = hessfun; 
            
            % mean model
            fun = @(x) loglik_varmean_matrix_var( x, X, ones(T,1), Y);
            predfun = @(x) preds_varmean_matrix_var( x, X, ones(T,1));
            
            
            if ~doconstrained 
                [paramsnow,fval,exitflag,output,grad,hessian] = fminunc(fun, initsnow(1:(LB + 1)), optim_opts);
            
            else
                x0 = initsnow(1:(LB + 1));
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb  = [];
                ub  = [];
                
                % Define the non-linear constraint function
                nonlcon = @(x) varconstraint2( x, X, ones(T,1));
                
                % Run the constrained optimizer
                [paramsnow,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optim_opts);
            end
            
            mparams_all( s, :) = paramsnow;
            predsm(s,:) = predfun(paramsnow);
        end
    end
    
    %% Calculate Residuals
    resm = zeros(nsims, T);
    resv = zeros(nsims, T);
    for i = 1:nsims
        resm(i,:) = YALL(i,:) - predsm(i,:);
        resv(i,:) = YALL(i,:) - predsv(i,:);
    end
    
    rss_v = sum( resv.^2, 2);
    rss_m = sum( resm.^2, 2);
    
    figure(3); clf
    plot(rss_v, rss_m, 'o');
    refline(1,0);
    xlabel('V\_M');
    ylabel('M');
    
    %% Plot Parameter Recovery
    
    figure(1); clf;
    
    subplot(2,2,1);
    histogram(vparams_all(:,1));
    title( sprintf('True Beta0 = %2.2f', beta(1)));
    
    subplot(2,2,2);
    histogram(vparams_all(:,2));
    title( sprintf('True Beta1 = %2.2f', beta(2)));
    
    subplot(2,2,3);
    histogram(vparams_all(:,3));
    title( sprintf('True Sigma0 = %2.2f', sigma(1)));
    
    subplot(2,2,4);
    histogram(vparams_all(:,4));
    title( sprintf('True Sigma1 = %2.2f', sigma(2)));
    
    filename = fullfile( 'images', 'exponential_recovery_predictions');
    print(filename, '-depsc');
end
%% Simulate Data From Exponential Model

if exponential
    
    nsims = 1000;
    
    doplots = 0;
    
    beta =  [ 0; -1 ];
    sigma = [ -1; 1] ;
    initsnow = [beta;sigma];
    X = XS{1}(:, ismember(designlabels', {'2-back', 'Intercept'})); % the 2 back indicator
    T = size(X, 1);
    % X(:,2) = X(:,2) > 0;
    in2back = X > 0;
    
    if doconstrained == 0
        max_iter = 1000;
        optim_opts = optimoptions('fminunc','Algorithm','trust-region',...
            'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
            'FiniteDifferenceType', 'central', ...
            'Diagnostics', 'off', 'MaxIterations', max_iter); %, 'Display', 'iter-detailed');
    else
        max_iter = 1000;
        optim_opts = optimoptions( 'fmincon','Algorithm','interior-point',...
            'SpecifyObjectiveGradient',true, 'MaxIterations', max_iter, 'Display', 'iter-detailed');
    end
    
    
    params_all = zeros( nsims, length( [beta; sigma]));
    YALL = zeros(nsims, size(X,1));
    predmall = zeros(nsims, T);
    predvall = zeros(nsims, T);
    for s= 1:nsims
        
        Y= normrnd(X*beta, sqrt(exp(X*sigma)));
        YALL(s, :) = Y;
        
        % high pass filter
        addpath( fullfile('/Users/Garren/spm12/external/fieldtrip/external/signal'));
        
        fun = @(x) loglik_varmean_matrix_logtransform_var( x, X, X , Y );
        predfun = @(x) preds_varmean_matrix_logtransform( x, X, X );
        
        if ~doconstrained
            [paramsnow,fval,exitflag,output,grad,hessian] = fminunc(fun, initsnow, optim_opts);
        else
            x0 = initsnow;
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb  = [];
            ub  = [];
            
            % Define the non-linear constraint function
            nonlcon = @(x) varconstraint( x, X, X);
            
            % Run the constrained optimizer
            [paramsnow,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optim_opts);
        end
        
        params_all( s, : ) = paramsnow;
        [ predm , predv ] = predfun( paramsnow );
        
        predmall( s, : ) = predm;
        predvall( s, : ) = predv;
        
    end
    
    %% Plot Parameter Recovery
    
    figure(2); clf;
    
    subplot(2,2,1);
    histogram(params_all(:,1));
    title( sprintf('True Beta0 = %2.2f', beta(1)));
    
    subplot(2,2,2);
    histogram(params_all(:,2));
    title( sprintf('True Beta1 = %2.2f', beta(2)));
    
    subplot(2,2,3);
    histogram(params_all(:,3));
    title( sprintf('True Sigma0 = %2.2f', sigma(1)));
    
    subplot(2,2,4);
    histogram(params_all(:,4));
    title( sprintf('True Sigma1 = %2.2f', sigma(2)));
    
    set(gcf, 'Name', 'Exponential Parameter Recovery');
    filename = fullfile( 'images', 'exponential_recovery_parameters');
    print(filename, '-depsc');
    
    %% Plot Data Versus Predictions
    lw = 2;
    h = figure( 300 ); clf;
    h.Position = [77 355 634 349];
    set( gcf , 'Name' , 'Model Predictions Overall + Data (BIC model Selection)' );
    
    % Predictions from model for these cases
    ypreds_mean = nanmedian( predmall );
    ypreds_var  = nanmedian( predvall );
    
    % get data statistics
    meany = mean( YALL);
    vary = var( YALL);
    
    % Plot the means
    subplot( 1 , 2 , 1);
    plot( 1:T , X* 0.2 + 1, '-' , 'LineWidth' , lw ); hold on;
    plot( 1:T , ypreds_mean  , 'b-' , 'LineWidth' , lw ); hold on;
    plot( 1:T , meany , 'r-' , 'LineWidth' , lw );
    % title( titlestr );
    % ylim( [ -0.5 1.5 ] );
    ylabel( 'Mean' );
    xlim( [ 0 T ] );
    
    
    subplot( 1,2,2);
    plot( 1:T ,X * 0.2 + max( vary) + 0.5, '-' , 'LineWidth' , lw ); hold on;
    plot( 1:T , ypreds_var, 'b-' , 'LineWidth' , lw ); hold on;
    plot( 1:T , vary , 'r-' , 'LineWidth' , lw );
    % title( titlestr );
    ylabel( 'Std Dev' );
    xlim( [ 0 T ] );
    
    set(gcf, 'Name', 'Exponential Predictions Recovery');
    filename = fullfile( 'images', 'exponential_recovery_predictions');
    print(filename, '-depsc');
end