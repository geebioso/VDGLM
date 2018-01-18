function [] = sample_null_distribution( )

%% bootstrap.m

% The goal of this function is to compute the empirical distribution of beta and
% sigma parameters. The function will randomly sample (with replacement)
% each time series from each subject 1000 times, compute parameter
% estimates, and compute the variance of each parameter estimate

% Step 1: Resample all the data
% Step 2: Compute parameter estimates
% STep 3: Post-hoc computation of distributions

%% Options
whsim = 7; 
NB = 1000; % number of bootstrap samples
doparallel = 1;
isHPC = 0; 

%% Open Parallel Pool 
if doparallel
    % set number of workers to minimum of number of sims, or number of cores - 1
    nWorkers = feature('numCores');
    
%     core_info = evalc('feature(''numcores'')');
%     [a,b] = regexp( core_info, 'MATLAB detected: \d logical cores'); 
%     core_info = core_info(a:b); 
%     nWorkers= regexp(core_info,'\d','Match');
%     nWorkers = str2num( nWorkers{1}); 
    
    p = gcp('nocreate'); % If pool, do not create new one.
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    
    if ~(poolsize == nWorkers)
        delete(gcp('nocreate'))
        
        parpool( nWorkers);
    end
    fprintf('Parallel using %d cores\n', nWorkers);   
else
   nWorkers = 0;  
   fprintf('Running in serial\n'); 
end

 

%% Set Simulation
[ whs, K, seed, Tremove, doconstrained, prewhiten, ...
    var_log_transform, TukN, roifile, designfile, combdesignfile, runmodels, ...
    output_directory, multivariate, dotest] = set_analysis_options(whsim, isHPC);

fprintf('Running Simulation %d\n', whsim);
%% Set random seed
rng( seed );

%% Load the ROI timecourse data and Design if not already in memory
if ~exist('tcn', 'var')
    [tcn, XS, design, designlabels, T, NS, R ] = load_data_and_design( roifile, ...
        designfile, combdesignfile, multivariate, whs, Tremove, dotest);
end

fprintf( '%d subjects, %d regions\n', NS, R); 

%% Optimization Options
if doconstrained == 0
    fprintf('Running Unconstrained\n'); 
    max_iter = 1000;
    optim_opts = optimoptions('fminunc','Algorithm','trust-region',...
        'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
        'FiniteDifferenceType', 'central', ...
        'Diagnostics', 'off', 'MaxIterations', max_iter,...
        'HessianFcn','objective'); 
else
    fprintf('Running Constrained\n'); 
    max_iter = 1000;
    % I specify the hessian function later on in the code 
    optim_opts = optimoptions('fmincon','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
        'Diagnostics', 'off', 'MaxIterations', max_iter, ...
        'SpecifyConstraintGradient', true);
end
% Temporary
%optim_opts.Display = 'iter';


% save filename
if dotest
    filenm = [ output_directory filesep sprintf( 'paramswm_cc_whs%d_test_sampnull' ,whsim)];
    % optim_opts.Display = 'iter';
    optim_opts.Display = 'none';
else
    optim_opts.Display = 'none';
    filenm = [ output_directory filesep sprintf( 'paramswm_cc_whs%d_sampnull' ,whsim)];
end

fprintf('Save filename = %s\n', filenm); 

%% Determine the designs for each model
% Number of models
M = size( runmodels , 1 );

% Cell array with information for each model
models = cell( M,1 );
badchol = []; % for prewhitening store if any of the pre-whitening matrices can't be cholesky decomposed
for m=1:M
    
    code = runmodels{ m };
    
    % Store the model description
    model.code = code;
    model.meaneffects = code{ 1 };
    model.vareffects  = code{ 2 };
    model.description = code{ 3 };
    
    fprintf( 'Creating design for model %d of %d (%s)\n' , m , M , model.description );
    
    [ ismem , meancolsnow ] = ismember( model.meaneffects , designlabels );
    if any( ismem == 0 )
        error( 'Could not find match to some strings' );
    end
    
    [ ismem , varcolsnow ] = ismember( model.vareffects , designlabels );
    if any( ismem == 0 )
        error( 'Could not find match to some strings' );
    end
    
    model.meancols = meancolsnow;
    model.varcols  = varcolsnow;
    P1 = length( meancolsnow );
    P2 = length( varcolsnow );
    
    % Save the design for this model
    ud = unique( [ meancolsnow varcolsnow ] );
    model.designlabels = designlabels( ud );
    
    % Determine appropriate initial conditions
    initsmean = zeros( 1 , P1 );
    % Set all parameters to 0.1
    initsmean( : ) = 0.1;
    % Set intercept to 0
    wh1 = find( ismember( model.meaneffects , 'Intercept' ));
    if ~isempty( wh1 )
        initsmean( wh1 ) = 0;
    end
    
    initsvar = zeros( 1 , P2 );
    % Set all parameters to 0.1
    initsvar( : ) = 0.1;
    % But set intercept to 1
    wh2 = find( ismember( model.vareffects , 'Intercept' ));
    if ~isempty( wh2 )
        initsvar( wh2 ) = 1;
    end
    
    model.inits = [ initsmean initsvar ]';
    
    % output which parameters the current model contains
    str1 = cell( P1 , 1 );
    for i=1:P1
        str1{ i } = sprintf( 'Beta%d (%s)' , i , designlabels{ meancolsnow( i ) } );
        fprintf( '\tParameter: %s\n' , str1{ i } );
    end
    str2 = cell( P2 , 1 );
    for i=1:P2
        str2{ i } = sprintf( 'Sigma%d (%s)' , i , designlabels{ varcolsnow( i ) } );
        fprintf( '\tParameter: %s\n' , str2{ i } );
    end
    
    model.paramlabels = [ str1; str2 ];
    model.whsim = whsim;
    
    % Create a flag to indicate which parameter is for mean and which is for variance
    model.ismeanparam = zeros( 1,P1+P2 );
    model.ismeanparam( 1:P1 ) = 1;
    model.isvarparam = zeros( 1,P1+P2 );
    model.isvarparam( P1+1:P1+P2 ) = 1;
    
    % Create a flag to indicate which parameter is for intercept
    model.isinterceptparam = zeros( 1,P1+P2 );
    model.isinterceptparam( [ wh1 P1+wh1 ] ) = 1;
    
    % indicate if we have done prewhitening
    if prewhiten
        model.description = [model.description '(prewhite)'];
    end
    models{ m } = model;
end

%% Storage
% Store results in cell arrays (matlab complains when using parallization with arrays)
allparams = cell( NS , R , M, NB );
allbic    = zeros( NS , R , M, NB );


%% Fit the models
tic;
parfor (i=1:NS, nWorkers)
    sub_bic = zeros(R, M, NB);
    sub_params = cell( R, M, NB );
    
    fprintf( 'Fitting models ..... working on subject %d of %d\n', i, NS);
    YS = double( squeeze( tcn( : , i , : )));
    
    X = XS{i};
    % Loop over regions/voxels
    for j=1:R
        if rem(j, 100) == 0
            fprintf( '\tregion %d\n' , j );
        end
        
        Y = YS(:,j);
        YB = zeros(T, NB);
        for b = 1:NB
            % YB(:,b) = datasample(Y, T); 
            perm_ii = randperm(T); 
            YB(:,b) = Y(perm_ii); 
        end
        
        %% Bootstrap Resampling
        
        for b = 1:NB
            YALL = YB(:,b);
            %%  Prewhiten
            % based on Temporal Autocorrelation .... Woolrich 2001
            if prewhiten
                % compute raw autocorrelation
                if TukN < 0
                    TukN = round(2*sqrt(T));
                end
                A = autocorr_woolrich(YALL,TukN - 1 );
                A = A(1:end);
                
                % Tukey Taper
                power = (1:T)';
                Tuk = 0.5*( 1 + cos( (pi*power)/TukN));
                Tuk(1:TukN) = Tuk(1:TukN).*A;
                Tuk((TukN+1):end) = 0;
                
                % compute sample autocovariance from Tukey Taper
                autocov = Tuk*var(YALL); % see Woolrich 2001 Temporal Autocorrelation
                S = toeplitz( autocov );
                
                % compute cholesky decomposition, if unable to compute
                % decomposition, don't prewhiten ( L = eye(size(S)) )
                [L, p] = chol(S, 'lower');
                if p ~= 0
                    fprintf('S is not positive definite ( s=%d, r=%d), setting autocovariance to identity matrix\n', i, j);
                    badchol = [badchol; [i j]];
                    L = eye(size(S));
                end
                
                % do pre-whitening
                YALL_pre = L\YALL;
                XS_pre = L\XS{i};
                
            end
            
            
            for m=1:M
                %fprintf( 'i=%d j=%d f=%d m=%d\n' , i , j , f , m );
                
                % Which columns to include for the mean effect
                meancolsnow = models{ m }.meancols;
                
                % Which columns to include for the variance effect
                varcolsnow  = models{ m }.varcols;
                
                % Initial parameters
                initsnow = models{ m }.inits;
                
                
                if var_log_transform
                    % Define the (negative) likelihood to minimize
                    fun = @(x) loglik_varmean_matrix_logtransform( x, X( : , meancolsnow ), X( : , varcolsnow ) , Y);
                else
                    % Define the (negative) likelihood to minimize
                    fun = @(x) loglik_varmean_matrix_var( x, X( : , meancolsnow ), X( : , varcolsnow ) , Y);
                    
%                     % define hessian function 
%                     hessfun = @(x, lambda) hessianfcn(x, X( : , meancolsnow ), X( : , varcolsnow ), Y, lambda); 
%                     optim_opts.HessianFcn = hessfun; 
                end
                
                if ~(length(unique(Y)) == 1)
                    % Run the unconstrained optimizer
                    if (doconstrained==0)
                        [paramsnow,fval,exitflag,output,grad,hessian] = fminunc(fun, initsnow, optim_opts);
                    else
                        x0 = initsnow;
                        A = [];
                        vb = [];
                        Aeq = [];
                        beq = [];
                        lb  = [];
                        ub  = [];
                        
                        % Define the non-linear constraint function
                        nonlcon = @(x) varconstraint2( x, X( : , meancolsnow ), X( : , varcolsnow ));
                        
                        % Run the constrained optimizer
                        [paramsnow,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,vb,Aeq,beq,lb,ub,nonlcon,optim_opts);
                    end
                else
                    paramsnow = nan(size(initsnow));
                end
                
                % compute BIC
                P = length( initsnow );
                bicnow = log(T)*P + 2*fun(paramsnow);
                
                % Only store the parameters and predictions when training data = full data
                % allparams{i,j,m,b}  = paramsnow';
                sub_params{j,m,b}  = paramsnow';
                
                if size(paramsnow',1) ~= 1
                    error('objective function should produce column vector');
                end
                % allbic(i,j,m,b)     = bicnow;
                sub_bic(j,m,b) = bicnow;
            end
        end
        
    end
    allparams{i} = sub_params;
    allbic(i,:,:,:) = sub_bic;
end
runtime = toc;

fprintf('Total Run Time = %2.2f for NS=%d,R=%d,M=%d,NB=%d,%d cores\n', runtime, ...
    NS, R, M, NB, nWorkers);
% Total Run Time = 536.15 for NS=2,R=1,M=5,NB=1000,2 cores
% Total Run Time = 920.58 for NS=2,R=1,M=5,NB=1000,0 cores


%% Convert cell arrays to arrays
for m=1:M
    % Get number of parameters
    P = length( models{m}.inits );
    mat_allparams = zeros( NS, P, R, NB ); 
    for i = 1:NS
       sub_params = allparams{i}; 
       sub_params = reshape( cell2mat( sub_params(:,m,:)), P, R, NB); 
       mat_allparams(i,:,:,:) = sub_params; 
    end
    models{ m }.allparams = mat_allparams; 
    % reshape( cell2mat( allparams(:,:,m, :)) , NS , P , R, NB );
end

% Save Data
savedwhsim = whsim;
save( filenm , 'allbic', 'models' , 'savedwhsim', 'runtime', '-v7.3');

end