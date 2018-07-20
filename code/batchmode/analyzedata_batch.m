function [] = analyzedata_batch(whsim, dotest, isHPC, start_sub, end_sub, logging, logfile, set_up_directory_structure)

% This function will fit a set of models (VDGLM and GLM pairs) specified
% by whsim and the function set_analysis_options_v2.m. We fit models for
% subjects start_sub:1:end_sub. The function can be run on my personal
% computer or the UCI HPC. The test flag causes the model to run with a
% reduced number of subjects and regions.

% INPUT:
%   numeric whsim: which simulation
%   bool dotest(0/1): run in test mode?
%   bool isHPC (0/1): run on HPC or personal computer
%   numeric start_sub: start subject
%   numeric end_sub: end subject
%   char logging: type of logging e.g, {'ALL','TRACE','DEBUG','INFO','WARN','ERROR','FATAL','OFF'}
%   char logfile: log file, useless argument because I have set the loglevel to off, but log4.m requires it

% NOTE: running this on the HPC will change the input type to char, hence
%   the following if statment that makes the requisite inputs numeric

if ~isa(whsim, 'numeric')
    whsim = str2num(whsim);
    start_sub = str2num(start_sub);
    end_sub = str2num(end_sub);
    dotest = str2num(dotest);
    isHPC = str2num(isHPC);
end

%% SETUP

%To create the logger reference:
LOG = log4m.getLogger(logfile);
LOG.setCommandWindowLevel(LOG.(logging));
LOG.setLogLevel(LOG.OFF);

if ~isHPC
    addpath('..');
    addpath(fullfile(getenv('HOME'), 'Dropbox', 'MATLAButils'));
end

if start_sub > end_sub
    LOG.error('ERROR', 'subject order is reversed');
end

%% Set Simulation
[ opts, dotest] = set_analysis_options(whsim, isHPC, dotest, set_up_directory_structure, LOG);

%% Load the ROI timecourse data and Design
[ dat ] = load_data_and_design( opts, dotest, LOG );

%% Unpack Options and Data

% Options
whsim = opts.whsim;
whs = opts.whs;
K = opts.K;
seed = opts.seed;
Tremove = opts.Tremove;
doconstrained = opts.doconstrained;
prewhiten = opts.prewhiten;
var_log_transform = opts.var_log_transform;
TukN = opts.TukN;
multivariate = opts.multivariate;
roifile = opts.roifile;
designfile = opts.designfile;
combdesignfile = opts.combdesignfile;
runmodels = opts.runmodels;
output_directory = opts.output_directory;
addmotion = opts.addmotion;

% Data
tcn = dat.tcn;
design = dat.design;
designlabels = dat.designlabels;
T = dat.T;
R = dat.R;
motionX = dat.motionX;
scrubX = dat.scrubX;

% get information on this batch of subjects
subject_list = start_sub:end_sub;
NS = length(subject_list);
subjs = dat.subjs(start_sub:end_sub);

LOG.info('INFO', sprintf('\tdotest = %d', dotest));
LOG.debug('DEBUG', sprintf('\tR = %d', R));
LOG.debug('DEBUG', sprintf('\tT = %d', T));

if start_sub > size(tcn, 2)
    LOG.error('ERROR', ['Start subject index is too large. Try using a smaller'...
        'index or not running in test mode']);
end

%% Set random seed
rng( seed );

%% Optimization Options
if doconstrained == 0
    max_iter = 1000;
    optim_opts = optimoptions('fminunc','Algorithm','trust-region',...
        'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
        'FiniteDifferenceType', 'central', ...
        'Diagnostics', 'off', 'MaxIterations', max_iter);
else
    max_iter = 1000;
    
    optim_opts = optimoptions('fmincon','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
        'Diagnostics', 'off', 'MaxIterations', max_iter,...
        'SpecifyConstraintGradient', true);
end

% save filename
if dotest
    filenm = fullfile( output_directory, sprintf( 'paramswm_whs%d_batch_%d_to_%d_test', ...
        int64(whsim), start_sub, end_sub));
    % optim_opts.Display = 'iter';
    optim_opts.Display = 'none';
else
    optim_opts.Display = 'none';
    filenm = fullfile( output_directory, sprintf( 'paramswm_whs%d_batch_%d_to_%d' , ...
        int64(whsim), start_sub, end_sub));
end

%% Determine the designs for each model
% Number of models
M = size( runmodels , 1 );

% Cell array with information for each model
models = cell( M,1 );
badchol = []; % for prewhitening store if any of the pre-whitening matrices can't be cholesky decomposed
for m=1:M
    
    code = runmodels{ m };
    
    % Store the model description
    model = struct();
    model.code = code;
    model.meaneffects = code{ 1 };
    model.vareffects  = code{ 2 };
    model.description = code{ 3 };
    
    LOG.info( 'INFO', sprintf('Creating design for model %d of %d (%s)' , m , M , model.description{1} ));
    
    [ ismem , meancolsnow ] = ismember( model.meaneffects , designlabels );
    if any( ismem == 0 )
        LOG.error( 'ERROR', 'Could not find match to some strings' );
    end
    
    [ ismem , varcolsnow ] = ismember( model.vareffects , designlabels );
    if any( ismem == 0 )
        LOG.error( 'ERROR', 'Could not find match to some strings' );
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
        if var_log_transform
            initsvar( wh2 ) = 0;
        else
            initsvar( wh2 ) = 1;
        end
    end
    
    model.initsmean = initsmean;
    model.initsvar = initsvar;
    model.initsmotion = 0.1;
    
    % output which parameters the current model contains
    str1 = cell( P1 , 1 );
    for i=1:P1
        str1{ i } = sprintf( 'Beta%d (%s)' , i , designlabels{ meancolsnow( i ) } );
        LOG.debug( 'DEBUG', sprintf('\tParameter: %s' , str1{ i } ));
    end
    str2 = cell( P2 , 1 );
    for i=1:P2
        str2{ i } = sprintf( 'Sigma%d (%s)' , i , designlabels{ varcolsnow( i ) } );
        LOG.debug( 'DEBUG', sprintf('\tParameter: %s' , str2{ i } ));
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
    
    models{ m } = model;
end

%% Fit the models
% Store results in cell arrays (matlab complains when using parallization with arrays)
allparams = cell( NS , R , M );
motionparams = cell( NS , R , M );
allbic    = cell( NS , R , M );
alllls    = cell( NS , R , K+1 , M ); % all log likelihoods
allpredm  = cell( NS , R , M );       % all mean predictions at each time point
allpredv  = cell( NS , R , M );       % all variance predictions at each time point
allbadchol = zeros(NS, R, M); 

% Determine train/test partitions
whtrain_sets = cell( 1 , K+1 );
whtest_sets = cell( 1 , K+1 );
if K>0
    cv = cvpartition( T ,'kfold',K);
    for f=1:K
        %         whtrain = find( cv.training( f ));
        %         whtest  = find( cv.test( f ));
        whtrain = cv.training( f );
        whtest  = cv.test( f );
        whtrain_sets{ f } = whtrain;
        whtest_sets{ f } = whtest;
    end
end

% Add one fold with all data used for training AND testing
whtrain_sets{ K+1 } = true(T,1);
whtest_sets{ K+1 } = true(T,1);

% Design matrix. Have to pre-assign design for each subject so we can
% process each subject in parallel
% Loop over subjects

tic;
for s=1:NS
    
    i = subject_list(s);
    
    LOG.info( 'INFO', sprintf('Fitting models ..... working on subject %d of [%d:%d]\n', i, start_sub, end_sub));
    YS = double( squeeze( tcn( : , i , : )));
    
    if addmotion
        motionnow = double(motionX{i});
        scrubnow  = scrubX{i};
    else
        motionnow = [];
        scrubnow  = [];
    end
    
    % Loop over regions/voxels
    for j=1:R
        if rem(j, 10) == 0
            LOG.info( 'INFO', sprintf('\tregion %d\n' , j ));
        end
        YALL = YS(:,j);
        
        %%  FIT MODELS
        
        [ paramsnow, ismotionparam, bicnow, lloutofsample, predm, predv, badchol ] = fit_models_cv(...
            models, design, motionnow, scrubnow, YALL, var_log_transform, ...
            doconstrained, TukN, prewhiten, optim_opts, whtrain_sets, ...
            whtest_sets, LOG, i, j);
        
        if any( badchol ) 
            LOG.warn('WARN', sprintf('S is not positive definite subject %d, region %d', i, j));
        end
        
        % Save Results for Each Fold
        for f = 1:K+1
            for m = 1:M
                
                if (f==K+1)
                    % Only store the parameters and predictions when training data = full data
                    allparams{s,j,m}  = paramsnow{m}(~ismotionparam{m});
                    motionparams{s,j,m} = paramsnow{m}(ismotionparam{m});
                    allbic{s,j,m}     = bicnow{m};
                    
                    % get around scrubbing
                    if addmotion
                        temp_predm = NaN(T,1);
                        temp_predv = NaN(T,1);
                        temp_predm(~scrubnow) = predm{m};
                        temp_predv(~scrubnow) = predv{m};
                    else
                        temp_predm = predm{m};
                        temp_predv = predv{m};
                    end
                    
                    allpredm{s,j,m}   = temp_predm';
                    allpredv{s,j,m}   = temp_predv';
                    allbadchol(s,j,:) = badchol; 
                else
                    % Store the out-of-sample likelihood for the f-th partition
                    alllls{s,j,f,m} = lloutofsample{f,m};
                end
            end
            
        end
    end
end
runtime = toc;

LOG.info('INFO', sprintf('Total Run Time = %2.2f', runtime));

%% Convert cell arrays to arrays
for m=1:M
    % Get number of parameters
    P = length( models{m}.initsmean ) + length(models{m}.initsvar);
    models{ m }.allparams = reshape( cell2mat( allparams(:,:,m)) , NS , P , R );
    models{ m }.motionparams = motionparams(:, :, m );
    
    predsmodelm  = reshape( cell2mat( allpredm(:,:,m))  , NS , T , R );
    predsmodelv  = reshape( cell2mat( allpredv(:,:,m))  , NS , T , R );
    predsmodelm = permute( predsmodelm , [ 2 1 3 ] );
    predsmodelv = permute( predsmodelv , [ 2 1 3 ] );
    models{ m }.allpredsm = single( predsmodelm );
    models{ m }.allpredsv = single( predsmodelv );
end

% Calculate the total out-of-sample log-likelihoods across folds
M = size(alllls, 4);
allllsm = nansum( cell2mat( alllls(:,:,1:K,:)) , 3 );
allllsm = reshape(allllsm, [ NS, R, M]);

allbicm = cell2mat( allbic );

% Determine Best Model for Subject + ROI combination
bestmodelCV  = NaN( NS , R );
bestmodelBIC = NaN( NS , R );
for s=1:NS
    for j=1:R
        [ a , whmin ] = min( allbicm(s,j,:) );
        if isnan(a)
            bestmodelBIC( s , j ) = -1; % set NaN to modelid -1
        else
            bestmodelBIC( s , j ) = whmin;
        end
        if (K > 0)
            [ a , whmax ] = max( allllsm( s,j,:) );
            if isnan(a)
                bestmodelCV( s , j ) = -1; % set NaN to modelid -1
            else
                bestmodelCV( s , j ) = whmax;
            end
        end
    end
end

savedwhsim = whsim;
LOG.info('INFO', sprintf( 'saving file %s', filenm) );
save( filenm , 'allbicm', 'allllsm', 'alllls', 'bestmodelBIC' , 'bestmodelCV' ,...
    'M', 'K' , 'models' , 'savedwhsim', 'runtime', 'subjs', 'start_sub', ...
    'end_sub', 'allbadchol', '-v7.3');

