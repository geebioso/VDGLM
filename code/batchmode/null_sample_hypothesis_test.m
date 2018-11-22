function [] = null_sample_hypothesis_test(whsim, dotest, isHPC, start_sub, end_sub, logging, logfile, set_up_directory_structure, Nsamp)

% This function generated data from the mean model, adds autocorrelation
% similar to the empirical autocorrelation, and then estimates the VDGLM
% parameters. The goal is to find the false positive rate from null data

% INPUT:
%   numeric whsim: which simulation
%   bool dotest(0/1): run in test mode?
%   bool isHPC (0/1): run on HPC or personal computer
%   numeric start_sub: start subject
%   numeric end_sub: end subject
%   numeric Nsamp: number of samples to draw
%   char logging: type of logging e.g, {'ALL','TRACE','DEBUG','INFO','WARN','ERROR','FATAL','OFF'}
%   char logfile: log file, useless argument because I have set the loglevel to off, but log4.m requires it

% NOTE: running this on the HPC will change the input type to char, hence
%   the following if statment that makes the requisite inputs numeric

% TIMING: 
% 449.44s for test mode with 100 samples on my laptop 
% 52.49s for test mode with 10 samples on my laptop 
% 7.94s for analyzedata_batch_v2.m for one subject during test on my laptop
%   with 10 samples 
%    => expect 65 times slower execution 

% 4088.10s for analyzedata_batch_v2.m on the HPC for 1 subject 
%   =>  228,930 seconds = 64 hours for 1 subject on the HPC with 100 samples 
%   =>  32 hours for 1 subject with 50 samples 

% to get resource allocation on HPC: 
% (/usr/bin/time -v ./null_sample_hypothesis_test 26 1 1 1 1 50 INFO test_log.txt) &> null_sample_resource_usage.txt
% (/usr/bin/time bash test_time.sh) &> test.txt

% Make all inputs numeric
if ~isa(whsim, 'numeric')
    whsim = str2num(whsim);
    start_sub = str2num(start_sub);
    end_sub = str2num(end_sub);
    dotest = str2num(dotest);
    isHPC = str2num(isHPC);
    Nsamp = str2num(Nsamp);
end

if ~isHPC
    addpath('..');
    addpath(fullfile(getenv('HOME'), 'Dropbox', 'MATLAButils'));
end

LOG = log4m.getLogger(logfile);
LOG.setCommandWindowLevel(LOG.(logging));
LOG.setLogLevel(LOG.OFF);

if start_sub > end_sub
    LOG.error('ERROR', 'subject order is messed up');
end

%% Set Simulation
[ opts, dotest] = set_analysis_options(whsim, isHPC, dotest, set_up_directory_structure, LOG);

%% Load the ROI timecourse data and Design
[ dat ] = load_data_and_design( opts, dotest, LOG);

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

% append null to the start of the output_directory
pth = strsplit(output_directory, filesep);
pth{end} = ['null_' pth{end}];
output_directory = join(pth, filesep);
output_directory = output_directory{1};

% Data
tcn = dat.tcn;
design = dat.design;
designlabels = dat.designlabels;
T = dat.T;
R = dat.R;
motionX = dat.motionX;
scrubX = dat.scrubX;

subject_list = start_sub:end_sub;
NS = length(subject_list);
subjs = dat.subjs(start_sub:end_sub);

LOG.info('INFO', sprintf('\tdotest = %d', dotest));
LOG.debug('DEBUG', sprintf('\tR = %d', R));
LOG.debug('DEBUG', sprintf('\tT = %d', T));

%% Set random seed
rng( seed );

%% Optimization Options
if doconstrained == 0
    max_iter = 1000;
      optim_opts = optimoptions('fminunc','Algorithm','sqp',...
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
% Temporary
%optim_opts.Display = 'iter';


% save filename
if dotest
    filenm = fullfile( output_directory, sprintf( 'null_sample_whs%d_batch_%d_to_%d_test', ...
        int64(whsim), start_sub, end_sub));
    % optim_opts.Display = 'iter';
    optim_opts.Display = 'none';
else
    optim_opts.Display = 'none';
    filenm = fullfile( output_directory, sprintf( 'null_sample_whs%d_batch_%d_to_%d' , ...
        int64(whsim), start_sub, end_sub));
end

%% Determine the designs for each model
% We only want storage for the mean and var+mean models. Mean for data
% generation and Var+Mean for inference
to_run = {'Var+Mean', 'Mean'};
M = size( runmodels , 1 );
to_keep = false(M,1);
for m = 1:M
    if ismember( runmodels{m}{3}{1}, to_run)
        to_keep(m) = true;
    end
end
runmodels = runmodels(to_keep);
M = size( runmodels, 1 );

% Cell array with information for each model
models = cell( M,1 );
badchol = []; % for prewhitening: store if any of the pre-whitening matrices can't be cholesky decomposed
for m=1:M
    
    code = runmodels{ m };
    
    % Store the model description
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
        initsvar( wh2 ) = 1;
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
allparams = cell( NS , R , M, Nsamp );
motionparams = cell( NS , R , M, Nsamp );
allbic    = cell( NS , R , M, Nsamp );
alllls    = cell( NS , R , M , Nsamp, K+1 ); % all log likelihoods
allbadchol = zeros(NS, R, M); 

% Determine train/test partitions
whtrain_sets = cell( 1 , K+1 );
whtest_sets = cell( 1 , K+1 );
if K>0
    cv = cvpartition( T ,'kfold',K);
    for f=1:K
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
    
    % motion regressors
    if addmotion
        motionnow = double(motionX{i});
        scrubnow  = scrubX{i};
    else
        motionnow = [];
        scrubnow  = [];
    end
    
    % create design matrices
    % we are only running the VDGLM model, so we can hard code this by
    % subject (compared to analyzedata_batch_v2.m where Xm must be flexible to
    % account for each model)
    meancols =  models{1}.meancols;
    varcols = models{1}.varcols;
    Xm = [design( :, meancols ), motionnow];
    Xv = design( :, varcols );
    
    ismotionparam = cell(2,1);
    ismotionparam{1} = [false(1, length(meancols)), true(1, size(motionnow,2)), false(1, length(varcols))];
    ismotionparam{2} = [false(1, length(meancols)), true(1, size(motionnow,2)), false(1, 1) ];
    
    %% Scrubbing
    for f = 1:K+1
        whtrain_sets{f}(scrubnow) = [];
        whtest_sets{f}(scrubnow)  = [];
    end
    
    YS(scrubnow,:) = [];
    Xm(scrubnow,:) = [];
    Xv(scrubnow,:) = [];
    motionnow(scrubnow, :) = [];
    
    % Loop over regions/voxels
    for j=1:R
        if rem(j, 10) == 0
            LOG.info( 'INFO', sprintf('\tregion %d\n' , j ));
        end
        Y = YS(:,j);
        
        [ paramsnow, ismotionparam, bicnow, lloutofsample, badchol ] = fit_models_cv_null_sample(...
            models, Xm, ismotionparam, Xv, Y, var_log_transform, doconstrained, TukN, ...
            prewhiten, optim_opts, whtrain_sets, whtest_sets, LOG, i, j, Nsamp);
        
        allbadchol(s,j,:) = badchol; 
        
        for f=1:K+1
            for m = 1:M
                for n = 1:Nsamp
                    if (f==K+1)
                        % Only store the parameters and predictions when training data = full data
                        allparams{s,j,m,n}  = paramsnow{m, n}(~ismotionparam{m});
                        motionparams{s,j,m,n} = paramsnow{m, n}(ismotionparam{m});
                        allbic{s,j,m,n}     = bicnow{m, n};
                    else
                        % Store the out-of-sample likelihood for the f-th partition
                        alllls{s,j,m,n,f} = lloutofsample{f,m, n};
                    end
                end
            end
        end
    end
end
runtime = toc;

LOG.info('INFO', sprintf('Total Run Time = %2.2f\n', runtime));

%% Convert cell arrays to arrays
for m=1:M
    
    % Get number of parameters
    P = length(models{m}.initsmean) +  length(models{m}.initsvar);
    models{ m }.allparams = zeros( NS, R, Nsamp, P );
    models{ m }.motionparams = cell( NS, R, Nsamp );
    
    for s = 1:NS
        for r = 1:R
            for k = 1:Nsamp
                models{ m }.allparams(s,r,k,:) = allparams{s,r,m,k};
                models{ m }.motionparams = motionparams(:, :, m );
            end
        end
    end
    
end

% model comparison measures
allllsm = zeros(NS, R, Nsamp, M);
allbicm = zeros(NS, R, Nsamp, M);
for s = 1:NS
    for r = 1:R
        for k = 1:Nsamp
            for m = 1:M
                allllsm(s,r,:,m) = sum( cell2mat(alllls(s,r,m,1:K)), 4 );
                allbicm(s,r,:,m) = allbic{s,r,m};
            end
        end
    end
end

% Determine Best Model for Subject + ROI combination
bestmodelCV  = NaN( NS , R, Nsamp);
bestmodelBIC = NaN( NS , R, Nsamp );
for s=1:NS
    for j=1:R
        for k = 1:Nsamp
            [ a , whmin ] = min( allbicm(s,j,k,:) );
            if isnan(a)
                bestmodelBIC( s , j , k ) = -1; % set NaN to modelid -1
            else
                bestmodelBIC( s , j , k ) = whmin;
            end
            if (K > 0)
                [ a , whmax ] = max( allllsm( s,j,k,:) );
                if isnan(a)
                    bestmodelCV( s , j, k ) = -1; % set NaN to modelid -1
                else
                    bestmodelCV( s , j, k ) = whmax;
                end
            end
        end
    end
end

save( filenm , 'M', 'K' , 'models' , 'runtime', 'subjs', 'allbicm', 'allllsm', ...
    'bestmodelCV', 'bestmodelBIC', 'start_sub', 'end_sub', 'allbadchol', '-v7.3');


