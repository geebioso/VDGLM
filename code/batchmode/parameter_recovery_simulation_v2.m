
%% Options

generating_whsim = 26;
whsim = 26;
isHPC = 0;
dotest = 1;
logfile = 'crap.txt';
logging = 'INFO';
start_sub = 1;
end_sub = 5;


%To create the logger reference:
LOG = log4m.getLogger(logfile);
LOG.setCommandWindowLevel(LOG.(logging));
LOG.setLogLevel(LOG.OFF);

if ~isHPC
    addpath('..');
    addpath(fullfile(getenv('HOME'), 'Dropbox', 'MATLAButils'));
end

%% Set Simulation
[ opts, dotest] = set_analysis_options_v2(whsim, isHPC, dotest, LOG);

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
        'SpecifyObjectiveGradient',true, 'CheckGradients', true,...
        'FiniteDifferenceType', 'central', ...
        'Diagnostics', 'on', 'MaxIterations', max_iter, ...
        'Display', 'iter-detailed');
else
    max_iter = 1000;
    
    optim_opts = optimoptions('fmincon','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true, 'CheckGradients', false,...
        'Diagnostics', 'off', 'MaxIterations', max_iter,...
        'SpecifyConstraintGradient', true);
end


%% Load Results
[results_directory] = set_results_directory( isHPC );

[models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, sub_nums] = ...
    load_results(results_directory, generating_whsim, dotest, LOG, 'analyze');

generating_model = models{2}; % take the mean model
simulated_model = models{1}; % corresponding VDGLM

%% Simulate Data From Our Model, Prewhiten, and Recover Parameters

Nsims = 100;
doplots = 0;

meancolsnow = simulated_model.meancols;
varcolsnow = simulated_model.varcols;

P = size(simulated_model.allparams,2);
simulated_model.allparams = zeros(NS, P, R, Nsims);
simulated_model.allpredsm = zeros(T, NS, R, Nsims);
simulated_model.allpredsv = zeros(T, NS, R, Nsims);
YALL = zeros(T, NS, R, Nsims);

% set up design matrices
Xm = double(design(:, meancolsnow)); % if no motion, motionnow will be empty
Xv = double(design(:, varcolsnow ));
true_params = zeros(NS, P, R); 
        
for s = 1:NS
    for r = 1:R
        
        true_beta = double(generating_model.allparams( s, logical(generating_model.ismeanparam), r ))';
        if whsim == 36
             true_sigma = [0; normrnd(-0.3,0.2,length(varcolsnow) - 1, 1)];   
             true_variance = exp(Xv*true_sigma); 
        else
             true_sigma = [1; normrnd(0,exp(0.2),length(varcolsnow) - 1, 1)];
             true_variance = Xv*true_sigma; 
        end
        true_params(s,:,r) = [true_beta; true_sigma]; 
        
        for n = 1:Nsims
            
            % generate data
            if whsim == 36
                Y= normrnd( Xm*true_beta, sqrt(exp(Xv*true_sigma)) );
            else
                Y= normrnd( Xm*true_beta, sqrt(Xv*true_sigma) );
            end
            YALL(:,s,r,n) = Y;
            
            % Initial parameters for VDGLM optimization
            initsnow = [simulated_model.initsmean, ... % mean experiment regressors
                simulated_model.initsvar]'; % variance parameters
            
%             initsnow = [true_beta', ... % mean experiment regressors
%                 true_sigma']'; % variance parameters
            
            initsnow = double(initsnow);
            
            % Define the Likelihoods and Prediction Functions
            if var_log_transform
                fun = @(x) loglik_varmean_matrix_logtransform_var( x,Xm, Xv, Y);
                
                % Define the (negative) out-of-sample likelihood
                funtest = @(x) loglik_varmean_matrix_logtransform_var( x ,Xm, Xv, Y);
                
                % Define the predicted activation function
                predfun = @(x) preds_varmean_matrix_logtransform( x,Xm, Xv);
            else
                fun = @(x) loglik_varmean_matrix_var( x,Xm, Xv, Y);
                
                % Define the (negative) out-of-sample likelihood
                funtest = @(x) loglik_varmean_matrix_var( x ,Xm, Xv, Y);
                
                % Define the predicted activation function
                predfun = @(x) preds_varmean_matrix_var( x,Xm, Xv);
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
                    nonlcon = @(x) varconstraint2( x, Xm, Xv);
                    
                    % Run the constrained optimizer
                    [depparams,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,optim_opts);
                end
            else
                depparams = nan(size(initsnow));
            end
            
            simulated_model.allparams(s,:,r,n) = depparams; 
            
        end
    end
    
end

%% Plot Parameter Recovery

param_differences = zeros(size(simulated_model.allparams)); 
for n = 1:Nsims
    param_differences(:,:,:,n) = true_params - simulated_model.allparams(:,:,:,n); 
end

for p = 1:P
    
    figure(p); clf; 
    subplot(1,2,1); 
    histogram(param_differences(:,p,:,:)); 
    title( 'Parameter Differences (true - inferred)'); 
    subplot(1,2,2); 
    temp1 = repmat(true_params(:,p,:), [1 1 1 Nsims]); 
    temp2 = simulated_model.allparams(:,p,:,:); 
    scatter( temp1(:), temp2(:));
    refline(1,0); 
    xlabel('true');
    ylabel('inferred');
    
    title( sprintf( 'True vs. inferred %s', simulated_model.paramlabels{p})); 
end





