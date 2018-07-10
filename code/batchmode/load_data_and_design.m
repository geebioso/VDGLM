 function [ dat ] = load_data_and_design( opts, dotest, LOG)


 % 12/21/17 GG: changed individual design matrices to use a single design matrix 
 
%% Unpack Options 

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
addmotion = opts.addmotion; 

%% Load the ROI timecourse data if not in memory already
LOG.info('INFO', sprintf('Loading data: %s' , roifile ));
load( roifile );

% if (whs == 8)
%     % For the Ohio State Well-Being data set, extract the WM task data for only the first 174 subjects
%     tc = tc( 1:174 , 9 );
% else
% if (whs == 9)
%     tc = tc(1:174, [5,9]); % extract the resting state and WM task data
% end

%% Create some constants
NS  = size( tc , 1 );     % number of subjects
R   = size( tc{1} , 2 );  % number of ROIs

%% Convert time-course data to array
LOG.info('INFO', 'Converting data and transforming to z-scores' );

if whs == 2 % then we are using HCP data and have run 1 and run 2 timecourses
    tc = tc(:,1);
     motionX = motionX(:,1); 
end

if whs==9
    tcn1 = cell2mat(tc(:,1)); % RS
    tcn2 = cell2mat(tc(:,2)); % WM
    
    T1 = size(tc{1,1}, 1);
    T2 = size(tc{1,2},1);
    
    tcn1 = reshape( tcn1, T1, NS, R);
    tcn2 = reshape( tcn2 , T2 , NS , R ); % Time x (Subject/ROI)
    
    tcn = [tcn1; tcn2];
    T = size(tcn, 1);
else
    T   = size( tc{1} , 1 );  % length of time series
    tcn = cell2mat( tc );
    tcn = reshape( tcn , T , NS , R ); % Time x (Subject/ROI)
end

%% Remove some initial time steps if needed
if Tremove>0
    tcn( 1:Tremove , : , : ) = [];
    T = size( tcn , 1 );
    
    for s = 1:NS
       motionX{s}( 1:Tremove , :) = [];  
    end
end

%% Split Motion and Scrubbing
if whs == 2
    tempX = cell(NS,1);
    scrubX = cell(NS,1);
    for s = 1:NS
        if s == 11
           x = 1;  
        end
        tempX{s} = motionX{s}(:,1:6);
        if size(motionX{s},2) > 6
            scrubX{s} = false(T,1); 
            for c = 7:size(motionX{s}, 2)
                temp = false(T,1); 
                [~, ii] = max( motionX{s}(:,c) ); 
               temp(ii) = true; 
               scrubX{s} = or( scrubX{s}, temp); 
            end
        else
            scrubX{s} = false(T, 1); 
        end
    end
    motionX = tempX; 
end
%% Convert each ROI / Subject data to Z-scores
stdnow = std( tcn , [] , 1 );
stdnow( stdnow == 0 ) = 1;
% rewrot:  tcn = ( tcn - mean( tcn , 1 )) ./ stdnow; to be compatible with 2016a for HPC
tcn = bsxfun( @rdivide, bsxfun( @minus, tcn, mean(tcn,1)), stdnow);



%% Create the convolved design
LOG.info('INFO', 'Creating convolved design matrix for all experimental variables' );
switch whs
    case {8,9} % OSU
        % Ohio State Well-Being
        load( roifile, 'mtx');
        
        design = mtx.WorkingMem;
        temp_d = diff(design(:,3));
        temp_d = ([temp_d; 0] + [0; temp_d]) / 2;
        design = [ones(T, 1), design, temp_d];
        
        design( 1:Tremove , : ) = [];
        
        designlabels = { 'Intercept' , 'Underlined' , '2-back' , 'Instruction', 'Fixation', 'T_Instruction'};
        % add resting state indicator
        if whs == 9
            rs = zeros(T, 1);
            rs(1:T1) = 1;
            
            % convolve
            hrf = spm_hrf(2);
            rs_conv = conv( rs, hrf, 'same');
            
            design  = [ zeros(T1, size(design, 2)); design];
            design  = [design rs_conv];
            
            designlabels = [designlabels, {'Rest'}];
        end
        
    case {0,1,2} % HCP
        
        % load(designfile); % loads rst
        design_dat = load(combdesignfile); % loads combined design
        
        rst = design_dat.rst; 
        design = [ ones(size(rst.mtx, 1), 1), rst.mtx];
        
        % compute temporal derivative for the fixation condition
        temp_d = diff(design(:,4)); % temporal derivative, one point less than X
        temp_d = ([temp_d; 0] + [0; temp_d]) / 2;
        
        design = [design temp_d];
        design(1:Tremove, :) =[];
        
        designlabels = { 'Intercept' , '0-back' , '2-back' , 'Instruction', 'Fixation', 'T_Instruction'};
end

%% Remove Subject 67 From HCP Data (HEAD MOTION)
switch whs
    case {0,1,2}
        if ~(whs == 2)
            tcn(:, 67, :) = [];
            NS = NS - 1;
        end
end

%% Multivariate

switch multivariate
    case 1
        %% Transform Data
        tcn = reshape(tcn, T*NS, 1, R); % tcn is
        %         isequal(tcn( T+1:2*T , 1, 2), tcn(:, 2,2))
        %
        %         ans =
        %
        %         logical
        %
        %         1
        % this means the second time series for the second region corresponds to
        % the time series from the second subject for the second region
        
        %% Transform Design
        design = repmat(design, NS, 1);
        
        %% Constants
        [T, NS, R] = size(tcn);
        
end

%% TEST
if dotest
    NS = 5;
    R = 3;
    tcn = tcn(:,1:NS, 1:R);
end

%% Load Into Data Structure 
dat = struct(); 
dat.tcn = tcn;  
dat.design = design; 
dat.designlabels = designlabels; 
dat.T = T; 
dat.NS = NS; 
dat.R = R; 
dat.subjs = subjs; 
if addmotion
    dat.motionX = motionX; 
    dat.scrubX = scrubX; 
else
    dat.motionX = NaN;  
    dat.scrubX = NaN; 
end


