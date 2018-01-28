function [Y, X, B, sigma2, L, badchol ] = solve_glm( X, Y, prewhiten, varargin) 
% TukN, prewhiten)

% This function returns the prewhitened OLS solution of the model
%   Y = Xm*B + e
% fit using only the data contained in trainnow. Prewhitening proceeds as
% follows:
%   1) compute the GLM OLS solution
%   2) estimate the autocovariance of the residuals
%   3) prewhiten using the autocovariance estimate (via Woolrich et al. 2001)
%   4) compute the prewhitened GLM OLS solution

% INPUT:
%   numeric Xm: the design matrix
%   numeric Y: the BOLD response
%   bool prewhiten: whether to prewhiten or not 
%   numeric varargin{1}: the size of the Tukey Taper

% OUTPUT:
%   numeric Y: response (prewhitened if varagin{2} is true) 
%   numeric X: design matrix (prewhitened if varagin{2} is true) 
%   numeric B: OLS mean solution (prewhitened if varargin{2} is true) 
%   numeric sigma2: OLS variance solution (prewhitened if varargin{2} is true) 
%   numeric L: estimated autocovariance if varargin{2} is true (else identity) 
%   bool badchol: is there a problem when we compute the cholesky decomposition
%       to get the square root of the autocovariance matrix?

if nargin == 4
    TukN = varargin{1}; 
end

%% 1) Solve GLM
T = size( X, 1);
P = size( X,2);

B = (X'*X) \ X' * Y;
res = Y - X*B;
df = T - P;

sigma2 = sum( res.^2 )/df;

res = Y - X * B; % compute residuals for prewhitening

if prewhiten
    %% 2) Residual Autocorrelation Estimation and Correction
    % compute raw autocorrelation
    if TukN < 0
        TukN = round(2*sqrt(T));
    end
    A = autocorr_woolrich(res,TukN -1 );
    
    % Tukey Taper
    %power = (1:T)';
    lag = (0:T-1)';
    Tuk = 0.5*( 1 + cosd( (pi*lag)/TukN ));
    %Tuk(1:TukN) = Tuk(1:TukN).*A;
    Tuk(1:TukN) = Tuk(1:TukN).*A;
    Tuk((TukN+1):end) = 0; % assumption from paper
    
    % compute sample autocovariance from Tukey Taper
    autocov = Tuk*var(res);
    S = toeplitz( autocov );
    
    % compute cholesky decomposition, if unable to compute
    % decomposition, don't prewhiten ( L = eye(size(S)) )
    [L, p] = chol(S, 'lower');
    if p ~= 0
        badchol = 1;
        L = eye(size(S));
    else
        badchol = 0;
    end
    
    %% 3) Pre-whitening
    Y = L\Y;
    X = L\X;
    
    %% 4) Compute GLM Parameters (with pre-whitened data)
    B = (X'*X) \ X' * Y;
    res = Y - X*B;
    sigma2 = sum( res.^2 )/df;

else % don't prewhiten 
    L = eye(T); 
    badchol = 0; 
    
end



