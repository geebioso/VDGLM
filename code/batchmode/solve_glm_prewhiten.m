function [Y_pre, Xm_pre, B_pre, sigma2_pre, B, sigma2, L, badchol ] = solve_glm_prewhiten( Xm, Y, TukN )

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
%   numeric TukN: the size of the Tukey Taper 

% OUTPUT: 
%   numeric Y_pre: prewhitened data
%   numeric X_pre: prewhitened design matrix
%   numeric B_pre: prewhitened OLS solution
%   numeric sigma2_pre: prewhitened variance estimate 
%   numeric B: OLS solution
%   numeric sigma2: variance estimate 
%   numeric L: estimated autocovariance 
%   bool badchol: is there a problem when we compute the cholesky decomposition
%       to get the square root of the autocovariance matrix? 


%% 1) Solve GLM 
[B, sigma2, ~] = solve_glm( Xm, Y ); 
res = Y - Xm * B; % compute residuals for prewhitening 

%% 2) Residual Autocorrelation Estimation and Correction
T = length(res); 

% compute raw autocorrelation
if TukN < 0
    TukN = round(2*sqrt(T));
end
A = autocorr_woolrich(res,TukN -1 );

% Tukey Taper
%power = (1:T)';
lag = (0:T-1)';
Tuk = 0.5*( 1 + cos( (pi*lag)/TukN ));
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
Y_pre = L\Y;
Xm_pre = L\Xm;

%% 4) Compute GLM Parameters (with pre-whitened data)
[B_pre, sigma2_pre, ~] = solve_glm( Xm_pre, Y_pre ); 

