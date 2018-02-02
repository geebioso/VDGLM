function [ A ] = autocorr_woolrich( res, M )

% computes the autocorrelation of the residuals according to eq. 11 in 
% Temporal Autocorrelation in Univariate Linear Modeling of FMRI Data
% Woolrich et al. 2001 

% INPUT: 
%   res: model residuals 
%   M: largest lag at which to compute autocorrelation 

% OUTPUT: 
%   A: autocorrelation at each lag 

varhat = var(res, 1); 
resbar = mean(res);
N = length(res); 

A = zeros(M,1); 
for m = 0:M
    
    t1 = res(1:(N-m)) - resbar; 
    t2 = res(m+1:N) - resbar; 
    % A(m+1) = (1/varhat)*sum(t1.*t2)/(N-m); % unbiased
    A(m+1) = (1/varhat)*sum(t1.*t2)/(N);  % biased 
end

end

