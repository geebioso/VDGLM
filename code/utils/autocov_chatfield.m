function [ A ] = autocov_chatfield( res, M )

% computes the autocovariance of the residuals according to eq. 4.1 in 
% Analyis of Time Series An Introdution Chatfield 1996 

% INPUT: 
%   res: model residuals 
%   M: largest lag at which to compute autocorrelation 

% OUTPUT: 
%   A: autocovariance at each lag 

varhat = var(res, 1); 
resbar = mean(res);
N = length(res); 

A = zeros(M,1); 
for m = 0:M
    
    t1 = res(1:(N-m)) - resbar; 
    t2 = res(m+1:N) - resbar; 
    A(m+1) = sum(t1.*t2)/N;  % Chatfield 1996 
    % A(m+1) = sum(t1.*t2)/(N - m);  % Jenkins and Watts 1968, smaller bias but higher MSE 
    
end

end

