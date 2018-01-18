function [B, sigma2, df] = solve_glm(X, Y)

% This function returns the OLS solution of the model
%   Y = Xm*B + e, e ~ N(0, sigma^2)

% INPUT:
%   Xm: the design matrix
%   Y: the BOLD response

% OUTPUT:
%   B: OLS solution
%   sigma2: variance estimate
%   df: degrees of freedom

%%
T = size( X, 1);
P = size(X,2);

% solve
B = (X'*X) \ X' * Y;
res = Y - X*B;
df = T - P;

sigma2 = sum( res.^2 )/df;

end