%function h = toepsolve(r,q);
%TOEPSOLVE  Solve Toeplitz system of equations.
%    Solves R*h = q, where R is the symmetric Toeplitz matrix
%    whos first column is r
%    Assumes all inputs are real
%    Inputs:
%       r - first column of Toeplitz matrix, length n
%       q - rhs vector, length n
%    Outputs:
%       h - length n solution
%
%   Algorithm from Roberts & Mullis, p.233
%
%   Author: T. Krauss, Sept 10, 1997