function [ predsm , predsv ] = preds_varmean_matrix_logtransform( params, Xm, Xv )
%% Compute Predicted Means and Stds of Variance Design Model
nbeta = size( Xm , 2 );
nsigma = size( Xv , 2 );
if nbeta+nsigma ~= length( params )
    error( 'Parameter mismatch' );
end

beta  = params(1:nbeta);
sigma = params(nbeta+1:end);

if size(Xm, 2) ~= size(beta, 1) % added this and switched from pointwise multiplication (only works in r2017) 
   beta = beta';  
end

if size(Xv, 2) ~= size(sigma, 1)
    sigma = sigma'; 
end


% predict the mean 
predsm  = Xm*beta; 
% predict the variance
predsv  = exp( Xv*sigma ); 


% Mark's old functions 
% predsv  = sum( exp( Xv*sigma ), 2 ); 
% predsm  = sum( Xm*beta , 2 ); 



