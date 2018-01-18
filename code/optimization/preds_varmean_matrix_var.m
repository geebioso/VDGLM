function [ predsm , predsv ] = preds_varmean_matrix_var( params, Xm, Xv )
%% Compute Predicted Means and Stds of Variance Design Model
nbeta = size( Xm , 2 );
nsigma = size( Xv , 2 );
if nbeta+nsigma ~= length( params )
    error( 'Parameter mismatch' );
end

beta  = params(1:nbeta);
sigma = params(nbeta+1:end);

predsm  = sum( Xm*beta , 2 );
% predicts the variance
predsv  = sum( Xv*sigma , 2 );





