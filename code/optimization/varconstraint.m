function [c,ceq] = varconstraint(params, Xm, Xv)
nbeta = size( Xm , 2 );
nsigma = size( Xv , 2 );
if nbeta+nsigma ~= length( params )
    error( 'Parameter mismatch' );
end

beta  = params(1:nbeta);
sigma = params(nbeta+1:end);

predsv  = sum( Xv * sigma , 2 );
count = sum( predsv < 0 );

% c = count;
% ceq = [];


c = []; 
ceq = count; 
