function [c,ceq, DC, DCeq] = varconstraint2(params, Xm, Xv)

% uses an inequality constraint rather than an equality constraint 
% uses the actual variance values rather than the count of whether they are
% greater than zero 

nbeta = size( Xm , 2 );
nsigma = size( Xv , 2 );
T = size(Xv, 1); 
if nbeta+nsigma ~= length( params )
    error( 'Parameter mismatch' );
end

beta  = params(1:nbeta);
sigma = params(nbeta+1:end);

predsv  = sum( Xv * sigma , 2 );
% count = sum( predsv < 0 );

c = -predsv; %  c <= 0 implies that -predsv <= 0 => predsv >= 0 
ceq = [];

% Gradient of the constraints:
if nargout > 2
    DC = -[zeros(nbeta, T) ;Xv'];
    DCeq= []; 
    
end