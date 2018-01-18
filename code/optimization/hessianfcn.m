function Hout = hessianfcn( params, Xm, Xv, Y, lambda)

nbeta = size( Xm , 2 );
nsigma = size( Xv , 2 );
T = size(Xv,1); 
if nbeta+nsigma ~= length( params )
    error( 'Parameter mismatch' );
end

beta  = params(1:nbeta);
sigma = params(nbeta+1:end);


yd  = Y - sum( Xm*beta , 2 );
yd2 = yd.^2;

vr  = sum( Xv * sigma , 2 );
vr( vr < 0 ) = 0.000001; % this constraint is violated

vr2 = vr.^2;


%% Hessian of objective

% Hbeta2 =  -Xm'*(Xm./vr); 
    Hbeta2 =  -Xm'*(Xm./repmat(vr, 1, size(Xm, 2))); 
    
    % Hbetasigma =  2*Xv'*( Xm.*yd./vr2); 
    Hbetasigma =  2*Xv'*( Xm.*repmat( yd./vr2, 1, size(Xm,2))); 
    
    % Hsigma2 = - Xv'*( Xv.*( (1./2*vr2) - yd2./vr2 )) ; 
    %Hsigma2 = - Xv'*( Xv.*( repmat( 1./2*vr2 - yd2./vr2, 1, size(Xv, 2)) )) ; 
    Hsigma2 =  Xv'*( Xv.*( repmat( 1./2*vr2 - yd2./vr2, 1, size(Xv, 2)) )) ; 
    
    H = [ [Hbeta2 Hbetasigma']; [Hbetasigma Hsigma2]]; 
    
%% Hessian of nonlinear inequality constraint
Hg = zeros(T, 1);  

%% Total Hessian 
Hout = H + lambda.ineqnonlin'*Hg;