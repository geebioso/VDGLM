function [p,g,H] = loglik_varmean_matrix_var( params, Xm, Xv, Y )
%% Compute Negative Log Likelihood of Variance Design Model
nbeta = size( Xm , 2 );
nsigma = size( Xv , 2 );
if nbeta+nsigma ~= length( params )
    error( 'Parameter mismatch' );
end

beta  = params(1:nbeta);
sigma = params(nbeta+1:end);

% if sum( sigma ) <= 0
%     warning( 'Sigma < 0' ); 
% end

yd  = Y - sum( Xm*beta , 2 );
yd2 = yd.^2;

vr  = sum( Xv * sigma , 2 );
vr( vr < 0 ) = 0.000001; % this constraint is violated

vr2 = vr.^2;
vr3 = vr.^3;

ps = 0.5*log(2*pi) + 0.5*log( vr ) + yd2./(2*vr); 
p = sum(ps); 

% if ~isreal( p )
%     %warning( 'Complex number' );
% end

%% Gradients 
if nargout > 1
    g = zeros(nbeta+nsigma,1); 
        
    % beta 
    % had to change formula to be compatible with 2016a 
    % pY   = -(Xm.*yd)./vr;  % 
    pY   = -(Xm.*repmat(yd, 1, size(Xm,2)))./repmat(vr, 1, size(Xm, 2)); 
    g(1:nbeta) = sum(pY,1);
       
    % sigma 
    % pY = Xv./(2*vr) - Xv.*yd2./(2*vr2); 
    pY = Xv./(2*repmat(vr, 1, size(Xv, 2))) - Xv.*repmat(yd2, 1, size(Xv, 2))./(2*repmat(vr2, 1, size(Xv,2))); 
    g(nbeta+1:nbeta+nsigma) = sum(pY,1);    
end

%% Compute Hessian 
% follows the math at
% https://www.sharelatex.com/project/59e1320f330afa2e57e9501f
% To match code to notation, Q = yd and Z = st
if nargout > 2
    
    % Hbeta2 =  -Xm'*(Xm./vr); 
    Hbeta2 =  Xm'*(Xm./repmat(vr, 1, size(Xm, 2))); 
    
    % Hbetasigma =  Xv'*( Xm.*yd./vr2); 
    Hbetasigma =  Xv'*( Xm.*repmat( yd./vr2, 1, size(Xm,2))); 
    
    % Hsigma2 = - Xv'*( Xv.*( (1./2*vr2) - yd2./vr2 )) ; 
    Hsigma2 = - Xv'*( Xv.*( repmat( 1./(2*vr2) - yd2./vr3, 1, size(Xv, 2)) )) ; 
    
    H = [ [Hbeta2 Hbetasigma']; [Hbetasigma Hsigma2]]; 
    
    % 
end 



