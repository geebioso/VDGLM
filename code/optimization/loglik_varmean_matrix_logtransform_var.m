function [p,g, H] = loglik_varmean_matrix_logtransform_var( params, Xm, Xv, Y )
%% Compute Negative Log Likelihood of Variance Design Model


% Using a single voxel, it looks like the transform doesn't greatly affect
% the mean parameters. However, all variances get set to 1. 

nbeta = size( Xm , 2 );
nsigma = size( Xv , 2 );
if nbeta+nsigma ~= length( params )
    error( 'Parameter mismatch' );
end

beta  = params(1:nbeta);
sigma = params(nbeta+1:end);

if size(Xm, 2) ~= size(beta, 1)
   beta = beta';  
end

if size(Xv, 2) ~= size(sigma, 1)
    sigma = sigma'; 
end

yd  = Y - sum( Xm*beta , 2 );
yd2 = yd.^2;

vr  = sum( Xv * sigma , 2 );

ps = 0.5*log(2*pi) + 0.5*vr + 0.5*yd2./exp(vr); 
p = sum(ps); 


%% Gradients 
if nargout > 1
    g = zeros(nbeta+nsigma,1); 
        
    % beta 
    pY   = -Xm.*yd./exp(vr);
    
    g(1:nbeta) = sum(pY, 1); 
      
    % sigma
    pY = 0.5*Xv.*( 1 - yd2./exp(vr)); 
   
    g(nbeta+1:nbeta+nsigma) = sum(pY,1);    

end

%% Compute Hessian 
% follows the math at
% https://www.overleaf.com/10517147vbhzwvjtsgzj#/39249751/
% To match code to notation, Q = yd and Z = st
% Keep in mind that signs are flipped because we are minimizing instead of
% maximizing 
if nargout > 2
    
    % Hbeta2 =  Xm'*(Xm./exp(vr)); 
    Hbeta2 =  Xm'*(Xm./repmat(exp(vr), 1, size(Xm, 2))); 
    
    % Hbetasigma =  Xv'*(Xm.*(yd./exp(vr))); 
    Hbetasigma =  Xv'*( Xm.*repmat( yd./exp(vr), 1, size(Xm,2))); 
    
    % Hsigma2 = 0.5*Xv'*( Xv.*(yd2./exp(vr))) ; 
    Hsigma2 = 0.5* Xv'*( Xv.*( repmat( yd2./exp(vr), 1, size(Xv, 2)) )) ; 
    
    H = [ [Hbeta2 Hbetasigma']; [Hbetasigma Hsigma2]]; 
    
end 


