function [p,g] = loglik_varmean_matrix_logtransform( params, Xm, Xv, Y )
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

st  = sum( Xv * sigma , 2 );

ps = 0.5*log(2*pi) + 0.5*st + 0.5*yd2./exp(st).^2; 
p = sum(ps); 


%% Gradients 
if nargout > 1
    g = zeros(nbeta+nsigma,1); 
        
    % beta 
    pY   = -0.5*Xm.*yd2./exp(st).^2;
    
    g(1:nbeta) = sum(pY, 1); 
      
    % sigma
    pY = Xv - Xv.*yd2./exp(st).^2; 
    g(nbeta+1:nbeta+nsigma) = sum(pY,1);    

end



