function [p,g] = loglik_varmean_matrix_2( params, Xm, Xv, Y )
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

% yd  = Y - sum( Xm.*beta , 2 );
yd = sum(Y - Xm*beta, 2); 
yd2 = yd.^2;

% st  = sum( Xv .* sigma , 2 );
st  = sum( Xv * sigma , 2 );
st( st < 0 ) = 0.000001;

st2 = st.^2;
st3 = st.^3;

ps = 0.5*log(2*pi) + log( st ) + yd2./(2*st2); 
p = sum(ps); 

% if ~isreal( p )
%     %warning( 'Complex number' );
% end

%% Gradients 
if nargout > 1
    g = zeros(nbeta+nsigma,1); 
        
    % beta 
    pY   = -(Xm.*yd)./st2; 
    g(1:nbeta) = sum(pY,1);
       
    % sigma 
    pY = Xv./st - Xv.*yd2./st3; 
    g(nbeta+1:nbeta+nsigma) = sum(pY,1);    
end



