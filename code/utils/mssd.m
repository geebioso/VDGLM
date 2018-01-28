function [mssd_out, v, mssd_over_v, std_mssd_over_v, mean_mssd_over_v] = mssd( x )

% this function computes the meas successive squared difference (MSSD).
% MSSD measures moment to moment changes in signal where the mean may be
% shifting.

% from http://onlinelibrary.wiley.com/store/10.1002/0471667196.ess2635.pub2/asset/ess2635.pdf?v=1&t=jcwju37j&s=76c5fff6ace51dde0f4d4f12c2aecfbb04a37c0a&systemMessage=Please+be+advised+that+we+experienced+an+unexpected+issue+that+occurred+on+Saturday+and+Sunday+January+20th+and+21st+that+caused+the+site+to+be+down+for+an+extended+period+of+time+and+affected+the+ability+of+users+to+access+content+on+Wiley+Online+Library.+This+issue+has+now+been+fully+resolved.++We+apologize+for+any+inconvenience+this+may+have+caused+and+are+working+to+ensure+that+we+can+alert+you+immediately+of+any+unplanned+periods+of+downtime+or+disruption+in+the+future.
% Successive Differences 

N = length(x);

% mssd
x_d = diff(x);
mssd_out = sum(x_d.^2);
mssd_out = mssd_out/(N-1);

allan = mssd_out/2;  % Allan variance

% variance
v = sum( (x - mean(x)).^2);
v = v/N;

% test for a trend 
R = allan/v; 

% var_R = 1/(N+2); % + order( n^{-3}) 
% ER = 1; 
% Z = ( R - 1 )*sqrt( ( N^2 - 1 )/( N - 2 ) ); 
z_alpha = 1.96; 
h = R < 1 + z_alpha/ sqrt( N + 0.5*( 1 + z_alpha^2) ); % there is a trend if true 

mssd_over_v = mssd_out/v; 
std_mssd_over_v = sqrt( 4*(N - 2)/(N^2 - 1) ); 
mean_mssd_over_v = 2; 