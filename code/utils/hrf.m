%% Haemodynamic Response Function
function [h] = hrf(alpha, beta, const, A, t)
term1 = (t.^(alpha(1)-1)*beta(1)^(alpha(1)).*exp(-beta(1)*t))/gamma(alpha(1));

term2 = (const*t.^(alpha(2)-1)*beta(2)^alpha(2).*exp(-beta(2)*t))/gamma(alpha(2));
h = A*(term1 - term2);
end


