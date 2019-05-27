%--------------------------------------------------------------------------
%Given the data x this function evaluates CDF of the standard nomal distribution
%--------------------------------------------------------------------------
function p = cdf_std_norm(x)

p = 0.5*(1+erf(x/sqrt(2)));  

%This can be done equivalently using the compementary error
%function as follows

% p = 0.5*erfc(-x/sqrt(2));