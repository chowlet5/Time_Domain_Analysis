%--------------------------------------------------------------------------
%Given the CDF p this function evaluates inverse of the standard normal CDF
%--------------------------------------------------------------------------
function phi_inv = isnc(p)

phi_inv = sqrt(2)*erfinv(2*p-1);  

%This can be equivalently done by using inverse of the comprelementary
%error function as follows:

% phi_inv = -sqrt(2)*erfcinv(2*p);  