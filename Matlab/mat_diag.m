%--------------------------------------------------------------------------
% This function is used to create a diagonal matrix of size nxn from a row
% vector of nx1
%--------------------------------------------------------------------------
function[B]=mat_diag(A,n)
for i=1:n
    B(i,i)=A(i,1);
end