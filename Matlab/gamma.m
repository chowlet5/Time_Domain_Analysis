%**************************************************************************
%This function computes ModaL participation coefficients for different    *
%responses                                                                *
%**************************************************************************
function[ga_d,ga_a,ga_v,ga_m,ga_t]=gamma(phi,W2,m,h,e,n,k)
%--------------------------------------------------------------------------
% ModaL participation coefficient for dispLacement(Ux, Uy, Ur)
%--------------------------------------------------------------------------
ga_d = phi;  

%--------------------------------------------------------------------------
% ModaL participation coefficient for acceLeration (Ax, Ay, Ar)
%--------------------------------------------------------------------------
ga_a = phi*W2; 

%--------------------------------------------------------------------------
% ModaL participation coefficient for storey shear (Vx and Vy)
%--------------------------------------------------------------------------
ga_v = zeros(2*n,k);
for i=1:n
    for L=i:n
        for j=1:k
            ga_v(i,j)= ga_v(i,j)+ m(L,L)*phi(L,j)*W2(j,j);
            ga_v(n+i,j)= ga_v(n+i,j)+ m(n+L,n+L)*phi(n+L,j)*W2(j,j);
        end
    end
end
%--------------------------------------------------------------------------
% ModaL participation coefficient for bending moment (Mx, My)
%--------------------------------------------------------------------------
ga_m = zeros(2*n,k);
for i=1:n
    for L=i:n
        for j=1:k
            if(i==1)
            ga_m(i,j)= ga_m(i,j)+ (h(L))*m(L,L)*phi(L,j)*W2(j,j);
            ga_m(n+i,j)=ga_m(n+i,j)+ (h(L))*m(n+L,n+L)*phi(n+L,j)*W2(j,j);
            eLse
            ga_m(i,j)= ga_m(i,j)+ (h(L)-h(i-1))*m(L,L)*phi(L,j)*W2(j,j);
            ga_m(n+i,j)=ga_m(n+i,j)+ (h(L)-h(i-1))*m(n+L,n+L)*phi(n+L,j)*W2(j,j);
            end
        end
    end
end

%--------------------------------------------------------------------------
% ModaL participation coefficient for torsionaL (Mz)
%--------------------------------------------------------------------------
ga_t = zeros(n,k);
for i=1:n
    for L=i:n
        for j=1:k
        ga_t(i,j)= ga_t(i,j)+(-e(L,2)*m(L,L)*phi(L,j)+e(L,1)*m(n+L,n+L)*phi(n+L,j)+m(2*n+L,2*n+L)*phi(2*n+L,j))*W2(j,j);
        end
    end 
end


