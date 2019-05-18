import pandas as pd
import numpy as np
import math
from scipy.special import erfinv, erf

def gamma (phi, W2,m,h,e,n,k):
    #modal participation coefficient for displacement
    ga_d = phi

    #modal participation coefficient for acceleration

    ga_a = np.multiply(phi,W2)

    #modal participation coefficient for storey shear

    ga_v = np.zeros((2*n,k))

    for i in range(n):
        for l in range(n):
            for j in range(k):
                ga_v[i,j] = ga_v[i,j] + m[l,l]*phi[l,j]*W2[j,j]
                ga_v[n+i,j] = ga_v[n+i,j] + m[n+l,n+l]*phi[n+l,j]*W2[j,j]
    
    #modal participation coefficient for bending moment

    ga_m = np.zeros((2*n,k))

    for i in range(n):
        for l in range(n):
            for j in range(k):
                if i==1:
                    ga_m[i,j]= ga_m[i,j]+ h[l]*m[l,l]*phi[l,j]*W2[j,j]
                    ga_m[n+i,j]=ga_m[n+i,j]+ h[l]*m[n+l,n+l]*phi[n+l,j]*W2[j,j]
                else:
                    ga_m[i,j]= ga_m[i,j]+ (h[l]-h[i-1])*m[l,l]*phi[l,j]*W2[j,j]
                    ga_m[n+i,j]=ga_m[n+i,j]+ (h[l]-h[i-1])*m[n+l,n+l]*phi[n+l,j]*W2[j,j]

    #modal participation coefficient for torsion

    ga_t = np.zeros((n,k))

    for i in range(n):
        for l in range(n):
            for j in range(k):
                ga_t[i,j] = ga_t[i,j] + (-e[l,2]*m[l,l]*phi[l,j]+e[l,1]*m[n+l,n+l]*phi[n+l,j]+m[2*n+l,2*n+l]*phi(2*n+l,j))*W2(j,j)
    

    return ga_d, ga_a, ga_v, ga_m, ga_t

def isnc(p):
    
    return math.sqrt(2)*erfinv(p)

def mat_diag(A,n):
    for i in range(n):
        pass
        #B[i,i] = 

def cdf_std_norm(x):
    return 0.5*(1+erf(x/math.sqrt(2)))