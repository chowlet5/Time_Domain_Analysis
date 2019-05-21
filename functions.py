import pandas as pd
import numpy as np
import math
from scipy.special import erfinv, erf
from scipy.signal import welch

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
    B = np.zeros((1,1))
    for i in range(n):
        
        B[i,i] = A[i]
    return B
def cdf_std_norm(x):
    return 0.5*(1+erf(x/math.sqrt(2)))

def q_explicit (M,P,C,K,k,N,dt):
    q = np.zeros((1,1))
    qd = np.zeros((1,1))
    qdd = np.zeros((1,1))
    for i in range(k):
        for j in range(N):
            if j ==1:
                q[i,j] = 0.0
                qd[i,j] = 0.0
                qdd[i,j] = (1/M([i,i]))*(P[i,j]-C[i,i]*qd[i,j]-K[i,i]*q[i,j])

def q_implicit (M,P,C,K,k,N,dt):

    ga = 0.5
    be = 0.166667

    dpp=0.0
    dq=0.0
    dqd=0.0
    dqdd=0.0

    q = np.zeros((k,N))
    qd = np.zeros((k,N))
    qdd = np.zeros((k,N))
    KP = np.zeros((k,k))
    a = np.zeros((k,k))
    b = np.zeros((k,k))

    for i in range(k):
        for j in range(k):
            KP[i,j] = K[i,j]+(ga/(be*dt))*(C[i,j]+(1/(be*(dt*dt)))*M[i,j])
            a[i,j] = (1/(be*dt))*M[i,j]+(ga/be)*C[i,j]
            b[i,j] = (1/(2*be))*M[i,j]+dt*((ga/(2*be))-1)*C[i,j]

    for i in range(k):
        for j in range(N):
            if j == 1:
                dpp = P[i,j]+a[i,i]*qd[i,j]+b[i,i]*qdd[i,j]

                dq = (1/KP[i,i])*dpp
                dqd = (ga/(be*dt))*dq-(ga/be)*qd[i,j]+dt*(1-(ga/(2*be)))*qdd[i,j]
                dqdd =(1/(be*dt*dt))*dq-(1/(be*dt))*qd[i,j]-(1/(2*be))*qdd[i,j]


                q[i,j]=q[i,j]+dq
                qd[i,j]=qd[i,j]+dqd
                qdd[i,j]=qdd[i,j]+dqdd
            else:
                q[i,j]=q[i,j-1]
                qd[i,j]=qd[i,j-1]
                qdd[i,j]=qdd[i,j-1]

                dpp=P[i,j]-P[i,j-1]+a[i,i]*qd[i,j]+b[i,i]*qdd[i,j]
                dq = (1/KP[i,i])*dpp
                dqd = (ga/(be*dt))*dq-(ga/be)*qd[i,j]+dt*(1-(ga/(2*be)))*qdd[i,j]
                dqdd =(1/(be*dt*dt))*dq-(1/(be*dt))*qd[i,j]-(1/(2*be))*qdd[i,j]
                
                q[i,j]=q[i,j-1]+dq
                qd[i,j]=qd[i,j-1]+dqd
                qdd[i,j]=qdd[i,j-1]+dqdd

def peak (data,dur_ratio):

    if not data:
        raise Exception('Time series input expected')
    
    num_CDF = 1000
    min_CDF = 0.0005
    max_CDF = 0.9995
    CDF_pk = np.linspace(min_CDF,max_CDF, num_CDF)
    
    sdata = data.shape
    max_est = np.zeros((sdata[0],1))
    min_est = np.zeros((sdata[0],1))

    for i in range(sdata[0]):

        X = data[i,:]
        n = len(X)
        avg_X = np.mean(X)
        sorted_X = sorted(X)

        CDF_X = (np.arange(1,n+1))/(n+1)
        snv = isnc(CDF_X)
        avg_snv = np.mean(snv)
        
        sigma = (np.sum(np.multiply(snv,sorted_X)) - n*avg_snv*avg_X)/(np.sum(np.power(snv,2))-n*avg_snv**2)
        mu = avg_X - sigma*avg_snv
        X_fit = mu + sigma*snv


        norm_PPCC = sigma *np.std(snv)/np.std(sorted_X)

        '''
        --------------------------------------------------------------------------
        Estimate the mean zero upcorssing rate of a process y(t) with standard
        normal probability distribution using the classical Rice(1954) results as
        follow:
        --------------------------------------------------------------------------
        Estimate the interval of integration in frequency domain. This is done
        by matching variances obtained from the timehistory data and integration
        of the spectra in the frequency domain
        Variance from time history
        '''
        stdX = np.std(X)
        varX = np.power(stdX,2)

        df = 65536
        fs = 5.12

        f, S_X = welch(X,fs,nfft=df, detrend=False)
        sf = (df*0.5)
        si = 0
        var_X = (fs/df)*np.sum(S_X[si:sf])
        temp = var_X
        si = 1
        var_X  = (fs/df)*np.sum(S_X[si:sf])

        while (abs(var_X-varX)<abs(temp-varX)) and (var_X>varX):
            temp = var_X
            si = si+1
            var_X = (fs/df)*np.sum(S_X[si:sf])

        var_X = temp

        si = si -1 
        
        numer = np.trapz(f[si:sf],np.multiply(np.power(f[si:sf],2),S_X[si:sf]))
        denom = np.trapz(f[si:sf],np.multiply(f[si:sf],S_X[si:sf]))

        nu_y = np.sqrt(numer/denom)

        y_pk = np.sqrt(2*np.log(np.divide(-dur_ratio*nu_y*3600,np.log(CDF_pk))))

        X_max = y_pk *sigma +mu
        X_min = -y_pk *sigma +mu

        pdf_pk = np.multiply(np.multiply(-y_pk,CDF_pk),np.log(CDF_pk))
        max_est[i] = np.trapz(y_pk,np.multiply(pdf_pk,X_max))
        min_est[i] = np.trapz(y_pk,np.multiply(pdf_pk,X_min))
    
    return max_est , min_est
def td_response(gama,q,n,k,N):

    mr = np.zeros((n,k,N))

    for i in range(n):
        for j in range(k):
            for l in range(N):
                mr[i,j,l] = gama[i,j]*q[j,l]
    
    tr = np.zeros((n,N))

    for i in range(n):
        for j in range(k):
            for l in range(N):
                tr[i,l] = tr[i,l] + gama[i,j]*q[j,l]
    
    return mr, tr

