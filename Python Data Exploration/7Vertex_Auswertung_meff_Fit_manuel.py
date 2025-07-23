# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 10:49:02 2023

@author: samue
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import linalg as sp
def linearFunc(x,intercept):
    y = intercept
    return y

modus = 2 #Modus 1: Variiere N0, N1=2. Modus 2: Variiere N1, N0=2**n_Gitter
n_Gitter = 5
n_masse = 40
M_min = 0.0
M_max = 4.0

u=1 #Anzahl Parameter Fit
m_eff_fit = np.zeros([n_Gitter,n_masse])
m_eff_fit_sd = np.zeros([n_Gitter,n_masse])

#%% Lese 2-Punkt Korrelation Daten und m_eff ein
N0 = 2**n_Gitter

propagator_long = np.zeros([n_Gitter,n_masse,2*N0-1])
propagator_long_sd =np.zeros([n_Gitter,n_masse,2*N0-1])
propagator = np.zeros([n_Gitter,n_masse,int(N0/2)+1])
propagator_sd = np.zeros([n_Gitter,n_masse,int(N0/2)+1])
m_eff = np.zeros([n_Gitter,n_masse,int(N0/2)])
m_eff_sd = np.zeros([n_Gitter,n_masse,int(N0/2)])
cov_long = np.zeros([n_Gitter,n_masse,2*N0-1,2*N0-1])
cov = np.zeros([n_Gitter,n_masse,int(N0/2)+1,int(N0/2)+1])
cov_m_eff = np.zeros([n_Gitter,n_masse,int(N0/2),int(N0/2)])

for k in range(n_Gitter):
    if modus==1:
        N0 = 2**n_Gitter
        N0dat = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N0dat = N0
        N1 = 2**(k+1)

    for i in range(n_masse):
        M = (i+1)*(M_max-M_min)/n_masse + M_min
        
        if modus==1:
            propagator_long[k,i,(N0-N0dat):(N0+N0dat-1)] = np.loadtxt("7Vertex_dt_number_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
            propagator_long_sd[k,i,(N0-N0dat):(N0+N0dat-1)] = np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
            cov_long[k,i,(N0-N0dat):(N0+N0dat-1),(N0-N0dat):(N0+N0dat-1)] = np.loadtxt("7Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        else:
            propagator_long[k,i] = np.loadtxt("7Vertex_dt_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
            propagator_long_sd[k,i] = np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
            cov_long[k,i] = np.loadtxt("7Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        
        c0 = propagator_long[k,i,N0-1]
        propagator_long[k,i] = propagator_long[k,i]/c0
        propagator_long_sd[k,i] = propagator_long_sd[k,i]/c0
        
        propagator[k,i,0] = propagator_long[k,i,N0-1]
        propagator[k,i,1:(int(N0dat/2)+1)] = \
                 (propagator_long[k,i,N0:(N0+int(N0dat/2))] \
                + np.flip(propagator_long[k,i,(N0+int(N0dat/2)-1):(N0+N0dat-1)])\
                + np.flip(propagator_long[k,i,(N0-int(N0dat/2)-1):(N0-1)])  \
                + propagator_long[k,i,(N0-N0dat):(N0-int(N0dat/2))])/2.0
        
        #A Matrix: Zum Berechnen von cov aus cov_long
        A = np.zeros([int(N0/2)+1,2*N0-1])
        A[0,N0-1] = 1.0/c0
        for j in range(int(N0/2)):
            A[j+1,N0-2-j] = 1.0/(2.0*c0)
            A[j+1,N0+j] = 1.0/(2.0*c0)
            A[j+1,j] = 1.0/(2.0*c0)
            A[j+1,2*N0-2-j] = 1.0/(2.0*c0)
        cov[k,i] = np.matmul(np.matmul(A,cov_long[k,i]),A.transpose())
        
        for l in range(int(N0/2)+1):
            propagator_sd[k,i,l] = np.sqrt(cov[k,i,l,l])
        
        m_eff[k,i,0:int(N0dat/2)] = -np.log(propagator[k,i,1:(int(N0dat/2)+1)] / propagator[k,i,0:int(N0dat/2)])
        #A Matrix: Zum Berechnen von cov_m_eff aus cov
        A = np.zeros([int(N0/2),int(N0/2)+1])
        for j in range(int(N0/2)):
            A[j,j] = 1.0 / propagator[k,i,j]
            A[j,j+1] = -1.0 / propagator[k,i,j+1]
        cov_m_eff[k,i] = np.matmul(np.matmul(A,cov[k,i]),A.transpose())
        
        for l in range(int(N0/2)):
            m_eff_sd[k,i,l] = np.sqrt(cov_m_eff[k,i,l,l])


# %% Effektive Masse Plot, fixe Masse und Gittergrösse (Modus 2)
N0=32
N1=2**2
M=2.0

i = int(np.log(N1)/np.log(2))-1
j = int((M-M_min)*n_masse/(M_max-M_min))-1
index = np.arange(int(N0/2))
plt.errorbar(index,m_eff[i,j],yerr=m_eff_sd[i,j],fmt=".-",label="M="+f'{M:.2f}')
plt.xlabel("t")
plt.ylabel("m_eff(t)")
plt.title("Effective Mass, N="+str(N0dat)+"x"+str(N1))
plt.legend(loc="upper right")
plt.show()

if len(np.where(m_eff[i,j] > 10**9)[0])==0:
    tmax= int(N0/2)-1
else:
    tmax = np.where(m_eff[i,j] > 10**9)[0][0]-1
tmax_tot = tmax
print("tmax =",tmax)

#%% Bestimme das AIC für verschiedene tmin
AIC = np.zeros(tmax-1)

for tmin in range(0,tmax-1):
    n=tmax-tmin

    a_fit,cov_fit=curve_fit(linearFunc,index[tmin:(tmax+1)],m_eff[i,j,tmin:(tmax+1)],\
                        sigma=cov_m_eff[i,j,tmin:(tmax+1),tmin:(tmax+1)],absolute_sigma=True)
    theorie = np.zeros(tmax-tmin+1)+a_fit[0]
    x = m_eff[i,j,tmin:(tmax+1)]-theorie
    cov = cov_m_eff[i,j,tmin:(tmax+1),tmin:(tmax+1)]
    logL = np.matmul( np.matmul(x.transpose(),sp.inv(cov)),x)/(n-u)
    AIC[tmin] = 2*u + n*np.log(logL)

index=np.arange(tmax-1)
plt.scatter(index,AIC)
plt.xlabel("tmin")
plt.ylabel("AIC")
plt.title("AIC für verschiedene tmin, tmax="+str(tmax)+", M="+f'{M:.2f}'+", N="+str(N0)+"x"+str(N1))
plt.show()

#%% Bestimme das AIC für verschiedene tmax (tmin manuell eingeben)
tmin = 1

AIC = np.zeros(tmax_tot-tmin-1)
for tmax in range(tmin+2,tmax_tot+1):
    n=tmax-tmin
    a_fit,cov_fit=curve_fit(linearFunc,index[tmin:(tmax+1)],m_eff[i,j,tmin:(tmax+1)],\
                        sigma=cov_m_eff[i,j,tmin:(tmax+1),tmin:(tmax+1)],absolute_sigma=True)
    theorie = np.zeros(tmax-tmin+1)+a_fit[0]
    x = m_eff[i,j,tmin:(tmax+1)]-theorie
    cov = cov_m_eff[i,j,tmin:(tmax+1),tmin:(tmax+1)]
    logL = np.matmul( np.matmul(x.transpose(),sp.inv(cov)),x)/(n-u)
    AIC[tmax-tmin-2] = 2*u + n*np.log(logL)

index=np.arange(tmax_tot-tmin-1)+tmin+2
plt.scatter(index,AIC)
plt.xlabel("tmax")
plt.ylabel("AIC")
plt.title("AIC für verschiedene tmax, tmin="+str(tmin)+", M="+f'{M:.2f}'+", N="+str(N0)+"x"+str(N1))
plt.show()

#%% Fit m_eff (tmax manuell eingeben)
tmax=9

n=tmax-tmin
a_fit,cov_fit=curve_fit(linearFunc,index[tmin:(tmax+1)],m_eff[i,j,tmin:(tmax+1)],\
                    sigma=cov_m_eff[i,j,tmin:(tmax+1),tmin:(tmax+1)],absolute_sigma=True)
m_eff_fit[i,j] = a_fit[0]
m_eff_fit_sd[i,j] = np.sqrt(cov_fit[0,0])
print("m_eff = ",a_fit[0])

#%% Falls man keinen Fit machen kann

m_eff_fit[i,j] = float('NaN')
m_eff_fit_sd[i,j] = 0.0

#%% Plotte m_eff_fit, sobald man alle Daten hat

for k in range(n_Gitter):
    N1=2**(k+1)
    index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
    plt.errorbar(index_Masse,m_eff_fit[k],yerr=m_eff_fit_sd[k],fmt=".-",label=("N1="+str(N1)))
plt.xlabel("M_bar")
plt.ylabel("m_eff")
plt.title("Effective Mass vs. Bar Mass, N="+str(N0)+"xN1")
plt.legend()
plt.show()