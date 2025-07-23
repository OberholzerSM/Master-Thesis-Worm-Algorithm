# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 08:59:08 2023

@author: samue
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import matrix_power

#%% Z Abweichung Theorie für verschiedene Gittergrössen
n_Gitter = 5
M_min = 0.0
M_max = 4.0
n_masse = 40
N1=2
Z = np.zeros([n_Gitter,n_masse])
Z_abs = np.zeros([n_Gitter,n_masse])
Z_plus = np.zeros([n_Gitter,n_masse])
Z_minus = np.zeros([n_Gitter,n_masse])
n_plus = np.zeros((n_Gitter,n_masse))
n_minus = np.zeros((n_Gitter,n_masse))
n_plus_sd = np.zeros((n_Gitter,n_masse))
n_minus_sd = np.zeros((n_Gitter,n_masse))

for i in range(n_Gitter):
    N0 = 2**(i+1)
    for j in range(n_masse):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        
        T = np.loadtxt("163Vertex_T_M="+f'{M:.2f}'+".dat")
        Z[i,j] = np.trace(matrix_power(T,N0))
        
        T_abs = np.loadtxt("163Vertex_T_abs_M="+f'{M:.2f}'+".dat")
        Z_abs[i,j] = np.trace(matrix_power(T_abs,N0))
        
        Z_plus[i,j] = (Z[i,j]+Z_abs[i,j])/(2.0*Z_abs[i,j])
        Z_minus[i,j] = (-Z[i,j]+Z_abs[i,j])/(2.0*Z_abs[i,j])
        
        sign_number = np.loadtxt(
            "163Vertex_sign_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        n_data = int(sum(sign_number))
        n_plus[i,j] = np.copy(sign_number[0])/n_data
        n_minus[i,j] = np.copy(sign_number[1])/n_data
        
        sign_sd = np.loadtxt(
            "163Vertex_sign_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        n_plus_sd[i,j] = np.copy(sign_sd[0])/n_data
        n_minus_sd[i,j] = np.copy(sign_sd[1])/n_data
              
index = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
for i in range(n_Gitter):
    N0 = 2**(i+1)
    plt.errorbar(index,n_plus[i]-Z_plus[i],yerr=n_plus_sd[i],fmt=".-")
    plt.plot(index,np.zeros(n_masse),color="grey")
    plt.xlabel("M")
    plt.ylabel("Deviation from Theory")
    plt.title("Ratio Z and Z_abs, Deviation from Theory, N="+str(N0)+"x2")
    plt.show()

#%% Propagator Abweichung Theorie für Fixe Gittergrösse
n_masse = 10
M_min = 2.0
M_max = 4.0
k0=5
k1=1
N0 = 2**k0
N1 = 2**k1
tmin = N0-1
tmax = 2*N0-1

Z = np.zeros(n_masse)
theory = np.zeros([5,n_masse,2*N0-1])
S_source = np.zeros([5,36,36])
S_sink = np.zeros([5,36,36])
S0 = np.zeros([5,36,36])

# Transfermatrix für fixe Gittergrösse
for j in range(n_masse):
    M = (j+1)*(M_max-M_min)/n_masse + M_min
    T = np.loadtxt("163Vertex_T_M="+f'{M:.2f}'+".dat")
    S_source_tot = np.loadtxt("163Vertex_S_source_M="+f'{M:.2f}'+".dat")
    S_sink_tot = np.loadtxt("163Vertex_S_sink_M="+f'{M:.2f}'+".dat")
    S0_tot = np.loadtxt("163Vertex_S0_M="+f'{M:.2f}'+".dat")
    
    Z[j] = np.trace(matrix_power(T,N0))
    for i in range(5):
        S_source[i] = S_source_tot[(i*36):((i+1)*36)]
        S_sink[i] = S_sink_tot[(i*36):((i+1)*36)]
        S0[i] = S0_tot[(i*36):((i+1)*36)]
        
        theory[i,j,N0-1] = N0*np.trace(np.matmul(S0[i],matrix_power(T,N0-1)))/Z[j]
        for k in range(N0-1):
            matrix1 = np.matmul(S_source[i],matrix_power(T,k))
            matrix2 = np.matmul(matrix1,S_sink[i])
            matrix3 = np.matmul(matrix2,matrix_power(T,N0-2-k))
            theory[i,j,N0-2-k] = (N0-1-k)*np.trace(matrix3)/Z[j]
            theory[i,j,N0+k]=theory[i,j,N0-2-k]
    

# 2-Punkt Korrelation
propagator = np.zeros([5,n_masse,2*N0-1])
fehler = np.zeros([5,n_masse,2*N0-1])
index = np.arange(2*N0-1) - 0.5*(2*N0-1-1)
Z_est = np.zeros(n_masse)
Z_fehler = np.zeros(n_masse)

for j in range(n_masse):
    
    M = (j+1)*(M_max-M_min)/n_masse + M_min
    
    sign = np.loadtxt("163Vertex_sign_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    Z_est[j] = sign[0]-sign[1]
    
    sigma = np.loadtxt("163Vertex_sign_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    Z_fehler[j] = np.sqrt(sigma[0]**2 + sigma[1]**2)
    
    cov = np.loadtxt("163Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    for i in range(5):
        propagator[i,j] = np.loadtxt("163Vertex_dt_number_Meson"+str(i+1)+"_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        for t in range(2*N0-1):
            fehler[i,j,t] = np.sqrt( cov[t+i*(2*N0-1),t+i*(2*N0-1)] )*(10**5)
        propagator[i,j]=propagator[i,j]/Z_est[j] 
        fehler[i,j]= np.sqrt( (fehler[i,j]/Z_est[j])**2 + (propagator[i,j]*Z_fehler[j]/(Z_est[j]**2))**2 )

for i in range(5):
    for j in range(n_masse):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[tmin:tmax],propagator[i,j,tmin:tmax]-theory[i,j,tmin:tmax],\
                     yerr=fehler[i,j,tmin:tmax],fmt=".-",label="M="+f'{M:.2f}')
        plt.plot(index[tmin:tmax],np.zeros(2*N0-1)[tmin:tmax],color="grey")
    plt.xlabel("t")
    plt.ylabel("Deviation from Theory")
    plt.title("2-point correlation function, Deviation from Theory, Meson "+str(i+1)+", N="+str(N0)+"x"+str(N1))
    plt.legend(loc="upper right")
    plt.show()

#%% Propagator abs Abweichung Theorie für Fixe Gittergrösse
Z = np.zeros(n_masse)
theory = np.zeros([5,n_masse,2*N0-1])
S_source = np.zeros([5,36,36])
S_sink = np.zeros([5,36,36])
S0 = np.zeros([5,36,36])

for j in range(n_masse):
    M = (j+1)*(M_max-M_min)/n_masse + M_min
    T = np.loadtxt("163Vertex_T_abs_M="+f'{M:.2f}'+".dat")
    S_source_tot = np.loadtxt("163Vertex_S_source_abs_M="+f'{M:.2f}'+".dat")
    S_sink_tot = np.loadtxt("163Vertex_S_sink_abs_M="+f'{M:.2f}'+".dat")
    S0_tot = np.loadtxt("163Vertex_S0_abs_M="+f'{M:.2f}'+".dat")
    
    Z[j] = np.trace(matrix_power(T,N0))
    for i in range(5):
        S_source[i] = S_source_tot[(i*36):((i+1)*36)]
        S_sink[i] = S_sink_tot[(i*36):((i+1)*36)]
        S0[i] = S0_tot[(i*36):((i+1)*36)]
        
        theory[i,j,N0-1] = N0*np.trace(np.matmul(S0[i],matrix_power(T,N0-1)))/Z[j]
        for k in range(N0-1):
            matrix1 = np.matmul(S_source[i],matrix_power(T,k))
            matrix2 = np.matmul(matrix1,S_sink[i])
            matrix3 = np.matmul(matrix2,matrix_power(T,N0-2-k))
            theory[i,j,N0-2-k] = (N0-1-k)*np.trace(matrix3)/Z[j]
            theory[i,j,N0+k]=theory[i,j,N0-2-k]
    

# 2-Punkt Korrelation abs
propagator = np.zeros([5,n_masse,2*N0-1])
fehler = np.zeros([5,n_masse,2*N0-1])
index = np.arange(2*N0-1) - 0.5*(2*N0-1-1)
Z_est = np.zeros(n_masse)
Z_fehler = np.zeros(n_masse)

for j in range(n_masse):
    
    M = (j+1)*(M_max-M_min)/n_masse + M_min
    
    sign = np.loadtxt("163Vertex_sign_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    Z_est[j] = sign[0]+sign[1]
    
    sigma = np.loadtxt("163Vertex_sign_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    Z_fehler[j] = np.sqrt(sigma[0]**2 + sigma[1]**2)
    
    cov = np.loadtxt("163Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    for i in range(5):
        propagator[i,j] = np.loadtxt("163Vertex_dt_number_abs_Meson"+str(i+1)+"_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        for t in range(2*N0-1):
            fehler[i,j,t] = np.sqrt( cov[t+i*(2*N0-1),t+i*(2*N0-1)] )*(10**5)
        propagator[i,j]=propagator[i,j]/Z_est[j] 
        fehler[i,j]= np.sqrt( (fehler[i,j]/Z_est[j])**2 + (propagator[i,j]*Z_fehler[j]/(Z_est[j]**2))**2 )

for i in range(5):
    for j in range(n_masse):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[tmin:tmax],propagator[i,j,tmin:tmax]-theory[i,j,tmin:tmax],\
                     yerr=fehler[i,j,tmin:tmax],fmt=".-",label="M="+f'{M:.2f}')
        plt.plot(index[tmin:tmax],np.zeros(2*N0-1)[tmin:tmax],color="grey")
    plt.xlabel("t")
    plt.ylabel("Deviation from Theory")
    plt.title("2-point correlation function abs, Deviation from Theory, Meson "+str(i+1)+", N="+str(N0)+"x"+str(N1))
    plt.legend(loc="upper right")
    plt.show()

#%% Z Theorie für verschiedene Gittergrössen

n_Gitter = 128
M_min = 0.0
M_max = 2.0
n_masse = 10
N1=2
Z = np.zeros([n_Gitter,n_masse])
Z_abs = np.zeros([n_Gitter,n_masse])
Z_plus = np.zeros([n_Gitter,n_masse])
Z_minus = np.zeros([n_Gitter,n_masse])

for i in range(n_Gitter):
    N0 = i+1
    for j in range(n_masse):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        
        T = np.loadtxt("163Vertex_T_M="+f'{M:.2f}'+".dat")
        Z[i,j] = np.trace(matrix_power(T,N0))
        
        T_abs = np.loadtxt("163Vertex_T_abs_M="+f'{M:.2f}'+".dat")
        Z_abs[i,j] = np.trace(matrix_power(T_abs,N0))
        
        Z_plus[i,j] = (Z[i,j]+Z_abs[i,j])/(2.0*Z_abs[i,j])
        Z_minus[i,j] = (-Z[i,j]+Z_abs[i,j])/(2.0*Z_abs[i,j])
              
index = (1+np.arange(n_Gitter))
for j in range(n_masse):
    M = (j+1)*(M_max-M_min)/n_masse + M_min
    plt.plot(index,Z[:,j]/Z_abs[:,j],label="M="+f'{M:.2f}')
    plt.plot(index,np.zeros(n_Gitter),color="grey")
plt.xlabel("N0")
plt.ylabel("Z/Z_abs")
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.tight_layout()
plt.title("Ratio of Z and Z_abs (Transfermatrices), N=N0x2")
plt.show()

#%% Symmetrisierter Propagator Theorie
n_masse = 10
M_min = 2.0
M_max = 4.0
k0=5
k1=1
N0 = 2**k0
N1 = 2**k1

Z = np.zeros(n_masse)
theory_long = np.zeros([5,n_masse,2*N0-1])
theory = np.zeros([5,n_masse,int(N0/2)+1])
S_source = np.zeros([5,36,36])
S_sink = np.zeros([5,36,36])
S0 = np.zeros([5,36,36])

# Transfermatrix für fixe Gittergrösse
for j in range(n_masse):
    M = (j+1)*(M_max-M_min)/n_masse + M_min
    T = np.loadtxt("163Vertex_T_M="+f'{M:.2f}'+".dat")
    S_source_tot = np.loadtxt("163Vertex_S_source_M="+f'{M:.2f}'+".dat")
    S_sink_tot = np.loadtxt("163Vertex_S_sink_M="+f'{M:.2f}'+".dat")
    S0_tot = np.loadtxt("163Vertex_S0_M="+f'{M:.2f}'+".dat")
    Z[j] = np.trace(matrix_power(T,N0))
    for i in range(5):
        S_source[i] = S_source_tot[(i*36):((i+1)*36)]
        S_sink[i] = S_sink_tot[(i*36):((i+1)*36)]
        S0[i] = S0_tot[(i*36):((i+1)*36)]
        
        theory_long[i,j,N0-1] = N0*np.trace(np.matmul(S0[i],matrix_power(T,N0-1)))/Z[j]
        for k in range(N0-1):
            matrix1 = np.matmul(S_source[i],matrix_power(T,k))
            matrix2 = np.matmul(matrix1,S_sink[i])
            matrix3 = np.matmul(matrix2,matrix_power(T,N0-2-k))
            theory_long[i,j,N0-2-k] = (N0-1-k)*np.trace(matrix3)/Z[j]
            theory_long[i,j,N0+k]=theory_long[i,j,N0-2-k]
        
        theory[i,j,0] = theory_long[i,j,N0-1]
        theory[i,j,1:(int(N0/2)+1)]= \
                 (theory_long[i,j,N0:(N0+int(N0/2))] \
                + np.flip(theory_long[i,j,(N0+int(N0/2)-1):(2*N0-1)])\
                + np.flip(theory_long[i,j,(N0-int(N0/2)-1):(N0-1)])  \
                + theory_long[i,j,0:int(N0/2)])/2.0
        
        theory[i,j] = theory[i,j]/theory[i,j,0]

index = np.arange(int(N0/2)+1)
for i in range(5):
    for j in range(n_masse):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        plt.plot(np.log(theory[i,j]))
        plt.scatter(index,np.log(theory[i,j]),label="M="+f'{M:.2f}')
    plt.xlabel("t")
    plt.ylabel("ln[cS(t)/cS(0)]")
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.title("Logarithmic Symmetrized Correlation Function, Meson "+str(i+1)+", N="+str(N0)+"x"+str(N1))
    plt.show()
#%% Effektive Masse
t_min=1
index = np.arange(int(N0/2)-t_min)+t_min
names=["Meson 1/2","Meson 3/4","Meson 5"]
for i in range(3):
    for j in range(n_masse):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        if i==0:
            plt.plot(index,-np.log(theory[i,j,(t_min+1):(int(N0/2)+1)]/theory[i,j,t_min:int(N0/2)]))
            plt.scatter(index,-np.log(theory[i,j,(t_min+1):(int(N0/2)+1)]/theory[i,j,t_min:int(N0/2)]),label="M="+f'{M:.2f}')
        elif i==1:
            plt.plot(index,-np.log(theory[2,j,(t_min+1):(int(N0/2)+1)]/theory[2,j,t_min:int(N0/2)]))
            plt.scatter(index,-np.log(theory[2,j,(t_min+1):(int(N0/2)+1)]/theory[2,j,t_min:int(N0/2)]),label="M="+f'{M:.2f}')
        elif i==2:
            plt.plot(index,-np.log(theory[4,j,(t_min+1):(int(N0/2)+1)]/theory[4,j,t_min:int(N0/2)]))
            plt.scatter(index,-np.log(theory[4,j,(t_min+1):(int(N0/2)+1)]/theory[4,j,t_min:int(N0/2)]),label="M="+f'{M:.2f}')
    plt.xlabel("t")
    plt.ylabel("m_eff(t)")
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.title("Effective Mass, "+names[i]+", N="+str(N0)+"x"+str(N1))
    plt.show()