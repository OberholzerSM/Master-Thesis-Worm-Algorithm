# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:03:19 2022

@author: samue
"""

import numpy as np
from numpy.linalg import matrix_power
from numpy.linalg import multi_dot
import matplotlib.pyplot as plt
n_masse = 20
M_min = 0.0
M_max = 2.0
k0=5
k1=1
N0 = 2**k0
N1 = 2**k1

theory_matrix = np.zeros([n_masse,2*N0-1])
theory2_matrix = np.zeros([n_masse,int(N0/2)+1])
number_matrix = np.zeros([n_masse,2*N0-1])
error_matrix = np.zeros([n_masse,2*N0-1])
index=np.arange(2*N0-1)-N0+1

for i in range(n_masse):
    M = (i+1)*(M_max-M_min)/n_masse + M_min
    
    T = np.array([[M**4+1,0,0,0.5],[0,M**2,0.5,0],[0,0.5,M**2,0],[0.5,0,0,1.0]])
    S = np.array([[0,M**2+1,M**2+1,0],[M**2+1,0,0,1],[M**2+1,0,0,1],[0,1,1,0]])
    S0 = np.array([[4*M**2+8,0,0,2],[0,4,2,0],[0,2,4,0],[2,0,0,0]])
    
    theory_matrix[i,N0-1] = N0*np.trace(np.dot(S0,matrix_power(T,N0-1)))
    for j in range(N0-1):
        theory_matrix[i,N0+j] = (N0-1-j)*np.trace( multi_dot([S,matrix_power(T,j),S,matrix_power(T,N0-2-j)]) )
        theory_matrix[i,N0-2-j] = theory_matrix[i,N0+j]
    theory_matrix[i] = theory_matrix[i]/sum(theory_matrix[i])
    
    theory2_matrix[i,0] = theory_matrix[i,N0-1]
    theory2_matrix[i,1:(int(N0/2)+1)] =   (theory_matrix[i,0:int(N0/2)] +\
                            np.flip(theory_matrix[i,(int(N0/2)-1):(N0-1)])+\
                            theory_matrix[i,N0:int(3*N0/2)] +\
                            np.flip(theory_matrix[i,(int(3*N0/2)-1):(2*N0-1)]))/2.0
    
    number_matrix[i] = np.loadtxt("7Vertex_dt_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    error_matrix[i] = np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
    n_data=sum(number_matrix[i])
    number_matrix[i] = number_matrix[i]/n_data
    #error_matrix[i] = error_matrix[i]/n_data

#%% Plotte die Theoretische Korrelationsfunktion

i=7
M = (i+1)*(M_max-M_min)/n_masse + M_min
plt.plot(index,theory_matrix[i])
plt.scatter(index,theory_matrix[i],label="M="+f'{M:.2f}')
plt.xlabel("t")
plt.ylabel("c(t)")
plt.title("2-point Correlation Function, M="+f'{M:.2f}'+", N="+str(N0)+"x"+str(N1))
plt.show()

index2 = np.arange(int(N0/2)+1)
plt.plot(index2,theory2_matrix[i])
plt.scatter(index2,theory2_matrix[i],label="M="+f'{M:.2f}')
plt.xlabel("t")
plt.ylabel("cS(t)")
plt.title("Symmetrized 2-point Correlation Function, M="+f'{M:.2f}'+", N="+str(N0)+"x"+str(N1))
plt.show()

#%% Plot Abweichung Theorie
#a=N0+10
#b=N0+20
a=0
b=2*N0-1
for i in range(n_masse):
    M = (i+1)*(M_max-M_min)/n_masse + M_min
    plt.errorbar(index[a:b],number_matrix[i,a:b]-theory_matrix[i,a:b],yerr=error_matrix[i,a:b],fmt=".-",label="M="+f'{M:.2f}')
    plt.plot(index[a:b],np.zeros(2*N0-1)[a:b],color="grey")
    plt.xlabel("t")
    plt.ylabel("c(t)/c(0)")
    plt.title("Deviation from Theory 2-point Correlation Function, N="+str(N0)+"x"+str(N1))
    plt.legend()
    plt.show()