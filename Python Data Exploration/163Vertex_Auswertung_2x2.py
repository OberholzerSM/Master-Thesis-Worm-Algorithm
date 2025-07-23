# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 11:45:30 2022

@author: samue
"""

import numpy as np
import matplotlib.pyplot as plt
n_masse = 5

# %% Histogram Sign

theory = np.zeros((n_masse, 2))
number = np.zeros((n_masse, 2))
fehler = np.zeros((n_masse, 2))
Z = np.zeros(n_masse)
Z_est = np.zeros(n_masse)
Z_fehler = np.zeros(n_masse)

weights_sign = abs(np.loadtxt("163Vertex_weights_sign_N=2x2.dat"))
Z=weights_sign[:,0]-weights_sign[:,1]

for i in range(n_masse):
    M = float( (i+1)*2/n_masse )

    theory[i] = weights_sign[i]/sum(weights_sign[i])

    configsnumber = np.loadtxt(
        "163Vertex_sign_number_M="+f'{M:.2f}'+",N=2x2.dat")
    Z_est[i] = configsnumber[0]-configsnumber[1] 
    n_data = int(sum(configsnumber))
    number[i] = configsnumber / n_data

    sigma = np.loadtxt("163Vertex_sign_sd_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler[i] = sigma / n_data
    Z_fehler[i] = np.sqrt(sigma[0]**2 + sigma[1]**2)

index = (np.arange(n_masse)+1)*2.0/float(n_masse)
bar_width = 0.04

for i in range(2):
    if i == 0:
        sign = 1
    else:
        sign = -1
    plt.bar(index-0.5*bar_width, theory[:, i],
            bar_width, label="Theory", edgecolor="black")
    plt.bar(index+0.5*bar_width, number[:, i], bar_width,
            label="MCMC", color="Green", edgecolor="black")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Frequency")
    plt.title("Sign="+str(sign))
    plt.legend(loc=0)
    plt.show()

    plt.plot(index, np.zeros(n_masse), color="grey")
    plt.scatter(index, number[:, i]-theory[:, i],
                color='red', label="Deviation from Theory")
    plt.errorbar(index, number[:, i]-theory[:, i], yerr=fehler[:, i],
                 fmt="none", ecolor='orange', label="Statistical Error")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Deviation from Theory")
    plt.title("Sign="+str(sign))
    plt.show()

# %% Histogram Isospin

theory = np.zeros((n_masse, 5))
number = np.zeros((n_masse, 5))
fehler = np.zeros((n_masse, 5))

weights_isospin = np.loadtxt("163Vertex_weights_isospin_N=2x2.dat")

for i in range(n_masse):
    M = float( (i+1)*2/n_masse )

    theory[i] = weights_isospin[i]/sum(weights_isospin[i])

    configsnumber = np.loadtxt(
        "163Vertex_isospin_number_M="+f'{M:.2f}'+",N=2x2.dat")
    n_data = int(sum(configsnumber))
    number[i] = configsnumber / n_data

    sigma = np.loadtxt("163Vertex_Isospin_sd_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler[i] = sigma / n_data
fehler = abs(fehler)

index = (np.arange(n_masse)+1)*2.0/float(n_masse)
bar_width = 0.04

for i in range(5):
    plt.bar(index-0.5*bar_width, theory[:, i],
            bar_width, label="Theory", edgecolor="black")
    plt.bar(index+0.5*bar_width, number[:, i], bar_width,
            label="MCMC", color="Green", edgecolor="black")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Frequency")
    plt.title("Isospin="+str(i-2))
    plt.legend(loc=0)
    plt.show()

    plt.plot(index, np.zeros(n_masse), color="grey")
    plt.scatter(index, number[:, i]-theory[:, i],
                color='red', label="Deviation from Theory")
    plt.errorbar(index, number[:, i]-theory[:, i], yerr=fehler[:, i],
                 fmt="none", ecolor='orange', label="Statistical Error")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Deviation from Theory")
    plt.title("Isospin="+str(i-2))
    plt.show()

#%% 2-Punkt Korrelation
N0=2
N1=2
propagator = np.zeros([5,n_masse,2*N0-1])
theory = np.zeros([5,n_masse,2*N0-1])
fehler = np.zeros([5,n_masse,2*N0-1])
index = np.arange(2*N0-1) - 0.5*(2*N0-1-1)

theory_list = np.loadtxt("163Vertex_weights_dt_N=2x2.dat")
Z = np.loadtxt("163Vertex_Z_N=2x2.dat")[:,0]

for i in range(5):
    for j in range(n_masse):
        M = float( (j+1)*2/n_masse )
        theory[i,j] = theory_list[i+j*n_masse]/Z[j]
        propagator[i] = np.loadtxt("163Vertex_dt_number_Meson"+str(i+1)+"_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        fehler[i] = np.loadtxt("163Vertex_dt_sd_Meson"+str(i+1)+"_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        propagator[i,j] = propagator[i,j]/Z_est[j]
        fehler[i,j]= np.sqrt( (fehler[i,j]/Z_est[j])**2 + (propagator[i,j]*Z_fehler[j]/(Z_est[j]**2))**2 )
        
        plt.errorbar(index,propagator[i,j]-theory[i,j],yerr=fehler[i,j],fmt=".-",label="M="+f'{M:.2f}')
    plt.xlabel("t")
    plt.ylabel("c(t)")
    plt.title("Deviation from Theory, Wormtype "+str(i+1))
    plt.legend()
    plt.show()

#%% Transfermatrix Test
from numpy.linalg import matrix_power

Z_T = np.zeros(n_masse)
theory_T = np.zeros([n_masse,2*N0-1,5])
S_source = np.zeros([5,36,36])
S_sink = np.zeros([5,36,36])
S0 = np.zeros([5,36,36])

for i in range(n_masse):
    M = float( (i+1)*2/n_masse )
    T = np.loadtxt("163Vertex_T_M="+f'{M:.2f}'+".dat")
    S_source_tot = np.loadtxt("163Vertex_S_source_M="+f'{M:.2f}'+".dat")
    S_sink_tot = np.loadtxt("163Vertex_S_sink_M="+f'{M:.2f}'+".dat")
    S0_tot = np.loadtxt("163Vertex_S0_M="+f'{M:.2f}'+".dat")
    
    Z_T[i] = np.trace(matrix_power(T,N0))
    for j in range(5):
        S_source[j] = S_source_tot[(j*36):((j+1)*36)]
        S_sink[j] = S_sink_tot[(j*36):((j+1)*36)]
        S0[j] = S0_tot[(j*36):((j+1)*36)]
        
        theory_T[i,N0-1,j] = N0*np.trace(np.matmul(S0[j],matrix_power(T,N0-1)))/Z_T[i]
        for k in range(N0-1):
            matrix1 = np.matmul(S_source[j],matrix_power(T,k))
            matrix2 = np.matmul(matrix1,S_sink[j])
            matrix3 = np.matmul(matrix2,matrix_power(T,N0-2-k))
            theory_T[i,N0-2-k,j] = (N0-1-k)*np.trace(matrix3)/Z_T[i]
            theory_T[i,N0+k,j]=theory_T[i,N0-2-k,j]