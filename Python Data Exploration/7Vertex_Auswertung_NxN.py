# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 10:40:17 2022

@author: samue
"""

import numpy as np
import matplotlib.pyplot as plt
modus = 2 #Modus 1: Variiere N0, N1=2. Modus 2: Variiere N1, N0=2**n_Gitter
n_Gitter = 5
n_masse = 40
M_min = 0.0
M_max = 4.0

#%% Topologische Klassen Häufigkeit
top_number = np.zeros([n_Gitter,n_masse,4])
top_sd = np.zeros([n_Gitter,n_masse,4])
index = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
names=["Z_00","Z_10","Z_01","Z_11"]

for k in range(n_Gitter):
    if modus==1:
        N0 = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N1 = 2**(k+1)

    for i in range(n_masse):
        M = (i+1)*(M_max-M_min)/n_masse + M_min
        
        top_number[k,i] = np.loadtxt("7Vertex_topology_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        top_sd[k,i] = np.loadtxt("7Vertex_topology_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        n_data = int(sum(top_number[k,i]))
        top_number[k,i] = top_number[k,i] / n_data
        top_sd[k,i] = top_sd[k,i] / n_data

    for i in range(4):
        plt.errorbar(index,top_number[k,:,i],yerr=top_sd[k,:,i],fmt=".-",label=names[i])
    plt.xlabel("M_bar")
    plt.ylabel("Frequency")
    plt.legend(loc="upper left")
    plt.ylim([-0.1,1.1])
    plt.title("Frequency Topological Classes, N="+str(N0)+"x"+str(N1))
    plt.show()


#Kritischer Punkt
kritisch = np.zeros([n_Gitter,n_masse])
kritisch_sd = np.zeros([n_Gitter,n_masse])
for k in range(n_Gitter):
    if modus==1:
        N0 = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**5
        N1 = 2**(k+1)
    for i in range(n_masse):
        kritisch[k,i] = top_number[k,i,0]-top_number[k,i,1]-top_number[k,i,2]-top_number[k,i,3]
        kritisch_sd[k,i] = np.sqrt( top_sd[k,i,0]**2 + top_sd[k,i,1]**2 + top_sd[k,i,2]**2 + top_sd[k,i,3]**2 )
    if modus==1:
        plt.errorbar(index,kritisch[k],yerr=kritisch_sd[k],fmt=".-",label="N0="+str(N0))
    else:
        plt.errorbar(index,kritisch[k],yerr=kritisch_sd[k],fmt=".-",label="N1="+str(N1))
plt.plot(index,np.zeros(n_masse),color="grey")
plt.xlabel("M_bar")
plt.ylabel("Z_pp/Z")
plt.legend(loc="lower right")
plt.ylim([-1.1,1.1])
if modus==1:
    plt.title("Critical Point, N=N0x"+str(N1))
else:
    plt.title("Critical Point, N="+str(N0)+"xN1")
plt.show()

#%% Lese 2-Punkt Korrelation Daten und m_eff ein
N0 = 2**n_Gitter
n_data = 10**5

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
            propagator_long_sd[k,i,(N0-N0dat):(N0+N0dat-1)] = n_data*np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
            cov_long[k,i,(N0-N0dat):(N0+N0dat-1),(N0-N0dat):(N0+N0dat-1)] = n_data*np.loadtxt("7Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        else:
            propagator_long[k,i] = np.loadtxt("7Vertex_dt_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
            propagator_long_sd[k,i] = n_data*np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
            cov_long[k,i] = n_data*np.loadtxt("7Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        
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

# %% Masse Range für Plots
M_min_Plot = 0.0
M_max_Plot = 2.0
n_masse_Plot = 10

dM = (M_max_Plot-M_min_Plot)/n_masse_Plot
i_min = int( (M_min_Plot+dM-M_min)*n_masse/(M_max-M_min) )-1
i_max = int( (M_max_Plot-M_min)*n_masse/(M_max-M_min) )
di = int((i_max-i_min+1)/n_masse_Plot)

# %% Plot 2-Punkt Korrelation gesamt
for k in range(n_Gitter):
    if modus==1:
        N0 = 2**n_Gitter
        N0dat = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N0dat = N0
        N1 = 2**(k+1)

    index = np.arange(2*N0-1)-N0+1
    for i in range(i_min,i_max,di):
        M = (i+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[(N0-N0dat):(N0+N0dat-1)],propagator_long[k,i,(N0-N0dat):(N0+N0dat-1)],\
                     yerr=propagator_long_sd[k,i,(N0-N0dat):(N0+N0dat-1)],fmt=".-",label="M="+f'{M:.2f}')
    plt.xlabel("t")
    plt.ylabel("c(t)/c(0)")
    plt.title("2-point Correlation Function, N="+str(N0dat)+"x"+str(N1))
    plt.ylim([-0.05,1.05])
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.show()
# %% Logarithmic Plot 2-Punkt Korrelation gesamt
for k in range(n_Gitter):
    if modus==1:
        N0 = 2**n_Gitter
        N0dat = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N0dat = N0
        N1 = 2**(k+1)
    index = np.arange(2*N0-1)-N0+1
    for i in range(i_min,i_max,di):
        M = (i+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[(N0-N0dat):(N0+N0dat-1)],np.log(propagator_long[k,i,(N0-N0dat):(N0+N0dat-1)]),
                     yerr=propagator_long_sd[k,i,(N0-N0dat):(N0+N0dat-1)]/propagator_long[k,i,(N0-N0dat):(N0+N0dat-1)],fmt=".-",label="M="+f'{M:.2f}')
    plt.xlabel("t")
    plt.ylabel("ln[c(t)/c(0)]")
    plt.title("Logarithmic 2-point Correlation Function, N="+str(N0dat)+"x"+str(N1))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.ylim([-19.0,0.5])
    plt.show()
           
# %% Plot 2-Punkt Korrelation zusammengefasst
for k in range(n_Gitter):
    if modus==1:
        N0 = 2**n_Gitter
        N0dat = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N0dat = N0
        N1 = 2**(k+1)

    index = np.arange(int(N0/2)+1)
    for i in range(i_min,i_max,di):
        M = (i+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[0:(int(N0dat/2)+1)],propagator[k,i,0:(int(N0dat/2)+1)],\
                     yerr=propagator_sd[k,i,0:(int(N0dat/2)+1)],fmt=".-",label="M="+f'{M:.2f}')
    plt.xlabel("t")
    plt.ylabel("cS(t)/cS(0)")
    plt.title("Symmetrized 2-point Correlation Function, N="+str(N0dat)+"x"+str(N1))
    plt.ylim([-0.05,1.05])
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.show()
# %% Logarithmic Plot 2-Punkt Korrelation zusammengefasst
for k in range(n_Gitter):
    if modus==1:
        N0 = 2**n_Gitter
        N0dat = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N0dat = N0
        N1 = 2**(k+1)
    index = np.arange(int(N0/2)+1) 
    for i in range(i_min,i_max,di):
        M = (i+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[0:(int(N0dat/2)+1)],np.log(propagator[k,i,0:(int(N0dat/2)+1)]),
                     yerr=propagator_sd[k,i,0:(int(N0dat/2)+1)]/propagator[k,i,0:(int(N0dat/2)+1)],fmt=".-",label="M="+f'{M:.2f}')
    plt.xlabel("t")
    plt.ylabel("ln[cS(t)/cS(0)]")
    plt.title("Symmetrized Logarithmic 2-point Correlation Function, N="+str(N0dat)+"x"+str(N1))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.ylim([-19.0,0.5])
    plt.show()

# %% Effektive Masse Plot
for k in range(n_Gitter):
    if modus==1:
        N0 = 2**n_Gitter
        N0dat = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N0dat = N0
        N1 = 2**(k+1)
    index = np.arange(int(N0/2))

    for i in range(i_min,i_max,di):
        M = (i+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[0:(int(N0dat/2))],m_eff[k,i,0:(int(N0dat/2))],\
                     yerr=m_eff_sd[k,i,0:(int(N0dat/2))],fmt=".-",label="M="+f'{M:.2f}')

    plt.xlabel("t")
    plt.ylabel("m_eff(t)")
    plt.title("Effective Mass, N="+str(N0dat)+"x"+str(N1))
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.ylim([-1.0,4.0])
    plt.show()