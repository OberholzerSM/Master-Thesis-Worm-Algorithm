# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 10:42:03 2022

@author: samue
"""

import numpy as np
import matplotlib.pyplot as plt
modus=2
n_Gitter = 5
n_masse = 20
M_min = 2.0
M_max = 4.0

# %% Histogram Sign
ratio = np.zeros([n_Gitter,n_masse])
ratio_sd = np.zeros([n_Gitter,n_masse])
index = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
for k in range(n_Gitter):
    if modus==1:
        N0 = 2**(k+1)
        N1 = 2
    else:
        N0 = 2**n_Gitter
        N1 = 2**(k+1)
    sign_number = np.zeros((n_masse, 2))
    sign_sd = np.zeros((n_masse, 2))

    for i in range(n_masse):
        M = (i+1)*(M_max-M_min)/n_masse + M_min

        daten = np.loadtxt(
            "163Vertex_sign_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        n_data = int(sum(daten))
        sign_number[i] = daten / n_data

        daten = np.loadtxt("163Vertex_sign_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        sign_sd[i] = daten / n_data

    plt.errorbar(index,sign_number[:,0],yerr=sign_sd[:,0],fmt=".-",label="+1")
    plt.errorbar(index,sign_number[:,1],yerr=sign_sd[:,0],fmt=".-",label="-1")
    plt.plot(index,np.zeros(n_masse)+0.5,color="grey")
    plt.xlabel("M")
    plt.ylabel("Frequency")
    plt.title("Sign Frequency, N="+str(N0)+"x"+str(N1))
    plt.ylim([-0.1,1.1])
    plt.legend()
    plt.show()

    ratio[k] = (sign_number[:,0]-sign_number[:,1]) / (sign_number[:,0]+sign_number[:,1])
    ratio_sd[k] = 2.0*sign_sd[:,0]

for k in range(n_Gitter):
    if modus==1:
        N0 = 2**(k+1)
        plt.errorbar(index,ratio[k],yerr=ratio_sd[k],fmt=".-",label="N0="+str(N0))
    else:
        N1 = 2**(k+1)
        plt.errorbar(index,ratio[k],yerr=ratio_sd[k],fmt=".-",label="N1="+str(N1))
plt.plot(index,np.zeros(n_masse),color="grey")
plt.xlabel("M")
plt.ylabel("Z/Z_abs")
if modus==1:
    plt.title("Ratio Z and Z_abs, N=N0x"+str(N1))
else:
    plt.title("Ratio Z and Z_abs, N="+str(N0)+"xN1")
plt.ylim([-0.1,1.1])
plt.legend()
plt.show()

# %% 2-Punkt Korrelation Daten
N0 = 2**n_Gitter
n_data = 10**5
names=["Meson 1/2","Meson 3/4","Meson 5"]

propagator_long = np.zeros([n_Gitter,3,n_masse,2*N0-1])
propagator_long_cov = np.zeros([n_Gitter,n_masse,3*(2*N0-1),3*(2*N0-1)])
propagator_long_sd = np.zeros([n_Gitter,3,n_masse,2*N0-1])
propagator = np.zeros([n_Gitter,3,n_masse,int(N0/2)+1])
propagator_cov = np.zeros([n_Gitter,n_masse,3*(int(N0/2)+1),3*(int(N0/2)+1)])
propagator_sd = np.zeros([n_Gitter,3,n_masse,int(N0/2)+1])
m_eff = np.zeros([n_Gitter,3,n_masse,int(N0/2)])
m_eff_cov = np.zeros([n_Gitter,n_masse,3*int(N0/2),3*int(N0/2)])
m_eff_sd = np.zeros([n_Gitter,3,n_masse,int(N0/2)])

propagator_long_abs = np.zeros([n_Gitter,3,n_masse,2*N0-1])
propagator_long_cov_abs = np.zeros([n_Gitter,n_masse,3*(2*N0-1),3*(2*N0-1)])
propagator_long_sd_abs = np.zeros([n_Gitter,3,n_masse,2*N0-1])
propagator_abs = np.zeros([n_Gitter,3,n_masse,int(N0/2)+1])
propagator_cov_abs = np.zeros([n_Gitter,n_masse,3*(int(N0/2)+1),3*(int(N0/2)+1)])
propagator_sd_abs = np.zeros([n_Gitter,3,n_masse,int(N0/2)+1])
m_eff_abs = np.zeros([n_Gitter,3,n_masse,int(N0/2)])
m_eff_cov_abs = np.zeros([n_Gitter,n_masse,3*int(N0/2),3*int(N0/2)])
m_eff_sd_abs = np.zeros([n_Gitter,3,n_masse,int(N0/2)])

for k in range(n_Gitter):
    if modus==1:
        N0dat = 2**(k+1) #N0dat: Momentane Gittergrösse. N0: Maximale Gittergrösse
        N1 = 2
    else:
        N0dat = N0
        N1 = 2**(k+1)
    
    for j in range(n_masse):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        
        #Normale Daten
        propagator_long[k,0,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_Meson1_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")\
            +np.loadtxt("163Vertex_dt_number_Meson2_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        propagator_long[k,1,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_Meson3_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")\
            +np.loadtxt("163Vertex_dt_number_Meson4_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        propagator_long[k,2,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_Meson5_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        
        #Konvertiere die Kovarianzmatrix vom N0dat Format ins N0 Format
        #modus 2: daten_cov = cov
        daten_cov = \
            np.loadtxt("163Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        cov = np.zeros([5*(2*N0-1),5*(2*N0-1)])
        A = np.zeros([5*(2*N0-1),5*(2*N0dat-1)])
        for i in range(2*N0dat-1):
            for l in range(5):
                A[i+(N0-N0dat)+l*(2*N0-1),i+l*(2*N0dat-1)] = 1.0
        cov = np.matmul(np.matmul(A,daten_cov),A.transpose())
        
        #Normiere die Daten mit c(0)
        c0 = np.zeros(3)
        c0 = np.copy(propagator_long[k,:,j,N0-1])
        for i in range(3):
            propagator_long[k,i,j] =  propagator_long[k,i,j]/c0[i]

        #Kovarianz und Fehler Lange Daten für die verschiedenen Mesonen
        A = np.zeros([3*(2*N0-1),5*(2*N0-1)])
        for i in range(2*N0-1):
            A[i,i] = n_data/c0[0] #Meson 1
            A[i,i+(2*N0-1)] = n_data/c0[0]  #Meson 2
            A[i+(2*N0-1),i+2*(2*N0-1)] = n_data/c0[1]  #Meson 3
            A[i+(2*N0-1),i+3*(2*N0-1)] = n_data/c0[1]  #Meson 4
            A[i+2*(2*N0-1),i+4*(2*N0-1)] = n_data/c0[2]  #Meson 5
        propagator_long_cov[k,j] = np.matmul(np.matmul(A,cov),A.transpose())
        
        for i in range(3):
            for t in range(2*N0-1):
                propagator_long_sd[k,i,j,t] = np.sqrt(propagator_long_cov[k,j,t+i*(2*N0-1),t+i*(2*N0-1)])
            
        #Symmetrisiere die Daten
        for i in range(3):
            propagator[k,i,j,0] = propagator_long[k,i,j,N0-1]
            propagator[k,i,j,1:(int(N0dat/2)+1)] = \
                     (propagator_long[k,i,j,N0:(N0+int(N0dat/2))] \
                    + np.flip(propagator_long[k,i,j,(N0+int(N0dat/2)-1):(N0+N0dat-1)])\
                    + np.flip(propagator_long[k,i,j,(N0-int(N0dat/2)-1):(N0-1)])  \
                    + propagator_long[k,i,j,(N0-N0dat):(N0-int(N0dat/2))])/2.0
        
        #Kovarianzmatrix Symmetrisierte Daten
        A = np.zeros([3*(int(N0/2)+1),3*(2*N0-1)])
        for i in range(3):
            A[i*(int(N0/2)+1),(N0-1)+i*(2*N0-1)] = 1.0
            for t in range(int(N0/2)):
                A[(t+1)+i*(int(N0/2)+1),(N0-2-t)+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),(N0+t)+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),t+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),(2*N0-2-t)+i*(2*N0-1)] = 0.5
        propagator_cov[k,j] = np.matmul(np.matmul(A,propagator_long_cov[k,j]),A.transpose())
        
        for i in range(3):
            for t in range(int(N0/2)+1):
                propagator_sd[k,i,j,t] = np.sqrt(propagator_cov[k,j,t+i*(int(N0/2)+1),t+i*(int(N0/2)+1)])
        
        #Effektive Masse
        for i in range(3):
            m_eff[k,i,j,0:int(N0dat/2)] = -np.log(propagator[k,i,j,1:(int(N0dat/2)+1)]/propagator[k,i,j,0:int(N0dat/2)])
            
        #A Matrix: Zum Berechnen von cov_m_eff aus propagator_cov
        A = np.zeros([3*int(N0/2),3*(int(N0/2)+1)])
        for i in range(3):
            for t in range(int(N0/2)):
                A[t+i*int(N0/2),t+i*(int(N0/2)+1)] = 1.0 / propagator[k,i,j,t]
                A[t+i*int(N0/2),t+1+i*(int(N0/2)+1)] = -1.0 / propagator[k,i,j,t+1]
        m_eff_cov[k,j] = np.matmul(np.matmul(A,propagator_cov[k,j]),A.transpose())
        
        for i in range(3):
            for t in range(int(N0/2)):
                m_eff_sd[k,i,j,t] = np.sqrt(m_eff_cov[k,j,t+i*(int(N0/2)),t+i*(int(N0/2))])
        
        #abs Daten
        propagator_long_abs[k,0,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_abs_Meson1_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")\
            +np.loadtxt("163Vertex_dt_number_abs_Meson2_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        propagator_long_abs[k,1,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_abs_Meson3_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")\
            +np.loadtxt("163Vertex_dt_number_abs_Meson4_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        propagator_long_abs[k,2,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_abs_Meson5_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        
        #Konvertiere die Kovarianzmatrix vom N0dat Format ins N0 Format
        #modus 2: daten_cov = cov
        daten_cov = \
            np.loadtxt("163Vertex_dt_cov_abs_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        cov = np.zeros([5*(2*N0-1),5*(2*N0-1)])
        A = np.zeros([5*(2*N0-1),5*(2*N0dat-1)])
        for i in range(2*N0dat-1):
            for l in range(5):
                A[i+(N0-N0dat)+l*(2*N0-1),i+l*(2*N0dat-1)] = 1.0
        cov = np.matmul(np.matmul(A,daten_cov),A.transpose())
        
        #Normiere die Daten mit c(0)
        c0 = np.zeros(3)
        c0 = np.copy(propagator_long_abs[k,:,j,N0-1])
        for i in range(3):
            propagator_long_abs[k,i,j] = propagator_long_abs[k,i,j]/c0[i]

        #Kovarianz und Fehler Lange Daten für die verschiedenen Mesonen
        A = np.zeros([3*(2*N0-1),5*(2*N0-1)])
        for i in range(2*N0-1):
            A[i,i] = n_data/c0[0] #Meson 1
            A[i,i+(2*N0-1)] = n_data/c0[0]  #Meson 2
            A[i+(2*N0-1),i+2*(2*N0-1)] = n_data/c0[1]  #Meson 3
            A[i+(2*N0-1),i+3*(2*N0-1)] = n_data/c0[1]  #Meson 4
            A[i+2*(2*N0-1),i+4*(2*N0-1)] = n_data/c0[2]  #Meson 5
        propagator_long_cov_abs[k,j] = np.matmul(np.matmul(A,cov),A.transpose())
        
        for i in range(3):
            for t in range(2*N0-1):
                propagator_long_sd_abs[k,i,j,t] = np.sqrt(propagator_long_cov_abs[k,j,t+i*(2*N0-1),t+i*(2*N0-1)])
            
        #Symmetrisiere die Daten
        for i in range(3):
            propagator_abs[k,i,j,0] = propagator_long_abs[k,i,j,N0-1]
            propagator_abs[k,i,j,1:(int(N0dat/2)+1)] = \
                     (propagator_long_abs[k,i,j,N0:(N0+int(N0dat/2))] \
                    + np.flip(propagator_long_abs[k,i,j,(N0+int(N0dat/2)-1):(N0+N0dat-1)])\
                    + np.flip(propagator_long_abs[k,i,j,(N0-int(N0dat/2)-1):(N0-1)])  \
                    + propagator_long_abs[k,i,j,(N0-N0dat):(N0-int(N0dat/2))])/2.0
        
        #Kovarianzmatrix Symmetrisierte Daten
        A = np.zeros([3*(int(N0/2)+1),3*(2*N0-1)])
        for i in range(3):
            A[i*(int(N0/2)+1),(N0-1)+i*(2*N0-1)] = 1.0
            for t in range(int(N0/2)):
                A[(t+1)+i*(int(N0/2)+1),(N0-2-t)+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),(N0+t)+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),t+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),(2*N0-2-t)+i*(2*N0-1)] = 0.5
        propagator_cov_abs[k,j] = np.matmul(np.matmul(A,propagator_long_cov_abs[k,j]),A.transpose())
        
        for i in range(3):
            for t in range(int(N0/2)+1):
                propagator_sd_abs[k,i,j,t] = np.sqrt(propagator_cov_abs[k,j,t+i*(int(N0/2)+1),t+i*(int(N0/2)+1)])
        
        #Effektive Masse
        for i in range(3):
            m_eff_abs[k,i,j,0:int(N0dat/2)] =\
                -np.log(propagator_abs[k,i,j,1:(int(N0dat/2)+1)]/propagator_abs[k,i,j,0:int(N0dat/2)])
            
        #A Matrix: Zum Berechnen von cov_m_eff aus propagator_cov
        A = np.zeros([3*int(N0/2),3*(int(N0/2)+1)])
        for i in range(3):
            for t in range(int(N0/2)):
                A[t+i*int(N0/2),t+i*(int(N0/2)+1)] = 1.0 / propagator_abs[k,i,j,t]
                A[t+i*int(N0/2),t+1+i*(int(N0/2)+1)] = -1.0 / propagator_abs[k,i,j,t+1]
        m_eff_cov_abs[k,j] = np.matmul(np.matmul(A,propagator_cov_abs[k,j]),A.transpose())
        
        for i in range(3):
            for t in range(int(N0/2)):
                m_eff_sd_abs[k,i,j,t] = np.sqrt(m_eff_cov_abs[k,j,t+i*(int(N0/2)),t+i*(int(N0/2))])

# %% Masse Range für Plots
M_min_Plot = 2.0
M_max_Plot = 4.0
n_masse_Plot = 10

dM = (M_max_Plot-M_min_Plot)/n_masse_Plot
i_min = int( (M_min_Plot+dM-M_min)*n_masse/(M_max-M_min) )-1
i_max = int( (M_max_Plot-M_min)*n_masse/(M_max-M_min) )
di = int((i_max-i_min+1)/n_masse_Plot)

# %% 2-Punkt Korrelation total

for i in range(3):
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
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[(N0-N0dat):(N0+N0dat-1)],propagator_long[k,i,j,(N0-N0dat):(N0+N0dat-1)],\
                         yerr=abs(propagator_long_sd[k,i,j,(N0-N0dat):(N0+N0dat-1)]),fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("c(t)/c(0)")
        plt.title("2-point Correlation Function, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()
        
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[(N0-N0dat):(N0+N0dat-1)],propagator_long_abs[k,i,j,(N0-N0dat):(N0+N0dat-1)],\
                         yerr=abs(propagator_long_sd_abs[k,i,j,(N0-N0dat):(N0+N0dat-1)]),fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("c(t)/c(0)")
        plt.title("2-point Correlation Function abs, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()

# %% 2-Punkt Korrelation zusammengefasst

for i in range(3):
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
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[0:(int(N0dat/2)+1)],propagator[k,i,j,0:(int(N0dat/2)+1)],\
                         yerr=abs(propagator_sd[k,i,j,0:(int(N0dat/2)+1)]),fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("cS(t)/cS(0)")
        plt.title("Symmetrizied 2-point Correlation Function, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend()
        plt.show()
        
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[0:(int(N0dat/2)+1)],propagator_abs[k,i,j,0:(int(N0dat/2)+1)],\
                         yerr=abs(propagator_sd_abs[k,i,j,0:(int(N0dat/2)+1)]),fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("cS(t)/cS(0)")
        plt.title("Symmetrizied 2-point Correlation Function abs, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()

# %% Logarithmic 2-Punkt Korrelation zusammengefasst

for i in range(3):
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
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[0:(int(N0dat/2)+1)],np.log(propagator[k,i,j,0:(int(N0dat/2)+1)]),\
                         yerr=abs(propagator_sd[k,i,j,0:(int(N0dat/2)+1)]/propagator[k,i,j,0:(int(N0dat/2)+1)]),\
                             fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("ln[cS(t)/cS(0)]")
        plt.title("Logarithmic Symmetrizied 2-point Correlation Function, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()
        
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[0:(int(N0dat/2)+1)],np.log(propagator_abs[k,i,j,0:(int(N0dat/2)+1)]),\
                         yerr=abs(propagator_sd_abs[k,i,j,0:(int(N0dat/2)+1)]/propagator_abs[k,i,j,0:(int(N0dat/2)+1)])\
                             ,fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("ln[cS(t)/cS(0)]")
        plt.title("Logarithmic Symmetrizied 2-point Correlation Function abs, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()

# %% Effektive Masse

for i in range(3):
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
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[1:(int(N0dat/2))],m_eff[k,i,j,1:(int(N0dat/2))],\
                         yerr=abs(m_eff_sd[k,i,j,1:(int(N0dat/2))]),fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("m_eff(t)")
        plt.title("Effective Mass, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()
        
        for j in range(i_min,i_max,di):
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.errorbar(index[1:(int(N0dat/2))],m_eff_abs[k,i,j,1:(int(N0dat/2))],\
                         yerr=abs(m_eff_sd_abs[k,i,j,1:(int(N0dat/2))]),fmt=".-",label="M="+f'{M:.2f}')
        plt.xlabel("t")
        plt.ylabel("m_eff(t)")
        plt.title("Effective Mass abs, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()

# %% Effektive Masse Fixes Gitter
N0dat=32
N1=32

for i in range(3):
    index = np.arange(int(N0/2)-1)
    for j in range(i_min,i_max,di):
        M = (j+1)*(M_max-M_min)/n_masse + M_min
        plt.errorbar(index[1:(int(N0dat/2)-1)],m_eff[k,i,j,1:(int(N0dat/2)-1)],\
                     yerr=abs(m_eff_sd[k,i,j,1:(int(N0dat/2)-1)]),fmt=".-",label="M="+f'{M:.2f}')
    plt.plot(np.zeros(int(N0/2)-1),color="grey")
    plt.xlabel("t")
    plt.ylabel("m_eff(t)")
    plt.title("Effective Mass, "+names[i]+", N="+str(N0dat)+"x"+str(N1))
    if i==2:
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
    plt.show()