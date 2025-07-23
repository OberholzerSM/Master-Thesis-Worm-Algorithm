# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 10:13:29 2023

@author: samue
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import linalg as sp
import math 
def linearFunc(x,intercept):
    y = intercept
    return y

modus = 2 #Modus muss 2 sein
n_Gitter = 5
n_masse = 20
M_min = 2.0
M_max = 4.0
names=["Meson 1/2","Meson 3/4","Meson 5"]

u=1 #Anzahl Parameter Fit m_eff
m_eff_fit = np.zeros([3,n_Gitter,n_masse])
m_eff_fit_sd = np.zeros([3,n_Gitter,n_masse])

# %% 2-Punkt Korrelation Daten
N0 = 2**n_Gitter
n_data = 10**5
names=["Meson 1/2","Meson 3/4","Meson 5"]

propagator_long = np.zeros([3,n_Gitter,n_masse,2*N0-1])
propagator_long_cov = np.zeros([n_Gitter,n_masse,3*(2*N0-1),3*(2*N0-1)])
propagator_long_sd = np.zeros([3,n_Gitter,n_masse,2*N0-1])
propagator = np.zeros([3,n_Gitter,n_masse,int(N0/2)+1])
propagator_cov_tot = np.zeros([n_Gitter,n_masse,3*(int(N0/2)+1),3*(int(N0/2)+1)])
propagator_cov = np.zeros([3,n_Gitter,n_masse,(int(N0/2)+1),(int(N0/2)+1)])
propagator_sd = np.zeros([3,n_Gitter,n_masse,int(N0/2)+1])
m_eff = np.zeros([3,n_Gitter,n_masse,int(N0/2)])
m_eff_cov_tot = np.zeros([n_Gitter,n_masse,3*int(N0/2),3*int(N0/2)])
m_eff_cov = np.zeros([3,n_Gitter,n_masse,int(N0/2),int(N0/2)])
m_eff_sd = np.zeros([3,n_Gitter,n_masse,int(N0/2)])

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
        propagator_long[0,k,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_Meson1_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")\
            +np.loadtxt("163Vertex_dt_number_Meson2_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        propagator_long[1,k,j,(N0-N0dat):(N0+N0dat-1)] = \
            np.loadtxt("163Vertex_dt_number_Meson3_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")\
            +np.loadtxt("163Vertex_dt_number_Meson4_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        propagator_long[2,k,j,(N0-N0dat):(N0+N0dat-1)] = \
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
        c0 = np.copy(propagator_long[:,k,j,N0-1])
        for i in range(3):
            propagator_long[i,k,j] =  propagator_long[i,k,j]/c0[i]

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
                propagator_long_sd[i,k,j,t] = np.sqrt(propagator_long_cov[k,j,t+i*(2*N0-1),t+i*(2*N0-1)])
            
        #Symmetrisiere die Daten
        for i in range(3):
            propagator[i,k,j,0] = propagator_long[i,k,j,N0-1]
            propagator[i,k,j,1:(int(N0dat/2)+1)] = \
                     (propagator_long[i,k,j,N0:(N0+int(N0dat/2))] \
                    + np.flip(propagator_long[i,k,j,(N0+int(N0dat/2)-1):(N0+N0dat-1)])\
                    + np.flip(propagator_long[i,k,j,(N0-int(N0dat/2)-1):(N0-1)])  \
                    + propagator_long[i,k,j,(N0-N0dat):(N0-int(N0dat/2))])/2.0
        
        #Kovarianzmatrix Symmetrisierte Daten
        A = np.zeros([3*(int(N0/2)+1),3*(2*N0-1)])
        for i in range(3):
            A[i*(int(N0/2)+1),(N0-1)+i*(2*N0-1)] = 1.0
            for t in range(int(N0/2)):
                A[(t+1)+i*(int(N0/2)+1),(N0-2-t)+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),(N0+t)+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),t+i*(2*N0-1)] = 0.5
                A[(t+1)+i*(int(N0/2)+1),(2*N0-2-t)+i*(2*N0-1)] = 0.5
        propagator_cov_tot[k,j] = np.matmul(np.matmul(A,propagator_long_cov[k,j]),A.transpose())
        
        for i in range(3):
            propagator_cov[i,k,j] = propagator_cov_tot[k,j,(i*(int(N0/2)+1)):(i+1)*(int(N0/2)+1),(i*(int(N0/2)+1)):(i+1)*(int(N0/2)+1)]
            for t in range(int(N0/2)+1):
                propagator_sd[i,k,j,t] = np.sqrt(propagator_cov_tot[k,j,t+i*(int(N0/2)+1),t+i*(int(N0/2)+1)])
        
        #Effektive Masse
        for i in range(3):
            m_eff[i,k,j,0:int(N0dat/2)] = -np.log(propagator[i,k,j,1:(int(N0dat/2)+1)]/propagator[i,k,j,0:int(N0dat/2)])
            
        #A Matrix: Zum Berechnen von cov_m_eff aus propagator_cov_tot
        A = np.zeros([3*int(N0/2),3*(int(N0/2)+1)])
        for i in range(3):
            for t in range(int(N0/2)):
                A[t+i*int(N0/2),t+i*(int(N0/2)+1)] = 1.0 / propagator[i,k,j,t]
                A[t+i*int(N0/2),t+1+i*(int(N0/2)+1)] = -1.0 / propagator[i,k,j,t+1]
        m_eff_cov_tot[k,j] = np.matmul(np.matmul(A,propagator_cov_tot[k,j]),A.transpose())
        
        for i in range(3):
            m_eff_cov[i,k,j] = propagator_cov_tot[k,j,(i*int(N0/2)):((i+1)*int(N0/2)),(i*int(N0/2)):((i+1)*int(N0/2))]
            for t in range(int(N0/2)):
                m_eff_sd[i,k,j,t] = np.sqrt(m_eff_cov_tot[k,j,t+i*(int(N0/2)),t+i*(int(N0/2))])
        
       
#%% Plotte m_eff für eine gegebene Gittergrösse und Masse und bestimme t_max

k = 1
N0 = 32
N1 = 2**5
M = 2.2

i = int(np.log(N1)/np.log(2))-1
j = int((M-M_min)*n_masse/(M_max-M_min))-1
index = np.arange(int(N0/2))
plt.errorbar(index,m_eff[k,i,j],yerr=m_eff_sd[k,i,j],fmt=".-",label="M="+f'{M:.2f}')
plt.xlabel("t")
plt.ylabel("m_eff(t)")
plt.title("Effective Mass, "+names[k]+", N="+str(N0)+"x"+str(N1))
plt.legend(loc="upper right")
plt.show()

#Bestimme das grösste erlaubte tmax
if len(np.where(m_eff[k,i,j] > 10**9)[0])==0:
    tmax_tot = int(N0/2)-1
else:
    tmax_tot = np.where(m_eff[k,i,j] > 10**9)[0][0]-1

tmin_tot = 1
test = 0
while test==0:
    #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
    if math.isnan(m_eff_sd[k,i,j,tmax_tot]) or \
        m_eff_sd[k,i,j,tmax_tot] < np.finfo(float).eps or \
        len(np.where(np.isnan(m_eff[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0 or \
        len(np.where(np.isnan(m_eff_sd[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
        tmax_tot = tmax_tot-1
    else:
        test=1
print("tmax =",tmax_tot)

#%% Bestimme m_eff_fit für verschiedene Bereiche von t und mache ein gewichtetes Mittel via AIC

index = np.arange(int(N0/2))
tmin_tot = 1

for k in range(3):
    for i in range(n_Gitter):
        for j in range(n_masse):
            
            #Bestimme das grösste erlaubte tmax
            if len(np.where(m_eff[k,i,j] > 10**9)[0])==0:
                tmax_tot = int(N0/2)-1
            else:
                tmax_tot = np.where(m_eff[k,i,j] > 10**9)[0][0]-1

            tmin_tot = 1
            test = 0
            while test==0:
                #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
                if math.isnan(m_eff_sd[k,i,j,tmax_tot]) or \
                    (np.any(np.linalg.eigvals(m_eff_cov[k,i,j,tmin_tot:(tmax_tot+1),tmin_tot:(tmax_tot+1)]) \
                    < np.finfo(float).eps)) or \
                    len(np.where(np.isnan(m_eff[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0 or \
                    m_eff_sd[k,i,j,tmax_tot] < np.finfo(float).eps or \
                    len(np.where(np.isnan(m_eff[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0 or \
                    len(np.where(np.isnan(m_eff_sd[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
                    tmax_tot = tmax_tot-1
                else:
                    test=1
            
            #Falls tmax_tot gross genug ist, mache den Fit
            if tmax_tot > u+1+tmin_tot:
                n = tmax_tot - tmin_tot
                m_eff_prop = np.zeros(int(0.5*(n-2)*(n-1)))
                AIC = np.zeros(int(0.5*(n-2)*(n-1)))
                weights = np.zeros(int(0.5*(n-2)*(n-1)))
                
                counter=0
                for tmin in range(tmin_tot,tmax_tot-u-1):
                    for tmax in range(tmin+u+2,tmax_tot+1):
                        n=tmax-tmin
                        cov_reduced = m_eff_cov[k,i,j,tmin:(tmax+1),tmin:(tmax+1)]
                        a_fit,cov_fit=curve_fit(linearFunc,index[tmin:(tmax+1)],m_eff[k,i,j,tmin:(tmax+1)],\
                                            sigma=cov_reduced,absolute_sigma=True)
                        theorie = np.zeros(tmax-tmin+1)+a_fit[0]
                        x = m_eff[k,i,j,tmin:(tmax+1)]-theorie
                        logL = np.matmul( np.matmul(x.transpose(),sp.inv(cov_reduced)),x)/(n-u)
                        
                        AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                        m_eff_prop[counter] = a_fit[0]
                        counter += 1
                
                AIC = AIC - min(AIC)
                weights = np.exp(-0.5*AIC)
                
                m_eff_fit[k,i,j] = sum(weights*m_eff_prop)/sum(weights)
                m_eff_fit_sd[k,i,j] = np.sqrt( sum(weights*(m_eff_fit[k,i,j]-m_eff_prop)**2)/sum(weights) )
            else:
                m_eff_fit[k,i,j] = float('NaN')
                m_eff_fit_sd[k,i,j] = float('NaN')

# Plotte m_eff_fit, sobald man alle Daten hat
for k in range(3):
    for i in range(n_Gitter):
        N1=2**(i+1)
        index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
        plt.errorbar(index_Masse,m_eff_fit[k,i],yerr=m_eff_fit_sd[k,i],fmt=".-",label=("N1="+str(N1)))
    plt.xlabel("M_bar")
    plt.ylabel("m_eff")
    plt.title("Effective Mass vs. Bar Mass, "+names[k]+", N="+str(N0)+"xN1")
    plt.legend() 
    plt.show()

#%% Bestimme m_eff_fit für verschiedene Bereiche von t, ignoriere die Kovarianz

index = np.arange(int(N0/2))
tmin_tot = 1

for k in range(3):
    for i in range(n_Gitter):
        for j in range(n_masse):
            
            #Bestimme das grösste erlaubte tmax
            if len(np.where(m_eff[k,i,j] > 10**9)[0])==0:
                tmax_tot = int(N0/2)-1
            else:
                tmax_tot = np.where(m_eff[k,i,j] > 10**9)[0][0]-1

            tmin_tot = 1
            test = 0
            while test==0:
                #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
                if math.isnan(m_eff_sd[k,i,j,tmax_tot]) or \
                    m_eff_sd[k,i,j,tmax_tot] < np.finfo(float).eps or \
                    len(np.where(np.isnan(m_eff[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0 or \
                    len(np.where(np.isnan(m_eff_sd[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
                    tmax_tot = tmax_tot-1
                else:
                    test=1
            
            #Falls tmax_tot gross genug ist, mache den Fit
            if tmax_tot > u+1+tmin_tot:
                n = tmax_tot - tmin_tot
                m_eff_prop = np.zeros(int(0.5*(n-2)*(n-1)))
                AIC = np.zeros(int(0.5*(n-2)*(n-1)))
                weights = np.zeros(int(0.5*(n-2)*(n-1)))
                
                counter=0
                for tmin in range(tmin_tot,tmax_tot-u-1):
                    for tmax in range(tmin+u+2,tmax_tot+1):
                        n=tmax-tmin
                        a_fit,cov_fit=curve_fit(linearFunc,index[tmin:(tmax+1)],m_eff[k,i,j,tmin:(tmax+1)],\
                                            sigma=m_eff_sd[k,i,j,tmin:(tmax+1)],absolute_sigma=True,bounds=(0.0, np.inf))
                        theorie = np.zeros(tmax-tmin+1)+a_fit[0]
                        x = m_eff[k,i,j,tmin:(tmax+1)]-theorie
                        logL = sum( (x-theorie)**2/m_eff_sd[k,i,j,tmin:(tmax+1)] )/(n-u)
                        
                        AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                        m_eff_prop[counter] = a_fit[0]
                        counter += 1
                
                AIC = AIC - min(AIC)
                weights = np.exp(-0.5*AIC)
                
                m_eff_fit[k,i,j] = sum(weights*m_eff_prop)/sum(weights)
                m_eff_fit_sd[k,i,j] = np.sqrt( sum(weights*(m_eff_fit[k,i,j]-m_eff_prop)**2)/sum(weights) )
            else:
                m_eff_fit[k,i,j] = float('NaN')
                m_eff_fit_sd[k,i,j] = float('NaN')

# Plotte m_eff_fit, sobald man alle Daten hat
for k in range(3):
    for i in range(n_Gitter):
        N1=2**(i+1)
        index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
        plt.errorbar(index_Masse,m_eff_fit[k,i],yerr=m_eff_fit_sd[k,i],fmt=".-",label=("N1="+str(N1)))
    plt.xlabel("M_bar")
    plt.ylabel("m_eff")
    plt.title("Effective Mass vs. Bar Mass, "+names[k]+", N="+str(N0)+"xN1")
    plt.legend() 
    plt.show()

#%% Bestimme den kritischen Exponenten

index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
daten = m_eff_fit[1,n_Gitter-1]
daten_sd = m_eff_fit_sd[1,n_Gitter-1]

def polynom(M_bar,M_crit,A,a_crit):
    y = A*(M_bar-M_crit)**a_crit
    return y

a_fit,cov_fit = curve_fit(polynom,index_Masse,daten,bounds=(0.0, np.inf))
print("M_crit=",a_fit[0],"+-",np.sqrt(cov_fit[0,0]))
print("A=",a_fit[1],"+-",np.sqrt(cov_fit[1,1]))
print("a_crit=",a_fit[2],"+-",np.sqrt(cov_fit[2,2]))


#%% m_fit Parameter

u=2 #Anzahl Parameter Fit
N0 = 2**n_Gitter
def cosh(t,mp,c):
    A = c / (1 + np.exp(-mp*N0))
    B = c / (1 + np.exp(mp*N0))
    y = A*np.exp(-mp*t) + B*np.exp(mp*t)
    return y
m_fit = np.zeros([3,n_Gitter,n_masse])
m_fit_sd = np.zeros([3,n_Gitter,n_masse])
c0_fit = np.zeros([3,n_Gitter,n_masse])
c0_fit_sd = np.zeros([3,n_Gitter,n_masse])

#%% Plotte cs(t) für eine gegebene Gittergrösse und Masse und bestimme t_max

k=1
N0 = 32
N1 = 32
M = 2.1

i = int(np.log(N1)/np.log(2))-1
j = int((M-M_min)*n_masse/(M_max-M_min))-1
index = np.arange(int(N0/2)+1)
plt.errorbar(index,propagator[k,i,j],yerr=propagator_sd[k,i,j],fmt=".-",label="M="+f'{M:.2f}')
plt.xlabel("t")
plt.ylabel("cS(t)/cS(0)")
plt.title("cS(t), N="+str(N0)+"x"+str(N1))
plt.legend(loc="upper right")
plt.show()

#Bestimme das grösste erlaubte tmax
if len(np.where(propagator[k,i,j] > 10**9)[0])==0:
    tmax_tot = int(N0/2)-1
else:
    tmax_tot = np.where(propagator[k,i,j] > 10**9)[0][0]-1

tmin_tot = 1
test = 0
while test==0:
    #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
    if math.isnan(propagator_sd[k,i,j,tmax_tot]) or \
        (np.any(np.linalg.eigvals(propagator_cov[k,i,j,tmin_tot:(tmax_tot+1),tmin_tot:(tmax_tot+1)]) \
        < np.finfo(float).eps)) or \
        len(np.where(np.isnan(propagator[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
        tmax_tot = tmax_tot-1
    else:
        test=1
print("tmax = ",tmax_tot)

a_fit,cov_fit=curve_fit(cosh,index[tmin_tot:(tmax_tot+1)],propagator[k,i,j,tmin_tot:(tmax_tot+1)],\
                    p0=[0.0,1.0],\
                    sigma=propagator_sd[k,i,j,tmin_tot:(tmax_tot+1)],absolute_sigma=True,\
                    bounds=(0.0, np.inf))
print(a_fit[0])
#%% Bestimme m_fit für verschiedene Bereiche von t und mache ein gewichtetes Mittel via AIC (ignoriere Kovarianz)

index = np.arange(int(N0/2))
tmin_tot = 1

for k in range(3):
    for i in range(n_Gitter):
        for j in range(n_masse):
            
            #Bestimme das grösste erlaubte tmax
            if len(np.where(propagator[k,i,j] > 10**9)[0])==0:
                tmax_tot = int(N0/2)-1
            else:
                tmax_tot = np.where(propagator[k,i,j] > 10**9)[0][0]-1

            tmin_tot = 1
            test = 0
            while test==0:
                #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
                if math.isnan(propagator_sd[k,i,j,tmax_tot]) or \
                    propagator_sd[k,i,j,tmax_tot] < np.finfo(float).eps or \
                    len(np.where(np.isnan(propagator[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0 or \
                    len(np.where(np.isnan(propagator_sd[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
                    tmax_tot = tmax_tot-1
                else:
                    test=1
            
            #Falls tmax_tot gross genug ist, mache den Fit
            if tmax_tot > u+1+tmin_tot:
                n = tmax_tot - tmin_tot
                m_prop = np.zeros(int(0.5*(n-2)*(n-1)))
                c0_prop = np.zeros(int(0.5*(n-2)*(n-1)))
                AIC = np.zeros(int(0.5*(n-2)*(n-1)))
                weights = np.zeros(int(0.5*(n-2)*(n-1)))
                
                counter=0
                for tmin in range(tmin_tot,tmax_tot-u-1):
                    for tmax in range(tmin+u+2,tmax_tot+1):
                        n=tmax-tmin
                        if math.isnan(m_eff_fit[k,i,j]):
                            a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[k,i,j,tmin:(tmax+1)],\
                                                p0=[0.0,1.0],\
                                                sigma=propagator_sd[k,i,j,tmin:(tmax+1)],absolute_sigma=True,\
                                                bounds=(0.0, np.inf))
                        else:
                            a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[k,i,j,tmin:(tmax+1)],\
                                                p0=[m_eff_fit[k,i,j],1.0],\
                                                sigma=propagator_sd[k,i,j,tmin:(tmax+1)],absolute_sigma=True,\
                                                bounds=(0.0, np.inf))
                        theorie = cosh(index[tmin:(tmax+1)],a_fit[0],a_fit[1])
                        x = propagator[k,i,j,tmin:(tmax+1)]-theorie
                        logL = sum( (x-theorie)**2/propagator_sd[k,i,j,tmin:(tmax+1)] )/(n-u)
                        
                        AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                        m_prop[counter] = a_fit[0]
                        c0_prop[counter] = a_fit[1]
                        counter += 1
                
                AIC = AIC - min(AIC)
                weights = np.exp(-0.5*AIC)
                
                m_fit[k,i,j] = sum(weights*m_prop)/sum(weights)
                m_fit_sd[k,i,j] = np.sqrt( sum(weights*(m_fit[k,i,j]-m_prop)**2)/sum(weights) )
                c0_fit[k,i,j] = sum(weights*c0_prop)/sum(weights)
                c0_fit_sd[k,i,j] = np.sqrt( sum(weights*(c0_fit[k,i,j]-c0_prop)**2)/sum(weights) )
            else:
                m_eff_fit[k,i,j] = float('NaN')
                m_eff_fit_sd[k,i,j] = float('NaN')

# Plotte m_eff_fit, sobald man alle Daten hat
for k in range(3):
    for i in range(n_Gitter):
        N1=2**(i+1)
        index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
        plt.errorbar(index_Masse,m_fit[k,i],yerr=m_fit_sd[k,i],fmt=".-",label=("N1="+str(N1)))
    plt.xlabel("M_bar")
    plt.ylabel("m_p")
    plt.title("Physical Mass vs. Bar Mass, "+names[k]+", N="+str(N0)+"xN1")
    plt.legend() 
    plt.show()

#%% Bestimme m_fit für verschiedene Bereiche von t und mache ein gewichtetes Mittel via AIC (ignoriere Fehler)

index = np.arange(int(N0/2))
tmin_tot = 1

for k in range(3):
    for i in range(n_Gitter):
        for j in range(n_masse):
            
            #Bestimme das grösste erlaubte tmax
            if len(np.where(propagator[k,i,j] > 10**9)[0])==0:
                tmax_tot = int(N0/2)-1
            else:
                tmax_tot = np.where(propagator[k,i,j] > 10**9)[0][0]-1

            tmin_tot = 1
            test = 0
            while test==0:
                #Falls ungültige Messwerte vorliegen, reduziere t_max weiter
                if len(np.where(np.isnan(propagator[k,i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
                    tmax_tot = tmax_tot-1
                else:
                    test=1
            
            #Falls tmax_tot gross genug ist, mache den Fit
            if tmax_tot > u+1+tmin_tot:
                n = tmax_tot - tmin_tot
                m_prop = np.zeros(int(0.5*(n-2)*(n-1)))
                c0_prop = np.zeros(int(0.5*(n-2)*(n-1)))
                AIC = np.zeros(int(0.5*(n-2)*(n-1)))
                weights = np.zeros(int(0.5*(n-2)*(n-1)))
                
                counter=0
                for tmin in range(tmin_tot,tmax_tot-u-1):
                    for tmax in range(tmin+u+2,tmax_tot+1):
                        n=tmax-tmin
                        if math.isnan(m_eff_fit[k,i,j]):
                            a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[k,i,j,tmin:(tmax+1)],\
                                                p0=[0.0,1.0],\
                                                bounds=(0.0, np.inf))
                        else:
                            a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[k,i,j,tmin:(tmax+1)],\
                                                p0=[m_eff_fit[k,i,j],1.0],\
                                                bounds=(0.0, np.inf))
                        theorie = cosh(index[tmin:(tmax+1)],a_fit[0],a_fit[1])
                        x = propagator[k,i,j,tmin:(tmax+1)]-theorie
                        logL = sum((x-theorie)**2)/(n-u)
                        
                        AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                        m_prop[counter] = a_fit[0]
                        c0_prop[counter] = a_fit[1]
                        counter += 1
                
                AIC = AIC - min(AIC)
                weights = np.exp(-0.5*AIC)
                
                m_fit[k,i,j] = sum(weights*m_prop)/sum(weights)
                m_fit_sd[k,i,j] = np.sqrt( sum(weights*(m_fit[k,i,j]-m_prop)**2)/sum(weights) )
                c0_fit[k,i,j] = sum(weights*c0_prop)/sum(weights)
                c0_fit_sd[k,i,j] = np.sqrt( sum(weights*(c0_fit[k,i,j]-c0_prop)**2)/sum(weights) )
            else:
                m_eff_fit[k,i,j] = float('NaN')
                m_eff_fit_sd[k,i,j] = float('NaN')

# Plotte m_eff_fit, sobald man alle Daten hat
for k in range(3):
    for i in range(n_Gitter):
        N1=2**(i+1)
        index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
        plt.errorbar(index_Masse,m_fit[k,i],yerr=m_fit_sd[k,i],fmt=".-",label=("N1="+str(N1)))
    plt.xlabel("M_bar")
    plt.ylabel("m_p")
    plt.title("Physical Mass vs. Bar Mass, "+names[k]+", N="+str(N0)+"xN1")
    plt.legend() 
    plt.show()

#%% Plotte c(0), sobald man alle Daten hat
for k in range(3):
    for i in range(n_Gitter):
        N1=2**(i+1)
        index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
        plt.errorbar(index_Masse,c0_fit[k,i],yerr=c0_fit_sd[k,i],fmt=".-",label=("N1="+str(N1)))
    plt.xlabel("M_bar")
    plt.ylabel("m_p")
    plt.title("Normalization Constant, "+names[k]+", N="+str(N0)+"xN1")
    plt.legend() 
    plt.show()

#%% Bestimme den kritischen Exponenten

index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
daten = m_fit[1,n_Gitter-1]

def polynom(M_bar,M_crit,A,a_crit):
    y = A*(M_bar-M_crit)**a_crit
    return y

a_fit,cov_fit = curve_fit(polynom,index_Masse,daten,bounds=(0.0, np.inf))
print("M_crit=",a_fit[0],"+-",np.sqrt(cov_fit[0,0]))
print("A=",a_fit[1],"+-",np.sqrt(cov_fit[1,1]))
print("a_crit=",a_fit[2],"+-",np.sqrt(cov_fit[2,2]))

#%% Teste, ob die linearisierten Daten Sinn ergeben

def invcosh(x,m,C):
    a = C / (1.0 + np.exp(-m*N0))
    b = C / (1.0 + np.exp(m*N0))
    if b > 0:
        d = x / (2.0*b)
        e = (x**2/(4.0*b**2)) - (a/b)
    else:
        d = 0.0
        e = 0.0
    if e>0.0 and b>0.0 and (d+np.sqrt(e)) > 0.0:
        return np.log( d + np.sqrt(e) )
    else:
        return float("NaN")

def linearFunc2(x,slope,intercept):
    y = slope*x+intercept
    return y

m_linear_fit = np.zeros([3,n_Gitter,n_masse])
propagator_linear = np.zeros([3,n_Gitter,n_masse,int(N0/2)+1])
index = np.arange(int(N0/2)+1)
index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
for k in range(3):
    for i in range(n_Gitter):
        N1 = 2**(i+1)
        for j in range(n_masse):
            for t in range(int(N0/2)+1):
                propagator_linear[k,i,j,t] = invcosh(propagator[k,i,j,t],m_fit[k,i,j],c0_fit[k,i,j])
            a_fit,cov_fit=curve_fit(linearFunc2,index[1:(int(N0/2)+1)],propagator[k,i,j,1:(int(N0/2)+1)],\
                                bounds=(0.0, np.inf))
            m_linear_fit[k,i,j] = a_fit[0]
        plt.errorbar(index_Masse,m_linear_fit[k,i]-m_fit[k,i],yerr=m_fit_sd[k,i],fmt=".-",label=("N1="+str(N1)))
    plt.xlabel("M_bar")
    plt.ylabel("Deviation")
    plt.title("Deviation between the linearized Data and the cosh Fit, "+names[k])
    plt.show()

#%% Plotte linearisierte Daten

def invcosh(x,m,C):
    a = C / (1.0 + np.exp(-m*N0))
    b = C / (1.0 + np.exp(m*N0))
    d = x / (2.0*b)
    e = (x**2/(4.0*b**2)) - (a/b)
    if e>0.0 and (d+np.sqrt(e)) > 0:
        return np.log( d + np.sqrt(e) )
    else:
        return float("NaN")

propagator_linear = np.zeros([3,n_Gitter,n_masse,int(N0/2)+1])
index = np.arange(int(N0/2)+1)
for k in range(3):
    for i in range(n_Gitter):
        N1 = 2**(i+1)
        for j in range(n_masse):
            for t in range(int(N0/2)+1):
                propagator_linear[k,i,j,t] = invcosh(propagator[k,i,j,t],m_fit[k,i,j],c0_fit[k,i,j])
            M = (j+1)*(M_max-M_min)/n_masse + M_min
            plt.scatter(index,propagator_linear[k,i,j],label="M="+f'{M:.2f}')
            plt.xlabel("t")
            plt.ylabel("F[cS(t)]")
            plt.title("Linearized Data, "+names[k]+", N="+str(N0)+"x"+str(N1))
            plt.legend(loc="upper right")
            plt.show()