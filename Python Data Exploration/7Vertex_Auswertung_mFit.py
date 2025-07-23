# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 15:18:57 2023

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
            propagator_long_sd[k,i,(N0-N0dat):(N0+N0dat-1)] = np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
            cov_long[k,i,(N0-N0dat):(N0+N0dat-1),(N0-N0dat):(N0+N0dat-1)] = np.loadtxt("7Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0dat)+"x"+str(N1)+".dat")
        else:
            propagator_long[k,i] = np.loadtxt("7Vertex_dt_number_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
            propagator_long_sd[k,i] = np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
            cov_long[k,i] = np.loadtxt("7Vertex_dt_cov_M="+f'{M:.2f}'+",N="+str(N0)+"x"+str(N1)+".dat")
        
        c0 = propagator_long[k,i,N0-1]
        propagator_long[k,i] = propagator_long[k,i]/c0
        #propagator_long_sd[k,i] = propagator_long_sd[k,i]/c0
        
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
#%% Plotte m_eff für eine gegebene Gittergrösse und Masse und bestimme t_max

N0 = 32
N1 = 32
M = 2.8

i = int(np.log(N1)/np.log(2))-1
j = int((M-M_min)*n_masse/(M_max-M_min))-1
index = np.arange(int(N0/2))
plt.errorbar(index,m_eff[i,j],yerr=m_eff_sd[i,j],fmt=".-",label="M="+f'{M:.2f}')
plt.xlabel("t")
plt.ylabel("m_eff(t)")
plt.title("Effective Mass, N="+str(N0)+"x"+str(N1))
plt.legend(loc="upper right")
plt.show()

#Bestimme das grösste erlaubte tmax
if len(np.where(m_eff[i,j] > 10**9)[0])==0:
    tmax_tot = int(N0/2)-1
else:
    tmax_tot = np.where(m_eff[i,j] > 10**9)[0][0]-1

tmin_tot = 3
test = 0
while test==0:
    #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
    if math.isnan(m_eff_sd[i,j,tmax_tot]) or \
        (np.any(np.linalg.eigvals(cov_m_eff[i,j,tmin_tot:(tmax_tot+1),tmin_tot:(tmax_tot+1)]) < np.finfo(float).eps)):
        tmax_tot = tmax_tot-1
    else:
        test=1
print("tmax =",tmax_tot)

#Bestimme das grösste erlaubte tmax
if len(np.where(m_eff[i,j] > 10**9)[0])==0:
    tmax_tot = int(N0/2)-1
else:
    tmax_tot = np.where(m_eff[i,j] > 10**9)[0][0]-1

tmin_tot = 1
test = 0
while test==0:
    #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
    if math.isnan(m_eff_sd[i,j,tmax_tot]) or \
        len(np.where(np.isnan(m_eff_sd[i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
        tmax_tot = tmax_tot-1
    else:
        test=1
print("tmax =",tmax_tot)

#%% Bestimme m_eff_fit für verschiedene Bereiche von t und mache ein gewichtetes Mittel via AIC

index = np.arange(int(N0/2))
tmin_tot = 1

for i in range(n_Gitter):
    for j in range(n_masse):
        
        #Bestimme das grösste erlaubte tmax
        if len(np.where(m_eff[i,j] > 10**9)[0])==0:
            tmax_tot = int(N0/2)-1
        else:
            tmax_tot = np.where(m_eff[i,j] > 10**9)[0][0]-1
        
        test = 0
        while test==0:
            #Falls die Kovarianzmatrix nicht positiv definit ist, reduziere tmax_tot weiter
            if math.isnan(m_eff_sd[i,j,tmax_tot]) or \
                len(np.where(np.isnan(m_eff_sd[i,j,tmin_tot:(tmax_tot+1)]))[0])!=0 or \
                (np.any(np.linalg.eigvals(cov_m_eff[i,j,tmin_tot:(tmax_tot+1),tmin_tot:(tmax_tot+1)]) < np.finfo(float).eps)):
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
                    a_fit,cov_fit=curve_fit(linearFunc,index[tmin:(tmax+1)],m_eff[i,j,tmin:(tmax+1)],\
                                        sigma=cov_m_eff[i,j,tmin:(tmax+1),tmin:(tmax+1)],absolute_sigma=True)
                    theorie = np.zeros(tmax-tmin+1)+a_fit[0]
                    x = m_eff[i,j,tmin:(tmax+1)]-theorie
                    cov_reduced = cov_m_eff[i,j,tmin:(tmax+1),tmin:(tmax+1)]
                    logL = np.matmul( np.matmul(x.transpose(),sp.inv(cov_reduced)),x)/(n-u)
                    
                    AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                    m_eff_prop[counter] = a_fit[0]
                    counter += 1
            
            AIC = AIC - min(AIC)
            weights = np.exp(-0.5*AIC)
            
            m_eff_fit[i,j] = sum(weights*m_eff_prop)/sum(weights)
            m_eff_fit_sd[i,j] = np.sqrt( sum(weights*(m_eff_fit[i,j]-m_eff_prop)**2)/sum(weights) )
        else:
            #Bestimme das grösste erlaubte tmax
            if len(np.where(m_eff[i,j] > 10**9)[0])==0:
                tmax_tot = int(N0/2)-1
            else:
                tmax_tot = np.where(m_eff[i,j] > 10**9)[0][0]-1
            
            test = 0
            while test==0:
                if math.isnan(m_eff_sd[i,j,tmax_tot]) or \
                    len(np.where(np.isnan(m_eff_sd[i,j,tmin_tot:(tmax_tot+1)]))[0])!=0:
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
                        a_fit,cov_fit=curve_fit(linearFunc,index[tmin:(tmax+1)],m_eff[i,j,tmin:(tmax+1)],\
                                            sigma=m_eff_sd[i,j,tmin:(tmax+1)],absolute_sigma=True)
                        theorie = np.zeros(tmax-tmin+1)+a_fit[0]
                        x = m_eff[i,j,tmin:(tmax+1)]-theorie
                        logL = sum( (theorie-x)**2 / m_eff_sd[i,j,tmin:(tmax+1)] )/(n-u)
                        
                        AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                        m_eff_prop[counter] = a_fit[0]
                        counter += 1
                
                AIC = AIC - min(AIC)
                weights = np.exp(-0.5*AIC)
                
                m_eff_fit[i,j] = sum(weights*m_eff_prop)/sum(weights)
                m_eff_fit_sd[i,j] = np.sqrt( sum(weights*(m_eff_fit[i,j]-m_eff_prop)**2)/sum(weights) )
                
            else: 
                m_eff_fit[i,j] = float('NaN')
                m_eff_fit_sd[i,j] = float('NaN')

# Plotte m_eff_fit, sobald man alle Daten hat

for k in range(n_Gitter):
    N1=2**(k+1)
    index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
    plt.errorbar(index_Masse,m_eff_fit[k],yerr=m_eff_fit_sd[k],fmt=".-",label=("N1="+str(N1)))
plt.xlabel("M_bar")
plt.ylabel("m_eff")
plt.title("Effective Mass vs. Bar Mass, N="+str(N0)+"xN1")
plt.legend() 
plt.show()



#%% Bestimme den kritischen Exponenten

index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
a = 13
b = 17
daten = m_eff_fit[n_Gitter-1,a:b]

def polynom(M_bar,M_crit,A,a_crit):
    y = A*(M_bar-M_crit)**a_crit
    return y

a_fit,cov_fit = curve_fit(polynom,index_Masse[a:b],daten,bounds=(0.0, np.inf))
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
m_fit = np.zeros([n_Gitter,n_masse])
m_fit_sd = np.zeros([n_Gitter,n_masse])
c0_fit = np.zeros([n_Gitter,n_masse])
c0_fit_sd = np.zeros([n_Gitter,n_masse])

#%% Plotte cs(t) für eine gegebene Gittergrösse und Masse und bestimme t_max

N0 = 32
N1 = 32
M = 0.6

i = int(np.log(N1)/np.log(2))-1
j = int((M-M_min)*n_masse/(M_max-M_min))-1
index = np.arange(int(N0/2)+1)
plt.errorbar(index,propagator[i,j],yerr=propagator_sd[i,j],fmt=".-",label="M="+f'{M:.2f}')
plt.xlabel("t")
plt.ylabel("cS(t)/cS(0)")
plt.title("cS(t), N="+str(N0)+"x"+str(N1))
plt.legend(loc="upper right")
plt.show()

#Bestimme das grösste erlaubte tmax
if len(np.where(propagator[i,j]==0)[0])==0:
    tmax_tot = int(N0/2)
else:
    tmax_tot = np.where(propagator[i,j]==0)[0][0]-1

tmin_tot = 1
test = 0
while test==0:
    if math.isnan(propagator_sd[i,j,tmax_tot]) or\
        (np.any(np.linalg.eigvals(cov[i,j,tmin_tot:(tmax_tot+1),tmin_tot:(tmax_tot+1)]) < np.finfo(float).eps)):
        tmax_tot = tmax_tot-1
    else:
        test=1
print("t_range für positiv definite Kovarianz: ",[tmin_tot,tmax_tot])

#Bestimme das grösste erlaubte tmax ohne Kovarianz
if len(np.where(propagator[i,j]==0)[0])==0:
    tmax_tot = int(N0/2)
else:
    tmax_tot = np.where(propagator[i,j]==0)[0][0]-1
while test==0:
    if math.isnan(propagator_sd[i,j,tmax_tot]):
        tmax_tot = tmax_tot-1
    else:
        test=1
print("t_range für valide Daten ohne Kovarianz",[tmin_tot,tmax_tot])

#%% Bestimme m_fit für verschiedene Bereiche von t und mache ein gewichtetes Mittel via AIC

index = np.arange(int(N0/2)+1)
tmin_tot = 1

for i in range(n_Gitter):
    for j in range(n_masse):
        #Bestimme das grösste erlaubte tmax
        if len(np.where(propagator[i,j]==0)[0])==0:
            tmax_tot = int(N0/2)
        else:
            tmax_tot = np.where(propagator[i,j]==0)[0][0]-1

        test = 0
        while test==0:
            if math.isnan(propagator_sd[i,j,tmax_tot]) or\
                (np.any(np.linalg.eigvals(cov[i,j,tmin_tot:(tmax_tot+1),tmin_tot:(tmax_tot+1)]) < np.finfo(float).eps)):
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
                    if math.isnan(m_eff_fit[i,j]):
                        a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[i,j,tmin:(tmax+1)],\
                                            p0=[0.0,1.0],\
                                            sigma=cov[i,j,tmin:(tmax+1),tmin:(tmax+1)],absolute_sigma=True,\
                                            bounds=(0.0, np.inf))
                    else:
                        a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[i,j,tmin:(tmax+1)],\
                                            p0=[m_eff_fit[i,j],1.0],\
                                            sigma=cov[i,j,tmin:(tmax+1),tmin:(tmax+1)],absolute_sigma=True,\
                                            bounds=(0.0, np.inf))
                    
                    
                    theorie = cosh(index[tmin:(tmax+1)],a_fit[0],a_fit[1])
                    x = propagator[i,j,tmin:(tmax+1)]-theorie
                    cov_reduced = cov[i,j,tmin:(tmax+1),tmin:(tmax+1)]
                    logL = np.matmul( np.matmul(x.transpose(),sp.inv(cov_reduced)),x)/(n-u)
                    
                    AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                    m_prop[counter] = a_fit[0]
                    c0_prop[counter] = a_fit[1]
                    counter += 1
            
            AIC = AIC - min(AIC)
            weights = np.exp(-0.5*AIC)
            
            m_fit[i,j] = sum(weights*m_prop)/sum(weights)
            m_fit_sd[i,j] = np.sqrt( sum(weights*(m_fit[i,j]-m_prop)**2)/sum(weights) )
            c0_fit[i,j] = sum(weights*c0_prop)/sum(weights)
            c0_fit_sd[i,j] = np.sqrt( sum(weights*(c0_fit[i,j]-c0_prop)**2)/sum(weights) )
        else:
            m_fit[i,j] = float('NaN')
            m_fit_sd[i,j] = float('NaN')
            c0_fit[i,j] = float('NaN')
            c0_fit_sd[i,j] = float('NaN')

# Plotte m_fit, sobald man alle Daten hat

for k in range(n_Gitter):
    N1=2**(k+1)
    index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
    plt.errorbar(index_Masse,m_fit[k],yerr=m_fit_sd[k],fmt=".-",label=("N1="+str(N1)))
plt.xlabel("M_bar")
plt.ylabel("m_fit")
plt.title("Physical Mass vs. Bar Mass, N="+str(N0)+"xN1")
plt.legend() 
plt.show()

for k in range(n_Gitter):
    N1=2**(k+1)
    index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
    plt.errorbar(index_Masse,c0_fit[k],yerr=c0_fit_sd[k],fmt=".-",label=("N1="+str(N1)))
plt.xlabel("M_bar")
plt.ylabel("cS(0)")
plt.title("Normalisation constant, N="+str(N0)+"xN1")
plt.legend() 
plt.show()

#%% Bestimme m_fit für verschiedene Bereiche von t, ignoriere die Kovarianz und mache ein gewichtetes Mittel via AIC

index = np.arange(int(N0/2)+1)
tmin_tot = 1

for i in range(n_Gitter):
    for j in range(n_masse):
        #Bestimme das grösste erlaubte tmax
        if len(np.where(propagator[i,j]==0)[0])==0:
            tmax_tot = int(N0/2)
        else:
            tmax_tot = np.where(propagator[i,j]==0)[0][0]-1

        test = 0
        while test==0:
            if math.isnan(propagator_sd[i,j,tmax_tot]):
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
                    if math.isnan(m_eff_fit[i,j]):
                        a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[i,j,tmin:(tmax+1)],\
                                            p0=[0.0,1.0],\
                                            sigma=propagator_sd[i,j,tmin:(tmax+1)],absolute_sigma=True,\
                                            bounds=(0.0, np.inf))
                    else:
                        a_fit,cov_fit=curve_fit(cosh,index[tmin:(tmax+1)],propagator[i,j,tmin:(tmax+1)],\
                                            p0=[m_eff_fit[i,j],1.0],\
                                            sigma=propagator_sd[i,j,tmin:(tmax+1)],absolute_sigma=True,\
                                            bounds=(0.0, np.inf))
                    
                    
                    theorie = cosh(index[tmin:(tmax+1)],a_fit[0],a_fit[1])
                    logL = sum((propagator[i,j,tmin:(tmax+1)]-theorie)**2/propagator_sd[i,j,tmin:(tmax+1)])/(n-u)
                    
                    AIC[counter] = 2*u + n*np.log(logL) + 2*u*(u+1)/(n-u-1)
                    m_prop[counter] = a_fit[0]
                    c0_prop[counter] = a_fit[1]
                    counter += 1
            
            AIC = AIC - min(AIC)
            weights = np.exp(-0.5*AIC)
            
            m_fit[i,j] = sum(weights*m_prop)/sum(weights)
            m_fit_sd[i,j] = np.sqrt( sum(weights*(m_fit[i,j]-m_prop)**2)/sum(weights) )
            c0_fit[i,j] = sum(weights*c0_prop)/sum(weights)
            c0_fit_sd[i,j] = np.sqrt( sum(weights*(c0_fit[i,j]-c0_prop)**2)/sum(weights) )
        else:
            m_fit[i,j] = float('NaN')
            m_fit_sd[i,j] = float('NaN')
            c0_fit[i,j] = float('NaN')
            c0_fit_sd[i,j] = float('NaN')

# Plotte m_fit, sobald man alle Daten hat

for k in range(n_Gitter):
    N1=2**(k+1)
    index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
    plt.errorbar(index_Masse,m_fit[k],yerr=m_fit_sd[k],fmt=".-",label=("N1="+str(N1)))
plt.xlabel("M_bar")
plt.ylabel("m_fit")
plt.title("Physical Mass vs. Bar Mass, N="+str(N0)+"xN1")
plt.legend() 
plt.show()

for k in range(n_Gitter):
    N1=2**(k+1)
    index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
    plt.errorbar(index_Masse,c0_fit[k],yerr=c0_fit_sd[k],fmt=".-",label=("N1="+str(N1)))
plt.xlabel("M_bar")
plt.ylabel("cS(0)")
plt.title("Normalisation constant, N="+str(N0)+"xN1")
plt.legend() 
plt.show()

#%% Bestimme den kritischen Exponenten

index_Masse = (np.arange(n_masse)+1)*(M_max-M_min)/float(n_masse)+ M_min
a = 12
b = len(index_Masse)
daten = m_fit[n_Gitter-1,a:b]

def polynom(M_bar,M_crit,A,a_crit):
    y = A*(M_bar-M_crit)**a_crit
    return y

a_fit,cov_fit = curve_fit(polynom,index_Masse[a:b],daten,bounds=(0.0, np.inf))
print("M_crit=",a_fit[0],"+-",np.sqrt(cov_fit[0,0]))
print("A=",a_fit[1],"+-",np.sqrt(cov_fit[1,1]))
print("a_crit=",a_fit[2],"+-",np.sqrt(cov_fit[2,2]))
