# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 10:14:36 2022

@author: samue
"""

import numpy as np
import matplotlib.pyplot as plt

# %%

def chi(t,liste):
    
    summe1 = 0
    summe2 = 0
    summe3 = 0
    
    summe2_square = 0
    summe3_square = 0
    
    n = len(liste)-t
    for i in range(n):
        summe1 += liste[i]*liste[i+t]
        summe2 += liste[i]
        summe3 += liste[i+t]
        
        summe2_square += liste[i]**2
        summe3_square += liste[i+t]**2
    
    o2 = np.sqrt( n*summe2_square - summe2**2)
    o3 = np.sqrt( n*summe3_square - summe3**2)
    
    return( (n*summe1 - summe2*summe3) / (o2*o3) )

tau_list = np.zeros([20,4])

for i in range(20):
    
    M = float((i+1))/10.0
    liste = np.loadtxt("7Vertex_TopologyList_M="+f'{M:.2f}'+",N=2x2.dat")
    
    for j in range(4):
        
        tau = 1.0
        t = 1
        chi0 = chi(0,liste[:,j])
        
        while 6.0*tau > float(t):
            tau += chi(t,liste[:,j])/chi0
            t += 1
        
        tau_list[i,j] = tau

index = (np.arange(20)+1)/10
plt.plot(index,tau_list[:,0],label="Z_00")
plt.plot(index,tau_list[:,1],label="Z_10")
plt.plot(index,tau_list[:,2],label="Z_01")
plt.plot(index,tau_list[:,3],label="Z_11")
plt.ylim([0.0,3.2])
plt.xlabel("M")
plt.ylabel("tau_int")
plt.title("Autokorrelation Topologische Klassen (Pearson Formel)")
plt.legend(loc=0)
plt.show()

# %%

def chi2(t,liste):
    
    summe1 = 0
    summe2 = 0
    summe3 = 0
    
    n = len(liste)-t
    for i in range(n):
        summe1 += liste[i]*liste[i+t]
        summe2 += liste[i]
        summe3 += liste[i+t]

    return((summe1/n) - (summe2/n)*(summe3/n) )

tau_list = np.zeros([20,4])

for i in range(20):
    
    M = float((i+1))/10.0
    liste = np.loadtxt("7Vertex_TopologyList_M="+f'{M:.2f}'+",N=2x2.dat")
    
    for j in range(4):
        
        tau = 1.0
        t = 1
        chi0 = chi2(0,liste[:,j])
        
        while 6.0*tau > float(t):
            tau += chi2(t,liste[:,j])/chi0
            t += 1
        
        tau_list[i,j] = tau
        
index = (np.arange(20)+1)/10
plt.plot(index,tau_list[:,0],label="Z_00")
plt.plot(index,tau_list[:,1],label="Z_10")
plt.plot(index,tau_list[:,2],label="Z_01")
plt.plot(index,tau_list[:,3],label="Z_11")
plt.ylim([0.0,3.2])
plt.xlabel("M")
plt.ylabel("tau_int")
plt.title("Autokorrelation Topologische Klassen (Barkema Formel)")
plt.legend(loc=0)
plt.show()
