# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 11:05:00 2022

@author: samue
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:21:36 2022

@author: samue
"""
import numpy as np
import matplotlib.pyplot as plt

# %%

#Histogram über alle Konfigurationen

for i in range(20):
    M = float((i+1))/10.0
    
    theory = np.array([M**4,1.0,1.0,0.25,0.25,0.25,0.25,
              M**2,M**2,0.25,0.25,0.25,0.25,
              M**2,M**2,0.25,0.25,0.25,0.25,
              0.25,0.25,0.25,0.25])
    theory = theory / sum(theory)

    configsnumber = np.loadtxt("7Vertex_configsnumber_M="+f'{M:.2f}'+",N=2x2.dat")
    n_data = int(sum(configsnumber))
    configsnumber = configsnumber / n_data

    fehler = np.loadtxt("7Vertex_sd_configs_Block(BootstrapBlock)_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler = fehler /n_data
   
    index = np.arange(len(configsnumber))
    bar_width = 0.40

    plt.bar(index ,theory,bar_width,label="Theory",edgecolor="black")
    plt.bar(index + bar_width,configsnumber,bar_width,label="MCMC",color="Green",edgecolor="black")
    plt.xticks(np.array([0,5,10,15,20]) + 0.5*bar_width,np.array([0,5,10,15,20]))
    plt.xlabel("Configuration")
    plt.ylabel("Frequency")
    plt.title("M="+str(M)+", n_data="+str(n_data))
    plt.legend(loc='upper right')
    plt.show()

    plt.plot(index,np.zeros(23),color="grey")
    plt.scatter(index,configsnumber-theory,color = 'red',label="Deviation from Theory")
    plt.errorbar(index,configsnumber-theory,yerr=fehler,fmt="none",ecolor = 'orange',label="Statistical Error")
    plt.xticks(np.array([0,5,10,15,20]) + 0.5*bar_width,np.array([0,5,10,15,20]))
    plt.xlabel("Configuration")
    plt.ylabel("Error")
    plt.title("M="+str(M)+", n_data="+str(n_data))
    plt.show()


#%% Histogram topologische Klassen

names=["Z_00","Z_10","Z_01","Z_11"]

theory = np.zeros((20,4))
number = np.zeros((20,4))
fehler = np.zeros((20,4))

for i in range(20):
    M = float(2*(i+1))/20.0
    
    weights = np.array([M**4+1.0+1.0+0.25+0.25+0.25+0.25,
              M**2+M**2+0.25+0.25+0.25+0.25,
              M**2+M**2+0.25+0.25+0.25+0.25,
              0.25+0.25+0.25+0.25])
    theory[i] = weights / sum(weights)

    configsnumber = np.loadtxt("7Vertex_topology_number_M="+f'{M:.2f}'+",N=2x2.dat")
    n_data = int(sum(configsnumber))
    number[i] = configsnumber / n_data

    sigma = np.loadtxt("7Vertex_sd_topology_Block(BootstrapBlock)_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler[i] = sigma / n_data

index = (np.arange(20)+1)/10
bar_width = 0.04

for i in range(4):
    plt.bar(index-0.5*bar_width,theory[:,i],bar_width,label="Theory",edgecolor="black")
    plt.bar(index+0.5*bar_width,number[:,i],bar_width,label="MCMC",color="Green",edgecolor="black")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Frequency")
    plt.title(names[i]+", n_data="+str(n_data))
    plt.legend(loc=0)
    plt.show()

    plt.plot(index,np.zeros(20),color="grey")
    plt.scatter(index,number[:,i]-theory[:,i],color = 'red',label="Deviation from Theory")
    plt.errorbar(index,number[:,i]-theory[:,i],yerr=fehler[:,i],fmt="none",ecolor = 'orange',label="Statistical Error")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Deviation from Theory")
    plt.title(names[i]+", n_data="+str(n_data))
    plt.show()
    
#%% Kritischer Punkt

Z_theory = np.zeros(20)
Z_number = np.zeros(20)
Z_fehler = np.zeros(20)

for i in range(20):
    Z_theory[i] = theory[i,0]-theory[i,1]-theory[i,2]-theory[i,3]
    Z_number[i] = number[i,0]-number[i,1]-number[i,2]-number[i,3]
    Z_fehler[i] = np.sqrt( sum(fehler[i,:]**2) )

index = (np.arange(20)+1)/10

plt.plot(index,np.zeros(20),color="grey")
plt.plot(index,Z_theory,label="Theory")
plt.errorbar(index,Z_number,yerr=Z_fehler,ls="None",color="Green")
plt.scatter(index,Z_number,s=30,color="Green",label="MCMC")
plt.xticks(index[::2])
plt.xlabel("Mass")
plt.ylabel("Z_00-Z_10-Z_01-Z_11")
plt.title("Critical Point, n_data="+str(n_data))
plt.legend(loc=0)
plt.show()

plt.plot(index,np.zeros(20),color="grey")
plt.scatter(index,Z_number-Z_theory,color = 'red',label="Deviation from Theory")
plt.errorbar(index,Z_number-Z_theory,yerr=Z_fehler,fmt="none",ecolor = 'orange',label="Statistical Error")
plt.xticks(index[::2])
plt.xlabel("Mass")
plt.ylabel("Deviation from Theory")
plt.title("Critical Point Error, n_data="+str(n_data))
plt.show()

#%% Autokorrelation
index = (np.arange(20)+1)/10
tau = np.loadtxt("7Vertex_Autokorrelation_2x2.dat")
tau_sd = np.loadtxt("7Vertex_Autokorrelation_sd_2x2.dat")
plt.plot(index,tau,color="purple")
plt.scatter(index,tau,color="purple")
plt.errorbar(index,tau,yerr=tau_sd,color="purple")
plt.xlabel("Mass")
plt.ylabel("Autocorrelation")
plt.title("Autocorrelation for different Masses")
plt.show()

#%% Naiver Fehler vs. Block Fehler für Topologische Klassen
#Block 1: Bestimme den Fehler für jeden einzelnen Block via Bootstrap und summiere den Fehler vie Gauss auf
#Block 2: Bootstrappe die einzelnen Blöcke selbst und bestimme so den Fehler

names=["Z_00","Z_10","Z_01","Z_11"]

fehler_naiv = np.zeros((20,4))
fehler_Block = np.zeros((20,4))
fehler_Block2 = np.zeros((20,4))

for i in range(20):
    M = float(2*(i+1))/20.0
    configsnumber = np.loadtxt("7Vertex_topology_number_M="+f'{M:.2f}'+",N=2x2.dat")
    n_data = int(sum(configsnumber))
    
    sigma = np.loadtxt("7Vertex_sd_topology_naiv_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler_naiv[i] = sigma / n_data
    
    sigma = np.loadtxt("7Vertex_sd_topology_Block_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler_Block[i] = sigma / n_data
    
    sigma = np.loadtxt("7Vertex_sd_topology_Block(BootstrapBlock)_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler_Block2[i] = sigma / n_data

index = (np.arange(20)+1)/10
for i in range(4):
    plt.scatter(index,fehler_naiv[:,i],color="Red",label="Naiv")
    plt.scatter(index,fehler_Block[:,i],color="Orange",label="Block 1")
    plt.scatter(index,fehler_Block2[:,i],color="Yellow",label="Block 2")
    plt.xlabel("Mass")
    plt.ylabel("Statistical Error")
    plt.title("Naiv Error vs. Block Error for "+names[i])
    plt.legend(loc=0)
    plt.show()