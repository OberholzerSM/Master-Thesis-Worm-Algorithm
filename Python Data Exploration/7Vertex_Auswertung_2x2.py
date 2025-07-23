# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:21:36 2022

@author: samue
"""
import numpy as np
import matplotlib.pyplot as plt
n_masse = 20
M_min = 2.0
M_max = 4.0

#%% Histogram topologische Klassen

names=["Z_00","Z_10","Z_01","Z_11"]

theory = np.zeros((n_masse,4))
number = np.zeros((n_masse,4))
fehler = np.zeros((n_masse,4))

for i in range(n_masse):
    M = (i+1)*(M_max-M_min)/n_masse + M_min
    
    weights = np.array([M**8+1.0+1.0+0.0625+0.0625+0.0625+0.0625,
              M**4+M**4+0.0625+0.0625+0.0625+0.0625,
              M**4+M**4+0.0625+0.0625+0.0625+0.0625,
              0.0625+0.0625+0.0625+0.0625])
    theory[i] = weights / ( M**8 + 4*M**4 + 3)

    configsnumber = np.loadtxt("7Vertex_topology_number_M="+f'{M:.2f}'+",N=2x2.dat")
    n_data = int(sum(configsnumber))
    number[i] = configsnumber / n_data

    sigma = np.loadtxt("7Vertex_topology_sd_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler[i] = sigma


index = 2*(np.arange(n_masse)+1)/n_masse
bar_width = 0.8/n_masse

for i in range(4):
    plt.bar(index-0.5*bar_width,theory[:,i],bar_width,label="Theory",edgecolor="black")
    plt.bar(index+0.5*bar_width,number[:,i],bar_width,label="MCMC",color="Green",edgecolor="black")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Frequency")
    plt.title(names[i]+", n_data="+str(n_data))
    plt.legend(loc=0)
    plt.show()
 
    plt.plot(index,np.zeros(n_masse),color="grey")
    plt.scatter(index,number[:,i]-theory[:,i],color = 'red',label="Deviation from Theory")
    plt.errorbar(index,number[:,i]-theory[:,i],yerr=fehler[:,i],fmt="none",ecolor = 'orange',label="Statistical Error")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Deviation from Theory")
    plt.title(names[i]+", n_data="+str(n_data))
    plt.show()

#%% Kritischer Punkt

Z_theory = np.zeros(n_masse)
Z_number = np.zeros(n_masse)
Z_fehler = np.zeros(n_masse)

for i in range(n_masse):
    Z_theory[i] = theory[i,0]-theory[i,1]-theory[i,2]-theory[i,3]
    Z_number[i] = number[i,0]-number[i,1]-number[i,2]-number[i,3]
    Z_fehler[i] = np.sqrt( sum(fehler[i,:]**2) )

index = (np.arange(n_masse)+1)*(M_max-M_min)/n_masse + M_min

plt.plot(index,np.zeros(n_masse),color="grey")
plt.plot(index,Z_theory,label="Theory")
plt.errorbar(index,Z_number,yerr=Z_fehler,ls="None",color="Green")
plt.scatter(index,Z_number,s=30,color="Green",label="MCMC")
plt.xticks(index[::2])
plt.xlabel("Mass")
plt.ylabel("Z_00-Z_10-Z_01-Z_11")
plt.title("Critical Point, n_data="+str(n_data))
plt.legend(loc=0)
plt.show()

plt.plot(index,np.zeros(n_masse),color="grey")
plt.scatter(index,Z_number-Z_theory,color = 'red',label="Deviation from Theory")
plt.errorbar(index,Z_number-Z_theory,yerr=Z_fehler,fmt="none",ecolor = 'orange',label="Statistical Error")
plt.xticks(index[::2])
plt.xlabel("Mass")
plt.ylabel("Deviation from Theory")
plt.title("Critical Point Error, n_data="+str(n_data))
plt.show()

#%% Propagator Histogram (x-Achse: Masse)

names=["dt=-1","dt=0","dt=+1"]

theory = np.zeros((n_masse,3))
number = np.zeros((n_masse,3))
fehler = np.zeros((n_masse,3))

for i in range(n_masse):
    M = (i+1)*(M_max-M_min)/n_masse + M_min
    
    weights = np.array([4*M**4 + 8*M**2 + 8,
                        8*M**6 + 8*M**4 + 16*M**2 + 16,
                        4*M**4 + 8*M**2 + 8])
    theory[i] = weights / sum(weights)

    configsnumber = np.loadtxt("7Vertex_dt_number_M="+f'{M:.2f}'+",N=2x2.dat")
    n_data = int(sum(configsnumber))
    number[i] = configsnumber / n_data

    sigma = np.loadtxt("7Vertex_dt_sd_M="+f'{M:.2f}'+",N=2x2.dat")
    fehler[i] = sigma


index = (np.arange(n_masse)+1)*(M_max-M_min)/n_masse + M_min
bar_width = 0.8/n_masse

for i in range(3):
    plt.bar(index-0.5*bar_width,theory[:,i],bar_width,label="Theory",edgecolor="black")
    plt.bar(index+0.5*bar_width,number[:,i],bar_width,label="MCMC",color="Green",edgecolor="black")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Frequency")
    plt.title(names[i]+", n_data="+str(n_data))
    plt.legend(loc=0)
    plt.show()
 
    plt.plot(index,np.zeros(n_masse),color="grey")
    plt.scatter(index,number[:,i]-theory[:,i],color = 'red',label="Deviation from Theory")
    plt.errorbar(index,number[:,i]-theory[:,i],yerr=fehler[:,i],fmt="none",ecolor = 'orange',label="Statistical Error")
    plt.xticks(index[::2])
    plt.xlabel("Mass")
    plt.ylabel("Deviation from Theory")
    plt.title(names[i]+", n_data="+str(n_data))
    plt.show()