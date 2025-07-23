# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 08:59:44 2022

@author: samue
"""

import numpy as np
import matplotlib.pyplot as plt

N=2

theory_list = np.zeros([20,23])
fortran_list = np.zeros([20,23])
error_list = np.zeros([20,23])

for i in range(20):
    M = (i+1)*0.10
    theory = np.array([M**8,1.0,1.0,0.0625,0.0625,0.0625,0.0625,
              M**4,M**4,0.0625,0.0625,0.0625,0.0625,
              M**4,M**4,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625])
    theory_list[i] = theory / ( M**8 + 4*M**4 + 3) * 100.0
    
    configsnumber = np.loadtxt("7Vertex_configsnumber_M="+f'{M:.2f}'+",N="+str(N)+".dat")
    ndata = int(sum(configsnumber))
    fortran_list[i] = configsnumber / ndata * 100.0
    
    error = np.loadtxt("7Vertex_sd_M="+f'{M:.2f}'+",N="+str(N)+".dat")
    error_list[i] = error / ndata * 100.0

masses = (np.arange(20)+1)/10

#Autokorrelation
tau = np.loadtxt("7Vertex_Autokorrelation_"+str(N)+"x"+str(N)+".dat")
plt.plot(masses,tau,color="purple")
plt.scatter(masses,tau,color="purple")
plt.xlabel("Mass")
plt.ylabel("Auto Correlation Time")
plt.title("Auto Correlation Time, "+str(N)+"x"+str(N)+" Lattice")
plt.show()

for i in range(23):
    plt.plot(masses,theory_list[:,i],label="Theory")
    plt.scatter(masses,theory_list[:,i])
    plt.plot(masses,fortran_list[:,i],color="green",label="Fortran")
    plt.scatter(masses,fortran_list[:,i],color="green")
    plt.xlabel("Mass")
    plt.ylabel("P[%]")
    plt.title("Probability for the "+str(i+1)+"th configuration")
    plt.legend(loc=0)
    plt.show()

    plt.plot(masses,fortran_list[:,i]-theory_list[:,i],color="red",label="Deviation from Theory")
    plt.scatter(masses,fortran_list[:,i]-theory_list[:,i],color="red")
    plt.plot(masses,error_list[:,i],color="orange",label="Statistical Error")
    plt.scatter(masses,error_list[:,i],color="orange")
    plt.xlabel("Mass")
    plt.ylabel("Absolute Error [%]")
    plt.title("Absolute Errors for the "+str(i+1)+"th configuration")
    plt.legend(loc=0)
    plt.show()
    
    plt.plot(masses,abs(fortran_list[:,i]-theory_list[:,i])/theory_list[:,i]*100.0,color="red",label="Deviation from Theory")
    plt.scatter(masses,abs(fortran_list[:,i]-theory_list[:,i])/theory_list[:,i]*100.0,color="red")
    plt.plot(masses,error_list[:,i]/theory_list[:,i],color="orange",label="Statistical Error")
    plt.scatter(masses,error_list[:,i]/theory_list[:,i],color="orange")
    plt.xlabel("Mass")
    plt.ylabel("Relative Error [%]")
    plt.title("Relative Errors for the "+str(i+1)+"th configuration")
    plt.legend(loc=0)
    plt.show()