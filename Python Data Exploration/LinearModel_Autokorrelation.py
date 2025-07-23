# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 14:59:51 2022

@author: samue
"""

import numpy as np
import math
import matplotlib.pyplot as plt

#%%

liste = np.loadtxt("b_samplings.dat")
n = len(liste)

chi_list = np.zeros(3)

for t in range(3):
    summe1=0.0
    summe2=0.0
    summe3=0.0
    
    for i in range(n-t):
        summe1 += liste[i]*liste[i+t]
        summe2 += liste[i]
        summe3 += liste[i+t]

    sum1 = summe1/float(n-t)
    sum2 = summe2/float(n-t)
    sum3 = summe3/float(n-t)
    
    chi_list[t] =sum1 - sum2*sum3
    
    print(t)
    print(summe1,summe2,summe3)
    print(sum1,sum2,sum3)
    print(sum2*sum3)
    print(chi_list[t])
    print("")

plt.plot(chi_list/chi_list[0])
plt.show()

plt.plot(np.log(chi_list)-np.log(chi_list[0]))
plt.show()