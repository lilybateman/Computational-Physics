#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 3, Question 1. b)

"""

import numpy as np
from gaussxw import gaussxw #Using code professor provided from text, requires gaussxw.py
from matplotlib import pyplot as plt #for graphing

def P(u10, Ta, th): #defining probability function 
    delta = 4.3 + (0.145*Ta) + (0.00196*((Ta)**2)) #defining standard deviation 
    u2 = 11.2 + (0.365*Ta) + (0.00706*(Ta**2)) + (0.9*(np.log(th))) #defining mean of u 
    
    def f(u): #defining the function under the integral 
        return (np.exp(-((u2-u)**2)/(2*((delta)**2))))
    
    N = 100 #number of slices 
    
    a = 0 #the lower bound of integration

    b = u10 #the upper bound of integration 
    

    x,w = gaussxw(N) #code from text 
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w


    s = 0.0 #starting value of of sum
    for k in range(N):
        s += wp[k]*f(xp[k]) #increasing the sum


    return ((1/((2*np.pi)**(1/2))*delta)*s) #returning the value outside the integral multplied by the sum of the integration


Ta = np.arange(-40, 15) #making the values of Ta took the tempature range from -40 to 10 celsius 


P1 = [] #creating empty lists of probabilities for different parameters
P2 = []
P3 = []
P4 = []
P5 = []
P6 = []
P7 = []
P8 = []
P9 = []
for i in range(len(Ta)): #creating a for loop to fill the probability arrays
    P1.append(P(6, Ta[i], 24)) 
    P2.append(P(8, Ta[i], 48))
    P3.append(P(10, Ta[i], 72))
    P4.append(P(6, Ta[i], 48))
    P5.append(P(8, Ta[i], 72))
    P6.append(P(10, Ta[i], 24))
    P7.append(P(6, Ta[i], 72))
    P8.append(P(8, Ta[i], 24))
    P9.append(P(10, Ta[i], 48))
    

    
 
plt.plot(Ta, P1, label = '$u_{10} = 6, t_h = 24$', color = 'r', linestyle = 'solid') #making all the different lines with different parameters
plt.plot(Ta, P2, label = '$u_{10} = 8, t_h = 48$', color = 'g', linestyle = 'dashed')
plt.plot(Ta, P3, label = '$u_{10} = 10, t_h = 72$', color = 'b', linestyle = 'dotted')
plt.plot(Ta, P4, label = '$u_{10} = 6, t_h = 48$', color = 'r', linestyle = 'dashed')
plt.plot(Ta, P5, label = '$u_{10} = 8, t_h = 72$', color = 'g', linestyle = 'dotted')
plt.plot(Ta, P6, label = '$u_{10} = 10, t_h = 24$', color = 'b', linestyle = 'solid')
plt.plot(Ta, P7, label = '$u_{10} = 6, t_h = 72$', color = 'r', linestyle = 'dotted')
plt.plot(Ta, P8, label = '$u_{10} = 8, t_h = 24$', color = 'g', linestyle = 'solid')
plt.plot(Ta, P9, label = '$u_{10} = 10, t_h = 48$', color = 'b', linestyle = 'dashed')

plt.xlabel('Tempature ($C^o$)') #labelling plot
plt.ylabel('Probability of Blowing Snow ')
plt.title('Probability of Blowing Snow for Various $u_{10}$ and $t_h$ Parameters') 
plt.legend()
plt.show()

