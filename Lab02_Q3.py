#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 2 Q 3 

"""
import numpy as np 
import numpy as np 
from matplotlib import pyplot as plt
import scipy
from scipy import special

#all in nanometers

#b)
def q(u):
    alpha = np.pi/(2e-5)
    return (np.sin(alpha*u)**2)


#c)
N = 8

w = 0.1 #unit is meters 
f = 1 #unit is meters 
l = 5e-7 #unit is meters (this is lambda)

a = -w/2 #unit is meters 
b = w/2 #unit is meters 
h = (b-a)/N #unit is meters 


def I(x): #defining main function 
    
    def f(u):
        return (q(u))**(1/2)*np.exp(((1j*2*np.pi*x*u)/(l*f))) #defining function under the integral 

    g = f(a) + f(b) #as by defintion of simpson's method, the beginning of the integral value
    for k in range(1, N):
        if k % 2 == 0: #making it so if the index k is at is an even integer...
            g += 2*f(a+k*h) #then add two multplied by the function 
        else: g += 4*f(a+k*h) #but if the index is an odd integer, multiply it by four 


    intsum = (h/3)*g #as by defintion of simpson's method
    
    return np.abs((intsum))**2



x = (np.arange(101)-50)*20


I_arr = np.zeros([len(x)]) #making zero array to fill with value of I(r)
for i in range(len(x)): #for loop to fill up array with values
    I_arr[i] = I(x[i])

#d)
plt.imshow(I_arr) #plotting
plt.xlabel('Position on the Screen (m)')
plt.title('Diffraction Pattern') 
plt.show()



# e)

# i)

def q2(u):
    alpha = np.pi/(2e-5) # this 2e-5 value is what will get changed 1e-5, 2e-5
    return (np.sin(alpha*u)**2*np.sin((1/2)*alpha*u)**2) 


def I2(x):
    
    def f2(u):
        (q2(u))**(1/2)*np.exp(((1j*2*np.pi*x*u)/l*f))

    g = f2(a) + f2(b) #as by defintion of simpson's method, the beginning of the integral value
    for k in range(1, N):
        if k % 2 == 0: #making it so if the index k is at is an even integer...
            g += 2*f(a+k*h) #then add two multplied by the function 
        else: g += 4*f(a+k*h) #but if the index is an odd integer, multiply it by four 


    intsum = (h/3)*g #as by defintion of simpson's method
    
    return (intsum)**2