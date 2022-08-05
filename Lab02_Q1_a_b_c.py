#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 14:46:53 2020

@author: user
"""

import matplotlib.pyplot as plt
import numpy as np


#Let f(x) be our function
def f(x):
    return np.exp(-(x**2))

#x0 is the actual value of the derivative at x=0.5 (f'(0.5))
x0=0.778800783

###Forward difference scheme

##First we need to calculate the derivative f'(0.5) for different step sizes h.
##Let the values of h be in the list h.
h = np.array([10**-16,10**-15,10**-14,10**-13,10**-12,10**-11,10**-10,0.000000001, 0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1, 1])

    
#Let's calculate the derivative for each step size in order to later calculate the error in each derivative.
##To calculate the derivative, we can use the forward difference scheme
##We set the derivative as df   
df = np.zeros([len(h)])
##We apply f'(0.5)=((f(0.5+h)-f(0.5))/h) for all the 17 values of h 
for i in range(len(h)):
    df[i] = ((f(0.5+h[i]) - f(0.5))/h[i])
    
print('The values for the derivative for their respective h values are')
print(df)
    
    
##We then use the values of the derivative and the step size to find the error in each derivative.
err = np.zeros([len(df)])
##To obtain the error, we subtract the 'actual' value of the derivative from the calculated derivative at x=0.5
for i in range(len(df)):
    err[i] = np.abs(df[i]+x0)
     
print('The absolute value of the error values for their respective derivative values are') 
print(err)  

#We can graph the log of errors vs the step size    
plt.loglog(h, err,label='Forward difference scheme')
plt.xlabel('Step size (h)')
plt.ylabel('Log of Errors')
plt.title('The log of Error vs step size')
plt.legend()
plt.show()
