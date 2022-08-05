#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 2, Question 1 

"""



import matplotlib.pyplot as plt
import numpy as np


#Let f(x) be our function
def f(x):
    return np.exp(-(x**2))

#x0 is the actual value of the derivative at x=0.5 (f'(0.5))
x0=0.778800783

###Forward differnce scheme

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
    
    
##We then use the values of the derivative and the step size to find the error in each derivative.
err = np.zeros([len(df)])
##To obtain the error, we subtract the 'actual' value of the derivative from the calculated derivative at x=0.5
for i in range(len(df)):
    err[i] = np.abs(df[i]+x0)
    

###Central forward scheme

#Let's calculate the derivative for each step size in order to later calculate the error in each derivative.
##To calculate the derivative, we can use the forward difference scheme
##We set the derivative as df2
df2 = np.zeros([len(h)])
##We apply f'(0.5)=((f(0.5+h)-f(0.5))/h) for all the 17 values of h 
for i in range(len(h)):
    df2[i] = (f(0.5+(h[i]/2)) - f(0.5 -(h[i]/2)))/(h[i])
    

##We then use the values of the derivative and the step size to find the error in each derivative.
err2 = np.zeros([len(df2)])
##To obtain the error, we subtract the 'actual' value of the derivative from the calculated derivative at x=0.5
for i in range(len(df2)):
    err2[i] = np.abs(df2[i]+x0)
    

#Let's plot both curves on the same graph
#We'll plot log of errors vs the step size
plt.loglog(h, err,label='Forward difference scheme')
plt.loglog(h, err2,label='Central difference scheme')
plt.xlabel('Step size (h)')
plt.ylabel('Log of Errors')
plt.title('The log of Error vs step size')
plt.legend()
plt.show()


