#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Lab 2, Question 2. a)

"""

# a) i 

import numpy as np 
from numpy import e 
import scipy
from scipy import special
import time




N = 8 #number of slices, also this N will be changed throughout the code for different sections, when that occurs there will be comment of the new N

a = 0 #the lower bound of integration

b = 4 #the upper bound of integration 

h = (b-a)/N #the step size 


def G(x): #defining the first part of the Dawson function, the part before the integral 
    return e**(-x**2)

TimeStart_trap = time.time() #beginning timer to time integration process
def f(t): #defining the function under the integral in the dawson function 
    return e**(t**2)



#trapzoid
    

s = 0.5*f(a) + 0.5*f(b) #integration sum beginning point from the definition of trapezoid integration method
for k in range(1, N): #creating a for loop for the integration from 1 (because index starts at 0th element) to N, the number of slices
    s += f(a+k*h) #this is what will be added to the integration sum each iteration 
    
intsum_trap = h*s #this is the integration sum multiplied by the step as in the definition of trapezoidal method


D_trap = G(b)*intsum_trap #this is to get final value for the Dawson function at x = 4, by multplying the integration sum by the function out front of it

TimeEnd_trap = time.time() #ending the timer for the integration process

#print(D_trap)



#simpsons


TimeStart_simp = time.time()

g = f(a) + f(b) #as by defintion of simpson's method, the beginning of the integral value
for k in range(1, N):
    if k % 2 == 0: #making it so if the index k is at is an even integer...
        g += 2*f(a+k*h) #then add two multplied by the function 
    else: g += 4*f(a+k*h) #but if the index is an odd integer, multiply it by four 

intsum_simp = (h/3)*g #as by defintion of simpson's method


D_simp = G(b)*intsum_simp

TimeEnd_simp = time.time()

#print(D_simp)





#scipy special function 


TimeStart_scipy = time.time()

D_scipy = scipy.special.dawsn(b) # built in scipy function of the dawson function at x=4 which is b

TimeEnd_scipy = time.time()

#print(D_scipy)




#ii) 

#For this question N must be changed to 37144 for trapezoid and 270 for Simpson's 



#trapezoid


def df(x): #defining the function of the derivative of the dawson function 
    return 1 - 2*x*scipy.special.dawsn(x) #using the scipy dawson function as the true value

trap_err = (h**2/12)*(df(a) - df(b)) #defining the error for the trapazoid method as per the definition

#print(trap_err) 



#simpson


def dddf(x): # defining function of the third derivative of the dawson function 
    return x*D_simp+(1-2*x*scipy.special.dawsn(x)*2)*x

simp_err = (h**4/180)*(dddf(a) - dddf(b)) #defining the error for the simpson's method as per the definition

#print(simp_err) 



#print(TimeEnd_trap - TimeStart_trap) #printing out the differences between the start times and the end times of the various methods to get the total time they took 
#print(TimeEnd_simp - TimeStart_simp)
#print(TimeEnd_scipy - TimeStart_scipy)




#iii)
    
#For this question, N from the beginning is set to N = 31
 



#trapezoid


D_trapN = D_trap  #making the value of the dawson function for one N equal to the value calculated at the beginning for simplicities sake and for clarity in the code

s_2N = 0.5*f(a) + 0.5*f(b) #starting to define the same integration process as above, all comments are the same notes except now all N are multplied by two 
for k in range(1, (2*N)):
    h_2N = (b-a)/(2*N)
    s_2N += f(a+k*h_2N)
    
intsum_trap2N = h_2N*s_2N 

D_trap2N = G(b)*intsum_trap2N


E2_trap = (1/3)*(np.abs(D_trap2N - D_trapN)) #practical error estimation for trapazoidal method

#print(E2_trap)




#simpson


D_simpN = D_simp

g_2N = f(a) + f(b) 
for k in range(1, (2*N)):
    h_2N = (b-a)/(2*N)
    if k % 2 == 0:
        g_2N += 2*f(a+k*h_2N)
    else: g_2N += 4*f(a+k*h_2N)

intsum_simp2N = (h_2N/3)*g_2N

D_simp2N = G(b)*intsum_simp2N
    

E2_simp = (1/15)*(np.abs(D_simp2N - D_simpN)) #practical error estimation for trapazoidal method


#print(E2_simp)






