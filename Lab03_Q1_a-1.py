#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Lab 3, Question 1. a)

"""

# 1 a) i 

import numpy as np 
from numpy import e 
import scipy
from scipy import special
from gaussxw import gaussxw #Using code given by professor in text, requires gaussxw.py
from matplotlib import pyplot as plt


N_arr = (2**3, 2**4, 2**5, 2**6, 2**7, 2**8, 2**9, 2**10, 2**11) #making array of N values, in powers of 3 since the first and last value of bounds given were both powers of 2


a = 0 #the lower bound of integration

b = 4 #the upper bound of integration 


def G(x): #defining the first part of the Dawson function, the part before the integral 
    return e**(-x**2)


def f(t): #defining the function under the integral in the dawson function 
    return e**(t**2)


#trapzoid method taken from my work in last assignment
trap = [] #creating empty list to fill

for N in N_arr: #creating a for loop to do multiple slice values
    
    
    h = (b-a)/N #the step size 

    s = 0.5*f(a) + 0.5*f(b) #integration sum beginning point from the definition of trapezoid integration method
    for k in range(1, N): #creating a for loop for the integration from 1 (because index starts at 0th element) to N, the number of slices
        s += f(a+k*h) #this is what will be added to the integration sum each iteration 
    
    intsum_trap = h*s #this is the integration sum multiplied by the step as in the definition of trapezoidal method


    trap.append(G(b)*intsum_trap) #this is to get final value for the Dawson function at x = 4, by multplying the integration sum by the function out front of it






#simpsons method taken from my work in last assignment


simp = []

for N in N_arr: 

    
    h = (b-a)/N #the step size 

    g = f(a) + f(b) #as by defintion of simpson's method, the beginning of the integral value
    for k in range(1, N):
        if k % 2 == 0: #making it so if the index k is at is an even integer...
            g += 2*f(a+k*h) #then add two multplied by the function 
        else: g += 4*f(a+k*h) #but if the index is an odd integer, multiply it by four 

    intsum_simp = (h/3)*g #as by defintion of simpson's method


    D_simp = G(b)*intsum_simp
    
    simp.append(D_simp)
  

#print(simp)




#Gaussian Method using code from text
gauss = [] 
for N in N_arr: 
    
    

    # Calculate the sample points and weights, then map them
    # to the required integration domain
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w

    # Perform the integration
    s = 0.0
    for k in range(N):
        s += wp[k]*f(xp[k])
        

    gauss.append(s*G(b))

#print(gauss)





#scipy special function 


D_scipy = scipy.special.dawsn(b) # built in scipy function of the dawson function at x=4 which is b



#print(D_scipy)



#ii) 

#For this question N must be changed to 37144 for trapezoid and 270 for Simpson's 



#relative error 


trap_err = [] #taking the error as the difference between the integration values and the true value then divided by the true value
for i in range(len(trap)):
    trap_err.append(np.abs((trap[i] - D_scipy)/D_scipy))

simp_err = []
for i in range(len(simp)):
    simp_err.append(np.abs((simp[i] - D_scipy)/D_scipy))


gauss_err = []
for i in range(len(gauss)):
    gauss_err.append(np.abs((gauss[i] - D_scipy)/D_scipy))



# Now doing error values with 2N and N method for each integration method
#trapezoid
    
E2_trap = []

for N in N_arr:   

    s = 0.5*f(a) + 0.5*f(b) 
    for k in range(1, N):
        h = (b-a)/N
        s += f(a+k*h)
    
    intsum_trap = h*s 

    trapN = G(b)*intsum_trap
    
    
    
    s_2N = 0.5*f(a) + 0.5*f(b) #starting to define the same integration process as above, all comments are the same notes except now all N are multplied by two 
    for k in range(1, (2*N)):
        h_2N = (b-a)/(2*N)
        s_2N += f(a+k*h_2N)
    
    intsum_trap2N = h_2N*s_2N 

    trap2N = G(b)*intsum_trap2N
    

    E2_trap.append((1/3)*(np.abs(trap2N - trapN))) #practical error estimation for trapazoidal method
    

#print(E2_trap)




#simpson

E2_simp = []

for N in N_arr: 
    

    g = f(a) + f(b) 
    for k in range(1, N):
        h = (b-a)/N
        if k % 2 == 0:
            g += 2*f(a+k*h)
        else: g += 4*f(a+k*h)

    intsum_simp = (h/3)*g

    simpN = G(b)*intsum_simp
    

    g_2N = f(a) + f(b) 
    for k in range(1, (2*N)):
        h_2N = (b-a)/(2*N)
        if k % 2 == 0:
            g_2N += 2*f(a+k*h_2N)
        else: g_2N += 4*f(a+k*h_2N)

    intsum_simp2N = (h_2N/3)*g_2N

    simp2N = G(b)*intsum_simp2N
    

    E2_simp.append((1/15)*(np.abs(simp2N - simpN))) #practical error estimation for trapazoidal method


#print(E2_simp)



#gaussian 


E2_gauss = []
for N in N_arr: 
    
    x,w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w

    # Perform the integration
    s_N = 0.0
    for k in range(N):
        s_N += wp[k]*f(xp[k])
        
    Gauss_N = s_N*G(b)

    x_2N,w_2N = gaussxw(2*N)
    xp_2N = 0.5*(b-a)*x_2N + 0.5*(b+a)
    wp_2N = 0.5*(b-a)*w_2N

    # Perform the integration
    s_2N = 0.0
    for k in range(N):
        s_2N += wp_2N[k]*f(xp_2N[k])
        
    Gauss_2N = s_2N*G(b)

    E2_gauss.append((np.abs(Gauss_2N - Gauss_N)))


#print(E2_gauss)


#Making all plots
plt.loglog(N_arr, trap_err, label = 'Trapezoid Method')
plt.loglog(N_arr, simp_err, label = "Simpson's Method")
plt.plot(N_arr, gauss_err, label = 'Gaussian Method')
plt.xlabel('Number of Slices N')
plt.ylabel('Relative Error')
plt.title('Integration Method Error Relative to Scipy Function') 
plt.legend()
plt.show()

plt.loglog(N_arr, E2_trap, label = 'Trapezoid Method')
plt.loglog(N_arr, E2_simp, label = "Simpson's Method")
plt.plot(N_arr, E2_gauss, label = 'Gaussian Method')
plt.xlabel('Number of Slices N')
plt.ylabel('Error')
plt.title('Integration Method Practical Error Estimation ') 
plt.legend()
plt.show()


