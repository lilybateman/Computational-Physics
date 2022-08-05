#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 4, Question 2, e)

"""

import numpy as np
from matplotlib import pyplot as plt #to make the plots



#defining constants
L = 5*10**(-10) # m
M = 9.1094*10**(-31) #kg
h_bar = 1.0545718 * 10**(-34) #J s
a = 1.602176634*10**(-18) #J

def H(m, n): #defining H function 
    if m != n:
        if (m % 2 == 0 and n % 2 == 0) or (m % 2 != 0 and n % 2 != 0): #putting in all condtions from lab handout
            return 0
        if  (m % 2 != 0 and n % 2 == 0) or (n % 2 != 0 and m % 2 == 0):
            return (-((8*a*m*n)/(np.pi**2*(m**2-n**2)**2)))*6.242e+18 #multplying by constant inorder to get the final answer in electron volts
            
    else:
        return ((1/2)*a+((np.pi**2*h_bar**2*m**2)/(2*M*L**2)))*6.242e+18

#H(m,n) is in terms of electron volts as desired


arr10 = np.zeros((10, 10)) #making a zero 10 by 10 array to load values into

for m in range(1, 10+1): #making nested for loop, using hint from the lab handout for range and index
    for n in range(1, 10+1):
        arr10[m-1, n-1] = H(m,n) #filling up array 
        
Energy_levels10 = np.linalg.eigvalsh(arr10) #these are the eigen values of the 10 by 10 array
#print(Energy_levels10) 


arr100 = np.zeros((100, 100)) #now doing the same process as above but for a 100 by 100 array instead
   
for m in range(1, 100+1):
    for n in range(1, 100+1):
        arr100[m-1, n-1] = H(m,n)


Energy_levels100 = np.linalg.eigvalsh(arr100)
#print(Energy_levels100)


e_values = np.linalg.eigh(arr100)[0] #these are the eigen values 
e_vectors = np.linalg.eigh(arr100)[1] #these are the eigen vectors


transpose_e_vectors = np.transpose(e_vectors) #transposing the eigen vectors 

psi_0 = transpose_e_vectors[: , 0] #slicing to take the eigen vector corresponding to a its eigen value which can be seen by index (ex. 0, 1, 2)
psi_1 = transpose_e_vectors[: , 1]
psi_2 = transpose_e_vectors[: , 2]



L=5 #changing units so plot can be in angstroms 

#defining the functions for each energy state
def psi0(x): #this is the ground state 
    p0 = 0 #start with intial value of 0 
    for i in range(100): #looping through 100 because that is shape of the input data
        p0 += psi_0[i]*np.sin(((i+1)*np.pi*x)/L) #now this will go through accumulating to make sure all the coefficiants are found by indexing through each entry in the eigen vector as part of the wave function equation in the lab handout
    return p0

def psi1(x): #same process as above but for the first excited state
    p1 = 0
    for i in range(100):
        p1 +=  psi_1[i]*np.sin(((i+1)*np.pi*x)/L)
    return p1

def psi2(x): #same process as above but for the second excited state
    p2 = 0
    for i in range(100):
        p2 += psi_2[i]*np.sin(((i+1)*np.pi*x)/L)
    return p2

x = np.arange(0, 5, 0.05) #this is in angstroms, intervals of 0.05 because want 100 values


psi0list = [] #making empty lists to fill
psi1list = []
psi2list = []

for i in range(len(x)): #filling lists with the function of each x value
    psi0list.append(psi0(x[i]))
    psi1list.append(psi1(x[i]))
    psi2list.append(psi2(x[i]))
 
psi0arr = np.array(psi0list) #turning lists into arrays
psi1arr = np.array(psi1list)
psi2arr = np.array(psi2list)


def f(x, phi): #function under integral
    return np.abs(x, phi)**2


N = 100 #number of slices
a = 0 #lower bound of integration
b = L #upper bound of integration
h = (b-a)/N #interval size


#Now going to normalize and get integration constants 
#using simpsons method integration from my code in lab 2 
g = f(a, psi_0) + f(b, psi_0) #as by defintion of simpson's method, the beginning of the integral value
for k in range(1, N):
    if k % 2 == 0: #making it so if the index k is at is an even integer...
        g += 2*f(a+k*h, psi_0) #then add two multplied by the function 
    else: g += 4*f(a+k*h, psi_0) #but if the index is an odd integer, multiply it by four 

A0 = (h/3)*g #as by defintion of simpson's method

g = f(a, psi_1) + f(b, psi_1) 
for k in range(1, N):
    if k % 2 == 0: 
        g += 2*f(a+k*h, psi_1) 
    else: g += 4*f(a+k*h, psi_1) 

A1 = (h/3)*g 

g = f(a, psi_2) + f(b, psi_2) 
for k in range(1, N):
    if k % 2 == 0: 
        g += 2*f(a+k*h, psi_2) 
    else: g += 4*f(a+k*h, psi_2)  

A2 = (h/3)*g 


#plotting all normalized wave functions 
plt.plot(x, 1/(A0)**(1/2)*psi0arr**2, label = "$\psi_{0}$")
plt.plot(x, 1/(A1)**(1/2)*psi1arr**2, label = "$\psi_{1}$")
plt.plot(x, 1/(A2)**(1/2)*psi2arr**2, label = "$\psi_{2}$")
plt.xlabel('Position x (Angstroms)')
plt.ylabel('Probability of Position')
plt.title('The Ground State and First Two Excited States') 
plt.legend()
plt.show()
