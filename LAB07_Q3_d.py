# -*- coding: utf-8 -*-
"""
Lab 7, Question 3. d)

"""

from numpy import array,arange, exp
from matplotlib import pyplot as plt
import numpy as np
from pylab import figure, xlabel, ylabel, title, show, plot

# Constants
m = 9.1093837015e-31 #mass of electron
hbar = 1.0545718176461565e-34 #Planck's constant over 2*pi
e = 1.602176634e-19 # Electron charge
L = 5.2918e-11 # Bohr radius
N = 1000
a = 4.9999999999999995e-11 #m
h = 0.002*a
E0 = 13.6 #eV

r0 = h
S0 = 1
l=0 #this can be changed for various l
n=1 #this can be changed for various n 


# Potential function
def V(r):
    e0 = 8.8541878128e-12
   
    return -(e**2)*(1/4*np.pi*e0*r) # -(e**2)*(1/4*np.pi*e0*r)

def f(j,r,E):

    R = j[0]
    S = j[1]
    fR = S/(r**2) # dR/dr = S/(r**2)
    fS = l*(l+1)*R + ((2*m*(r**2)/(hbar**2)) * (V(r)-E) * R) # dS/dr = l*(l+1)*R + ((2*m*(r**2)/hbar) * (V(r)-E) * R)
    return array([fR,fS],float)


# Calculate the wavefunction for a particular energy
def solve(E):
    R = 0.0
    S = 1.0
    j = array([R,S],float)

    j_list = [] #create empty list to load wave function values into 

    for r in arange(r0,20*a,h): 
        k1 = h*f(j,r,E)
        k2 = h*f(j+0.5*k1,r+0.5*h,E)
        k3 = h*f(j+0.5*k2,r+0.5*h,E)
        k4 = h*f(j+k3,r+h,E)
        j += (k1+2*k2+2*k3+k4)/6
        j_list.append(j[0]) #load values into list for wave function
    return j[0], j_list


# Main program to find the energy using the secant method

E1 = -15*e/n**2
E2 = -13*e/n**2
R2 = solve(E1)[0]

target = e/1000 # this is the target accuracy 
while abs(E1-E2)>target:
    R1,R2 = R2,solve(E2)[0]
    E1,E2 = E2,E2-R2*(E2-E1)/(R2-R1)

print("E =",E2/e,"eV") 


R = solve(E1)[1]

rpoints = arange(r0,20*a,h)


# normalize
# want to calculate the integral of the absolute value squared of R(r)

# calculate the values of R(r) ** 2
absoluteR = []
for i in range(len(R)):
    absoluteR.append( R[i]**2)



def simps_int(a, b, wavef): #using code from lab 04 but instead having it allow the input of a vector

    h = (b-a)/N

    result = wavef[0]+wavef[-1]
    odds = 0.
    for k in range(1, N, 2):
        odds += wavef[k]
    result += 4.*odds  
    evens = 0.
    for k in range(2, N, 2):
        evens += wavef[k]
    result += 2*evens  

    return h*result/3  



# integrate
Normalization_Constant = simps_int(r0, 20*a, absoluteR)

# divide each term in R by the normalizatio constant
final = []
for i in range(len(R)):
    final.append(R[i]/(Normalization_Constant))


plot(rpoints, final)
xlabel('r (m)')
ylabel('R(r): Normalized Wave function')
title('Wave function with l =0, n=1')
show()

#defining the theoretical equations

a0 = hbar**2/(m*e**2)*(1e+9) #scaled this to put it into meters 
def R1(r): #n=1, l=0
    return 2/a0**(3/2)*e**(-r/a0)
    
def R2(r): #n = 2, l = 0
    return 1/(2*2**(1/2)*a0**(3/1))*(2-r/a0)*e**(-r/(2*a0))
    
def R3(r): #n = 2, l = 1
    return 1/(2*6**(1/2)*a0**(3/2))*r/a0*e**(-r/(2*a0))
    


#normalize the first function R1(r) (change index to 2 or 3 to do other functions)
R1_arr = []

for i in range(len(rpoints)):
    R1_arr.append((R1(rpoints[i]))**2)
    
    
constant1 = simps_int(r0, 20*a, R1_arr)

#scale by the constant
final_R1 = []
for i in range(len(R)):
    final_R1.append(R1_arr[i]/(constant1))

#plot
plt.plot(rpoints, final_R1)
xlabel('r (m)')
ylabel('R(r): Normalized Wave function')
title('Theoretical Wave function with l =0, n=1')
plt.show()
