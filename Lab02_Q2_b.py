#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Lab 2, Question 2. b)

"""

# a)


import numpy as np 
from matplotlib import pyplot as plt
import scipy
from scipy import special


#defining constants

N = 1000 #the amount of slices 

a = 0 #begining of integral

b = np.pi #end of integral

h = (b-a)/N #step size


def J(m, x): #defining the Bessel Functions 
    def f(theta): #defining the function to be integrated 
        return np.cos(m*theta - x*np.sin(theta))
    
    g = f(a) + f(b) #initial sum of integral to build upon 
    for k in range(1, N): 
        if k % 2 == 0: #if k is even 
            g += 2*f(a+k*h)
        else: g += 4*f(a+k*h) #if k is not even 
    intsum = (h/3)*g #final sum of integral 

    return (1/np.pi)*intsum #what the entire J function should return 



x = np.array(np.arange(0, 20.1, 0.1)) #x values from 0 to 20 


J0 = [] #making empty list to fill in J values, this is J where m = 0
for i in x: #for loop to loop through x and then append this J to the J0 list
    J0.append(J(0, i))


J1 = [] #same as above and below but for m = 1 (below is m = 2)
for i in x:
    J1.append(J(1, i))
    

J2 = []
for i in x:
    J2.append(J(2, i))
    

plt.plot(x, J0, label = '$J_0$') #graphing the Bessel function from my code 
plt.plot(x, J1, label = '$J_1$')
plt.plot(x, J2, label = '$J_2$')
plt.xlabel('x')
plt.ylabel('J')
plt.title("Simpson's Integration Method with N = 10: Bessel Functions") 
plt.legend()
plt.show()




J0_scipy = [] #doing the same append to list process as before but now filling list with scipy function output for bessel functions
for i in x:
    J0_scipy.append(scipy.special.jv(0, i))


J1_scipy = []
for i in x:
    J1_scipy.append(scipy.special.jv(1, i))


J2_scipy = []
for i in x:
    J2_scipy.append(scipy.special.jv(2, i))
            

plt.plot(x, J0_scipy, label = '$J_0$') #graphing the scipy generating bessel functions 
plt.plot(x, J1_scipy, label = '$J_1$')
plt.plot(x, J2_scipy, label = '$J_2$')
plt.xlabel('x')
plt.ylabel('J')
plt.title('Scipy Integration Method: Bessel Functions') 
plt.legend()
plt.show()






#b)


def I(r): #defining intensity function of r 
    k = (2*np.pi)/500
    return (J(1, (k*r))/(k*r))**2



x = (np.arange(101)-50)*20 #making x and y lists, 101-50 will ensure 50 multplied by 20 will be provide values from -1000 to 1000, centered and sorts the indexing at 0 issue
y = (np.arange(101)-50)*20



r = np.zeros([len(x), len(y)]) #making 2d array of zeros for r to fill up 
Intensity = np.zeros([len(x), len(y)]) # making intensity array (also 2d) that will be filled up with respect to the generated r array

for i in range(len(x)): #creating a nested for loop for x has index i and y have j, then each element for the two dimensions of r will be filled up by looping through pythagorous theorem for x and y 
    for j in range(len(y)):
        r[i,j] = (((x[i]**2+y[j]**2))**(1/2))
        Intensity[i, j] = I(r[i,j]) # now indexing the intensity array to be at the same index as r to create its values

Intensity[np.isinf(Intensity)] = 0.25 # setting the point in which intensity would be equal to infinity to 0.25 as the hint from the textbook mentioned the limit of J(x) over x as it approaches 0 is 0.5 but it should be squared so it's o.25 






plt.imshow(Intensity, vmax=0.01, extent= (-1000, 1000, -1000, 1000)) # the v max is set to show the same number of rings as the figure in the textbook, the extent is from -1000 to 1000 nanometers for x and y as specified in text (but in the units of micrometers)
plt.colorbar() #setting the side bar to show intensity corresponding to what colour
plt.xlabel('Radius (nm)')
plt.title('Intensity Map') 
plt.legend()
plt.show()

