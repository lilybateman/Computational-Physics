#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 7, Question 2. b)

"""

from numpy import empty,array,arange 
import numpy as np
from matplotlib import pyplot as plt


G = 6.6738e-11 * ( 8760 * 60**2) ** 2 #scaling to have in m^3kg^-1year^-2 * 60 to be in minutes, *60 to be in hours, *8760 because there are that many hours in a year
M = 1.9891 * 10 ** 30  # mass of sun in kg
H = 1/52   # 1 week in years, interval length for Bulirsch-Stoer
delta = 1000  # in meters because want error within 1 km a year
a = 0 #years
b = 5 #years
x0 = 4.4368e12 #m
y0 = 0 #m
vx0 = 0 #m/yr
vy0 =  6.1218e3 * 8760 * 60**2 #scaling to have in m/year 


def f(r): #define the function to return the gravitional force 


    x, y, vx, vy = r[0], r[1], r[2], r[3] #assigning variables to indexes of the r array
    r = (x**2 + y**2)**(1/2) #distance between sun and pluto
    fx = -G * M * x / r ** 3 #x component of force

    fy = -G * M * y / r ** 3 #y componenet of force

    return np.array([vx, vy, fx, fy], float) 


tpoints = arange(a,b,H) #create an array of time to iterate through



r = array([x0, y0, vx0, vy0],float) #defining r with initial values of positions and velocities 

xpoints = [] #creating empty lists to fill the x and y positional data into 
ypoints = []
H = 1 #now assign an H that gives a reasonable run time for the plot
for t in tpoints:

    xpoints.append(r[0]) #adding the x and y positional values to their respective spots in the r array
    ypoints.append(r[1])
    #now using the Bulirsch-stoer method
    n=1 
    r1 = r + 0.5*H*f(r) 
    r2 =r + H*f(r1)

    R1 =empty([1,4],float)
    R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))
    
    error = 2*H*delta #saying what the target accuracy is

    while error > H*delta:
        n += 1 #increase n until the target accuracy defined above is reached 
        h = H/n
        r1 = r + 0.5*h*f(r)
        r2 = r + h*f(r1)
        for i in range(n-1):
            r1 += h*f(r2) 
            r2 += h*f(r1)

        R2 = R1 # Calculating extrapolation estimates
        R1 = empty([n,4] ,float)
        R1[0] = 0.5*(r1 + r2 + 0.5*h*f(r2)) 
        for m in range(1,n):
            epsilon= (R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1)
            R1[m] = R1[m-1] + epsilon 
        error= abs(epsilon[0])

    r = R1[n-1] #setting r to the most accurate esitmate found before going through next iteration 

#plot the positions to give the orbit    
plt.plot(xpoints,ypoints) 
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('The orbit of Pluto')
plt.show()