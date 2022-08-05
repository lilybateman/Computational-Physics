#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lab 7, Question 1.

"""

# PSEUDOCODE

# apply the adaptive step size method to RK4

# perform two steps of size h (obtaining x1)
# go back to t0 and perform one step of size 2*h (obtaining x2)
# use these to calculate p
# use the Euclidean error
#p = h*delta/((ex**2 + ey**2)**0.5)
#ex = 1/30 * (x1-x2) ($)
#ey = 1/30 * (y1-y2) ($)

# if p > 1 then we exceeed desired accuracy 
    # keep results and move onto time t+2*h but use h prime obtained from p to find new step size 
# elif p < 1 then we dont meet desired accuracy
    # dont keep result, instead repeat current step with new h prime

# only use x1, not x2 for calculating the final solution (x2 is only for updating the step size)

#PSEUDOCODE FINISHED

# import statements
from math import sin
from numpy import arange, array
from pylab import plot, xlabel, ylabel, show, xlim, figure, title
from time import time

# use RK4 to solve the equations

# want to solve dx/dt = a
# da/dt = -G * M * x / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))
# dx/dt = a
# dy/dt = b
# db/dt = -G * M * y / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))
# from t = 0 to t = 10, M = 10, G = 1, L = 2, IC : (x,y) = (1,0) @ t=0

# define constants
G = 1
L = 2
M = 10

# define the four fisrt order ODEs
def f(r, t):
    a = r[0]  # update x component velocity
    b = r[1]  # update y component velocity
    x = r[2]  # update x position
    y = r[3]  # update y position
    fx = a  # dx/dt = a
    fy = b  # dy/dt = b
    fb = -G * M * y / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))  # db/dt
    fa = -G * M * x / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))  # da/dt
    return array([fa, fb, fx, fy], float)


# TIME THE FIRST METHOD
start_firstmethod = time()

# define the interval
start = 0.0
stop = 10.0
N = 10000
h = (stop - start) / N

# initialize lists
tpoints = arange(start, stop, h)
apoints = []  # this is a list of x component velcoities at each time step
bpoints = []  # this is a list of y component velcoities at each time step
xpoints = []  # this is a list of x positions at each time step
ypoints = []  # this is a list of y positions at each time step

# this array starts with the initial x component velocity, y component velocity, x position, yposition, in that order and updates it at each time iteration
r = array([0.0, 1.0, 1.0, 0.0], float)

# use RK4 to solve the ODEs
for t in tpoints:
    
    apoints.append(r[0])  # update list of x component velcoities at each time step
    bpoints.append(r[1])  # update list of y component velcoities at each time step
    xpoints.append(r[2])  # update list of x positions at each time step
    ypoints.append(r[3])  # update list of y positions at each time step

    # RK4 equations
    k1 = h * f(r, t)
    k2 = h * f(r + (0.5 * k1), t + (0.5 * h))
    k3 = h * f(r + (0.5 * k2), t + (0.5 * h))
    k4 = h * f(r + k3, t + h)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
# END
end_firstmethod = time()

difference_firstmethod = end_firstmethod - start_firstmethod

# NEW ########################################################################
# TIME THE SECOND METHOD
start_secondmethod = time()

r = (0,1,1,0)
delta = 10e-6
times = []
h_list = []
t = 0
h = 0.01
plotx = [1] # initialize list of x points using two steps
ploty = [0] # initialize list of y points using two steps

while t <=10: # 11.98 to cover the same distance, check the x,y,v values between this and other method, maybe difference in velcoity is just due to accuracy
    a = r[0]
    b = r[1]
    x = r[2]
    y = r[3]
    # do two steps of size h
    r1 = array([a,b,x,y], float)
    for i in range(2):
        k1 = h * f(r1, t)
        k2 = h * f(r1 + (0.5 * k1), t + (0.5 * h))
        k3 = h * f(r1 + (0.5 * k2), t + (0.5 * h))
        k4 = h * f(r1 + k3, t + h)
        r1 += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        
    x1 = r1[2]
    y1 = r1[3]
        
    # start again and do one step with size 2h
    r2 = array([a,b,x,y], float)
    h2 = 2*h
    
    k1 = h2 * f(r2, t)
    k2 = h2 * f(r2 + (0.5 * k1), t + (0.5 * h2))
    k3 = h2 * f(r2 + (0.5 * k2), t + (0.5 * h2))
    k4 = h2 * f(r2 + k3, t + h2)
    r2 += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
    x2 = r2[2]
    y2 = r2[3]
    
    # Calculate rho
    ex = (1/30) * (x1-x2) 
    ey = (1/30) * (y1-y2)
    p = h*delta/((ex**2 + ey**2)**0.5)
    
    # Calculate new h value
    h_prime = h*((p)**(1/4))
    h = h_prime
    
    # if statements
    
    if p > 1:
        r = r1
        t +=2*h
        times.append(t)
        h_list.append(h)
        
        plotx.append(x1)
        ploty.append(y1)
        # iteration complete
        
    else:
        r = (a,b,x,y)
        # redo
        
        
# END
end_secondmethod = time()

difference_secondmethod = end_secondmethod - start_secondmethod



# plot the orbit of the ball bearing using the two different methods
figure(1)
plot(xpoints, ypoints)
plot(plotx, ploty, 'k.')
xlabel("x")
ylabel("y")
title('The orbit of the ball bearing')
show()

print('difference_firstmethod = ', difference_firstmethod)
print('difference_secondmethod = ', difference_secondmethod)

# FOR PART C

# create a plot of the size of the time step as a function of time 

# time values of adaptive step size method are in a list called times
dtpoints = array(times[1:]) - array(times[:-1])

figure(2)
plot(times[:-1], dtpoints)
xlabel('time')
ylabel('size of time step')
title('size of the time step verses time')
show()














