#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Lab 8, Question 3 

########## part a ###################

# import statements
from numpy import empty, exp, arange, copy, zeros, mean, linspace
from pylab import plot,xlabel,ylabel,show, figure, title

# constants
L = 1.0      # maximum distance in meters
J = 50       # number of grid spaces
a = L / J    # width of each grid space
h = 0.01     # time increment
epsilon = h/1000
g = 9.81     # gravity
H = 0.01     # in meters

# define different times we want to plot the graphs
t1 = 0.0  # in seconds
t2 = 1.0
t3 = 4.0

# define final time iteration
tfinal = t3 + epsilon

# Boundary and initial conditions

# initial u array
u = zeros(J + 1, float)
# array to hold new u values
up = zeros(J + 1, float)
uq = zeros(J + 1, float)

# initial eta array
eta = zeros(J + 1, float)
# array to hold new eta values
etap = zeros(J + 1, float)
etaq = zeros(J + 1, float)

eta_b = zeros(J + 1, float)

# constants used to define initial conditions of eta
A = 0.002
mew = 0.5
sigma = 0.05
xvals = linspace(0, L, J + 1)

constant = mean(A * exp(-(xvals - mew)**2 / sigma**2))
for i in range(len(xvals)): #loading values for eta 
    eta[i] = H + A * exp(-(xvals[i] - mew)**2 / sigma**2) - (constant)

# define initial time
t = 0.0
c = -h / (2 * a)

# iterate through time
while t < tfinal:
    
    # plot at t=0
    if t == 0.0:
        figure(1)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of Surface: $\eta(x, 0)$")
        title('Lax- Wendroff scheme: Shallow Water System at t = 0 s')
        show()

    # iterate through space
    for j in range(0, J + 1):  # 0, ..., 50
        if j == 0:
            etap[0] = eta[0] - c * (u[1] * (eta[1] - eta_b[1]) - u[0] * (eta[0] - eta_b[0]))  # check factor of 2, not sure about this either            up[0] = 0
            up[0] = 0
        # backwards difference scheme for eta at the end
        elif j == (J): #goes to J +1 so J +1 - 1 = J
            etap[j] = eta[j] - c * (u[j] * (eta[j] - eta_b[j]) - u[j - 1] * (eta[j - 1] - eta_b[j - 1]))
            up[j] = 0

        # for interior points

        else:
            one = (1 / 2) * (u[j + 1] + u[j]) + c * (g * eta[j + 1] + u[j + 1]**2 / 2 - g * eta[j] - u[j]**2 / 2) #this is the upper half of u  +
            two = (1 / 2) * (eta[j + 1] + eta[j]) + c * (u[j + 1] * (eta[j + 1] - eta_b[j + 1]) - u[j] * (eta[j] - eta_b[j])) #lower half of u  -
            etab1 = (eta_b[j + 1] + eta_b[j]) / 2

            three = (1 / 2) * (u[j - 1] + u[j]) + c * (g * eta[j] + u[j]**2 / 2 - g * eta[j - 1] - u[j - 1]**2 / 2) #this is the upper half of eta + 
            four = (1 / 2) * (eta[j - 1] + eta[j]) + c * (u[j] * (eta[j] - 0) - u[j - 1] * (eta[j - 1] - 0)) #this is the lower half oe ta -
            etab2 = (eta_b[j - 1] + eta_b[j]) / 2

            up[j] = u[j] - c * (g * two + (one**2 / 2) - g * four - (three**2) / 2) #now loading values into the main arrays for plotting 
            etap[j] = eta[j] - c * (one * (two - etab1) - three * (four - etab2))

    # use new values of u and eta for the next iteration
    eta = copy(etap)
    u = copy(up)

    
    # update time
    t += h

    # Make plots at the given times
    # plot at t=1
    if abs(t-t2)<epsilon:
        figure(2)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of Surface: $\eta(x, 1)$")
        title('Lax- Wendroff scheme: Shallow Water System at t = 1 s')
        show()

    # plot at t=4
    if abs(t - t3) < epsilon:
        figure(3)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of Surface: $\eta(x, 1)$")
        title('Lax- Wendroff scheme: Shallow Water System at t = 4 s')
        show()

########## part b ###################

# import statements
from numpy import empty, exp, arange, copy, zeros, mean, linspace
from pylab import plot, xlabel, ylabel, show, figure, title
import numpy as np

# constants
L = 1.0      # maximum distance in meters
J = 150       # number of grid spaces
a = L / J       # width of each grid space
h = 0.001      # time increment
epsilon = h / 5000
g = 9.81
eta_b = zeros(J + 1, float)

H = 0.01  # m
eta_bs = H - 4 * 10**(-4)  
alpha = 1 / (8 * np.pi)  # m^-1
x0 = 0.5  # m

t1 = 0.0  # seconds
t2 = 1.0
t3 = 2.0
t4 = 4.0

tfinal = t4 + epsilon #the time to iterate up until 

# Boundary and initial conditions
# Create arrays
u = zeros(J + 1, float)
up = zeros(J + 1, float)
uq = zeros(J + 1, float)

eta = zeros(J + 1, float)
etap = zeros(J + 1, float)
etaq = zeros(J + 1, float)


A = 2 * 10**(-4)  # m

sigma = 0.01  # m
xvals = linspace(0, L, J + 1)


constant = mean(A * exp(-(xvals)**2 / sigma**2))

for i in range(len(xvals)): #make a for loop to load values of eta 
    eta[i] = H + A * exp(-(xvals[i])**2 / sigma**2) - (constant)

for i in range(len(xvals)): #now make one for eta_b because it is a function now
    eta_b[i] = (eta_bs / 2) * (1 + np.tanh((xvals[i] - x0) * alpha))

t = 0.0
c = -h/(2*a)

# iterate through time
while t < tfinal:
    
    # plot at t=0
    if t == 0.0:
        figure(4)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of Surface: $\eta(x, 0)$")
        title('Lax- Wendroff scheme: Tsunami at t = 0 s')
        show()

    # iterate through space
    for j in range(0, J + 1):  # 0, ..., 50
        if j == 0:
            etap[0] = eta[0] - c * (u[1] * (eta[1] - eta_b[1]) - u[0] * (eta[0] - eta_b[0]))  # check factor of 2, not sure about this either            up[0] = 0
            up[0] = 0
        # backwards difference scheme for eta at the end
        elif j == (J): #goes to J +1 so J +1 - 1 = J
            etap[j] = eta[j] - c * (u[j] * (eta[j] - eta_b[j]) - u[j - 1] * (eta[j - 1] - eta_b[j - 1]))
            up[j] = 0

        # for interior points

        else:
            one = (1 / 2) * (u[j + 1] + u[j]) + c * (g * eta[j + 1] + u[j + 1]**2 / 2 - g * eta[j] - u[j]**2 / 2) #this is the upper half of u  +
            two = (1 / 2) * (eta[j + 1] + eta[j]) + c * (u[j + 1] * (eta[j + 1] - eta_b[j + 1]) - u[j] * (eta[j] - eta_b[j])) #lower half of u  -
            etab1 = (eta_b[j + 1] + eta_b[j]) / 2

            three = (1 / 2) * (u[j - 1] + u[j]) + c * (g * eta[j] + u[j]**2 / 2 - g * eta[j - 1] - u[j - 1]**2 / 2) #this is the upper half of eta + 
            four = (1 / 2) * (eta[j - 1] + eta[j]) + c * (u[j] * (eta[j] - 0) - u[j - 1] * (eta[j - 1] - 0)) #this is the lower half oe ta -
            etab2 = (eta_b[j - 1] + eta_b[j]) / 2

            up[j] = u[j] - c * (g * two + (one**2 / 2) - g * four - (three**2) / 2) #now loading values into the main arrays for plotting 
            etap[j] = eta[j] - c * (one * (two - etab1) - three * (four - etab2))


    # use new values of u and eta for the next iteration
    eta = copy(etap)
    u = copy(up)

    
    # update time
    t += h

    # Make plots at the given times
    # plot at t=1
    if abs(t-t2)<epsilon:
        figure(5)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of Surface: $\eta(x, 1)$")
        title('Lax- Wendroff scheme: Tsunami at t = 1 s')
        show()

    # plot at t=2
    if abs(t - t3) < epsilon:
        figure(6)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of Surface: $\eta(x, 2)$")
        title('Lax- Wendroff scheme: Tsunami at t = 2 s')
        show()

    # plot at t=4
    if abs(t - t4) < epsilon:
        figure(7)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of Surface: $\eta(x, 4)$")
        title('Lax- Wendroff scheme: Tsunami at t = 4 s')
        show()