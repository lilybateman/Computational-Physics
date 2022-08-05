# LAB 11 QUESTION 1 part b

# import statements
from math import sqrt, pi, log
from random import random, randrange
from numpy import linspace, cos, exp, sin
from pylab import plot, xlabel, ylabel, show, figure, title

# use simulated annealing starting at x = 2, y = 2
# at each step x = x + d, y = y + k
# where d and k are random numbers sampled from a gaussian distribution with mean zero and standard deviation 1
# find the global minimum of the function
# use an exponential cooling scheduele

# define the function we want to find the minimum for
def f(x, y):
    return x**2 - cos(4 * pi * x) + ((y-1)**2)

# define gaussian distribution
sigma = 1
def gaussian():
    r = sqrt(-2 * sigma * sigma * log(1 - random()))
    theta = 2*pi*random()
    n = r * cos(theta)
    m = r * sin(theta)  # check
    return n, m

# markov chain monte carlo simulation but lower the temperature: annealing

Tmax = 1  # maxmim temperature: pick such that beta*dfcn << 1
beta = 1 / (Tmax)  # constant
Tmin = 0.00001  # minimum temperature: we want to cool untill we hit
tou = 1000

# define x values
x_vals = []
x = 2
x_vals.append(x)

# define y values
y_vals = []
y = 2
y_vals.append(y)

# define time
time = 0
timevals = []
timevals.append(time)

# start at max temperature
T = Tmax

# while the temperature hasn't reached the minimum temperature yet
while T > Tmin:
    time += 1
    T = Tmax * exp(-time / tou)
    
    # choose a move uniformly at random from a gaussian distribution
    deltax, deltay = gaussian()

    fcn_previous = f(x, y)
    x = x + deltax
    y = y + deltay
    fcn = f(x, y)

    # find the difference in the function value (or 'energy')
    dfcn = fcn_previous - fcn
    beta = 1 / (T)

    # decide wether or not to accept the move
    if random() < (exp(-beta * dfcn)):
        # reject the new move and restore values to their previoud values
        print('reject move')
        fcn = fcn_previous
        x = x - deltax
        y = y - deltay

    elif fcn < fcn_previous:
        # accept the better move
        print('accept move')

    else:
        # reject the new move and restore values to their previoud values
        print('reject move')
        fcn = fcn_previous
        x = x - deltax
        y = y - deltay

    # add x, y, and time values to their respective lists
    x_vals.append(x)
    y_vals.append(y)
    timevals.append(time)

# plot the x and y values over time as they converge to the minimum
figure(1)
plot(timevals, x_vals, '.')
plot(timevals, y_vals, '.')
xlabel('Time')
ylabel('X and Y')
title('Path To Find The Global Minimum')
show()

