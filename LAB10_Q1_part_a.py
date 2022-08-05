'''Q1a : started with code from Brownian_start.py by Nico Grisouard, University of Toronto

Purpose of code:
Model a random walk with one particle completing 5000 steps
'''

# import statements
from random import random, seed, randrange
from numpy import arange, zeros
from pylab import plot, xlabel, ylabel, show, figure, title, xlim, ylim
from matplotlib.pyplot import imshow
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# helper function
def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    direction = randrange(4)  # pick a random number out of 0, 1, 2, or 3

    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y

font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)

# %% main program starts here -------------------------------------------------|

plt.ion()

Lp = 101  # size of domain
Nt = 5000  # number of steps

# arrays to record the trajectory of the particle
positionx = []
positiony = []

# define the starting point which is the middle of the grid
centre_point = (Lp - 1) // 2 
xp = centre_point
yp = centre_point

# main loop
for i in range(Nt):  # for every time step
    xpp, ypp = nextmove(xp, yp)  # randomly pick the next step

    if (0 < xpp < (Lp)) and (0 < ypp < (Lp)):  # make sure its contained within the walls of the grid
        positionx.append(xpp)  # add postition to list recording trajectory information
        positiony.append(ypp)
        xp = xpp  # update for the next step
        yp = ypp

    else:
        while (xpp > Lp) or (ypp > Lp) or (xpp < 0) or (ypp < 0):  # if a particle is on the wall, make it start on a different trajectory, so it doesnt leave the grid
            xpp, ypp = nextmove(xp, yp)

        positionx.append(xpp)  # add postition to list recording trajectory information
        positiony.append(ypp)
        xp = xpp  # update for the next step
        yp = ypp

# plot
figure(1)
title('Random walk for 1 particle')
xlabel('x')
ylabel('y')
xlim([-1, Lp])
ylim([-1, Lp])
plot(positionx, positiony)
show()