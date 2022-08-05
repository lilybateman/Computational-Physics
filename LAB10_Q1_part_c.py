'''Q1c : started with code from DLA_start.py by Nico Grisouard, University of Toronto

Purpose of code:
Animate DLA trajectory for just anchored articles. 
Particles are introduced one after another at the center of the grid and walk 
until a wall is hit or another anchored particle.
'''

# import
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from pylab import clf, plot, xlim, ylim, show, pause, draw
from random import random, seed, randrange
from numpy import arange, zeros
from pylab import plot, xlabel, ylabel, show, figure, title
from matplotlib.pyplot import imshow

# excersise 10.13 part b
# use 100 particles

# helper function which picks direction of the next step
def nextmove(x, y):
    """ randomly choose a direction
    1 = up, 2 = down, 3 = left, 4 = right"""
    direction = randrange(1, 5)  # pick a random number out of 1, 2, 3, or 4

    if direction == 1:  # move up
        y += 1
    elif direction == 2:  # move down
        y -= 1
    elif direction == 3:  # move right
        x += 1
    elif direction == 4:  # move left
        x -= 1
    else:
        print("error: direction isn't 1-4")

    return x, y

font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)

# %% main program starts here ------------------------------------------------|

plt.ion()

Lp = 101  # size of domain
N = 100  # number of particles

# array to represent whether each gridpoint has an anchored particle
# 1 means there is an anchored particle there
# 0 means there is no anchored particle there
anchored = np.zeros((Lp + 1, Lp + 1), dtype=int)

# list to represent x and y positions of anchored cells
anchored_points = [[], []]

# define the starting point which is the middle of the grid, each of the 1oo particles will be introduced here
centre_point = (Lp - 1) // 2

# set up animation of anchored points
plt.figure(2)
plt.title('DLA run for {} particles'.format(N))
plt.plot(centre_point, centre_point, '.r', markersize=10)
plt.xlim([-1, Lp])
plt.ylim([-1, Lp])
plt.xlabel('$x$ []')
plt.ylabel('$y$ []')

# initialize variables for storing final position of anchored particles
freezex = 0
freezey = 0

# initialize lists for storing the full trajectories for all the particles
totalx = []
totaly = []

# for each particle j
for j in range(N):
    
    # arrays to record the trajectory of the particle j
    positionx = []
    positiony = []

    # set first point as the center point
    xp = centre_point
    yp = centre_point

    xpp, ypp = nextmove(xp, yp)  # pick next direction

    # iterate while the particle is in the grid and is not next to an anchored cell, as well as untill the center becomes an anchored point
    while (0 < xp < (Lp)) and (0 < yp < (Lp)) and (anchored[xp + 1][yp] == 0) and (anchored[xp + 1][yp] == 0) and (anchored[xp][yp + 1] == 0) and (anchored[xp][yp - 1] == 0) and (anchored[50][50] == 0):

        positionx.append(xpp)  # add the move to the full trajectory of x positions
        positiony.append(ypp)  # add the move to the full trajectory of y positions

        xp = xpp  # update x position
        yp = ypp  # update y position

        freezex = xp  # store the final x position of the anchored particle
        freezey = yp  # store the final y position of the anchored particle

        xpp, ypp = nextmove(xp, yp)  # pick next direction

    anchored_points[0].append(freezex)  # add to list of anchored cells
    anchored_points[1].append(freezey)
    anchored[freezex][freezey] = 1  # record the precense of an anchored cell at this grid point

    totalx.append(positionx)  # add trajectory of particle j to list of trajectories of all 100 particles
    totaly.append(positiony)

    # animation
    plt.xlim([-1, Lp])
    plt.ylim([-1, Lp])
    plot(positionx, positiony)
    draw()
    pause(0.01)  # pause to allow a smooth animation
        

# plot of the final trajectory
figure(1)
for i in range(len(totalx)):
    plot(totalx[i], totaly[i])
title('DLA run for {} particles'.format(N))
xlabel('x')
ylabel('y')
plt.xlim([-1, Lp])
plt.ylim([-1, Lp])
show()








