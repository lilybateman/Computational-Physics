" Lab 6, Question 3. c)"


# import statements

from math import sin

import numpy as np

from numpy import arange, array, meshgrid, sqrt, ones, empty

from pylab import plot, xlabel, ylabel, show, xlim, figure, title, ylim

import matplotlib.pyplot as plt


#defining constants, portion given in the lab handout to create the initial conditions
N = 16 #number of particles

Lx = 4.0

Ly = 4.0

dx = Lx / sqrt(N)

dy = Ly/sqrt(N)

x_grid = arange(dx/2, Lx, dx)

y_grid = arange(dy/2, Ly, dy)

xx_grid, yy_grid = meshgrid(x_grid, y_grid)

x_initial = xx_grid.flatten()

y_initial = yy_grid.flatten()



#creating an r array of just zeros to fill, with the right dimensions
r_temp = array((([0, 0, 0, 0]),) * (N + 1), float)

r = r_temp[0:N] 





# creating a for loop to fill the initial conditions into their respective places: in the order a (vx),b (vy),x,y at t=0
for i in range(len(x_initial)):

    r[i][2] = x_initial[i] 

    r[i][3] = y_initial[i]

    r[i][0] = 0 #x component of velocity 

    r[i][1] = 0 #y component of velocity 



# define the interval

start = 0.0 #initial time 

stop = 10.0 #the end time

T = 1000 #this is the number of time steps 

h = (stop - start) / T #this is the change in time (time step) also known as dt 



#create the data set of times to loop through

tpoints = arange(start, stop, h)


#define the function to return accelerations, to be used in the verlet 
def f(r, tpoints):

    alvals = np.array((([0] * 4),) * (N), float) #create an array of zeros that is suited for the N number of particles

    for i in range(len(r)): # create a for loop for a particle i, to go through and assign each variable to its designated index as specified in the initial r array 

        fai = 0 #initial acceleration in the x direction 

        fbi = 0 #initial acceleration in the y direction 

        ai = r[i][0] #velocity x component for particle i 

        bi = r[i][1] #velocity y component for particle i 

        xi = r[i][2] #position x component for particle i 

        yi = r[i][3] #position y component for particle i 

        fxi = ai #velocity x for particle i 

        fyi = bi #velocity y for particle i 

        for j in range(len(r)): #create another for loop for a particle j 

        # for each particle i

            xj = r[j][2] #position x component for particle j 
 
            yj = r[j][3] #position y component for particle j 



            if i != j: #if the ith particle is not the same index as the jth particle



                x = xj - xi #then define the difference between their x positions as the final x for the acceleration

                y = yj - yi #and the difference between their y positions as the final y for acceleration

                fai_t = 2 * ((-6 * 2 * x * (x**2 + y**2)**-7) + 3 * 2 * x * (x**2 + y**2)**-4) #defining the acceleration for the x direction

                fbi_t = 2 * ((-6 * 2 * y * (x**2 + y**2)**-7) + 3 * 2 * y * (x**2 + y**2)**-4) #defining the acceleration for the y direction 

                fai += fai_t #add the new accelerations to the original acceleration of 0 

                fbi += fbi_t

        alvals[i][0] = fai #indexing to return an array with an order of fai, fbi, fxi, fyi

        alvals[i][1] = fbi

        alvals[i][2] = fxi

        alvals[i][3] = fyi



    return alvals #final array returned




v = 0 #the initial velocity 

v0 = v + 0.5 * h * f(r, tpoints) #defining the velocity for the verlet



#empty arrays to index values into for position and velocity in the verlet
xpoints = [[] for i in range((N))]

ypoints = [ [] for i in range((N)) ]

apoints = [ [] for i in range((N)) ]

bpoints = [ [] for i in range((N)) ]


#define a for loop to iterate through all the t values
for t in tpoints:

    for i in range(N): #create another for loop to create a list for each set of x, y, vx, and vy points

        apoints[i].append(r[i][0]) #they are all indexed at their respective place in r as specified at the beginning in the initial r array

        bpoints[i].append(r[i][1])

        xpoints[i].append(r[i][2])

        ypoints[i].append(r[i][3])



        
# below is the actual calculations from the verlet method as described in the lab handout
    r= r + h * v0 

    k= h * f(r, t + h)

    v= v0 + 0.5 * k

    v0= v0 + k

    for i in range(8): # repeating the position 8 times as stated in the lab handout

        r[i][2] = np.mod(r[i][2], Lx) #using indexes for x position 
        r[i][3] = np.mod(r[i][3], Lx) #using indexes for y position 

    # it is known the above code is wrong, cannot figure out how to repeat the position 8 times 
    



# below is plotting for question 3 c 
figure(1)

for i in range(N):

    plot(xpoints[i], ypoints[i], '.')



xlabel("x")

ylabel("y")

xlim(0, Lx) #from the bounds in handout

ylim(0, Ly)

title('The position of 16 particles')

show()



