" Lab 6, Question 3. b)"


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



# collect information at half steps

v0_apoints = [[] for i in range((N))]

v0_bpoints = [[] for i in range((N))]

v0_xpoints = [[] for i in range((N))]

v0_ypoints = [[] for i in range((N))]


#empty lists for energies 
PE = []

KE = []


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

        # the following are used for energy calculations 
        #at the half time step 

        v0_apoints[i].append(v0[i][0])  # x velocity component 

        v0_bpoints[i].append(v0[i][1])  # y velocity component

        v0_xpoints[i].append(v0[i][2])  # x position

        v0_ypoints[i].append(v0[i][3])  # y position



        
# below is the actual calculations from the verlet method as described in the lab handout
    r= r + h * v0 

    k= h * f(r, t + h)

    v= v0 + 0.5 * k

    v0= v0 + k



# Calculate Energies over time for each of the 16 particles

# Kinetic energy 

KE = [[] for i in range((N))]
 
#creating for loops to go through the particles and time steps 
for i in range(N):

    for time_step in range(len(v0_apoints[i])):

        KEtime = 1 / 2 * (v0_apoints[i][time_step]**2 + v0_bpoints[i][time_step]**2) #from the equation KE = 1/2mv**2

        KE[i].append(KEtime) #adding the new kinetic energy of iteration to the list

KE_total = [] #calculating the total KE 

for dt in range(len(tpoints)):

    Accumulator = 0

    for j in range(len(KE)):

        Accumulator += KE[j][dt]

    KE_total.append(Accumulator)

    

# Potential energy 

PE = [[] for i in range((N))]
#creating for loops to go through the particles and time steps 
for i in range(N):

    for time_step in range(len(v0_apoints[i])):

        Acc = 0

        for j in range(N):

            ri = (v0_xpoints[i][time_step]**2 + v0_ypoints[i][time_step]**2)**0.5

            rj = (v0_xpoints[j][time_step]**2 + v0_ypoints[j][time_step]**2)**0.5

            r = rj - ri #the difference between distance for a particle i and j 

            if r >= 0.0000001: #if this difference is greater than approximately 0 

                PEtime = (2) * ((1 / (r)**12) - (1 / (r))**6) # from equation for potential in the lab handout

                Acc += PEtime #accumulating the potential energy 

        PE[i].append(Acc)

#now get the total potential energy 
PE_total = []

for dt in range(len(tpoints)):

    Accumulator = 0

    for j in range(len(PE)):

        Accumulator += PE[j][dt]

    PE_total.append(Accumulator)


# below is plotting for question 3 a
# figure(1)

# for i in range(N):

#     plot(xpoints[i], ypoints[i], '.')



# xlabel("x")

# ylabel("y")

# xlim(-10, 20)

# ylim(-10, 20)

# title('The position of 16 particles')

# show()




#plotting for just the kinetic energy 
# figure(2)

# plt.plot(arange(T), KE_total, label='Kinetic')

# xlabel("Time (s)")

# ylabel("Energy")

# ylim(-100, 100)

# title('The energy of the N body system over time')

# show()


#plotting for just the potential energy 
# figure(3)

# plt.plot(arange(T), PE_total, label='Potential')

# xlabel("Time (s)")

# ylabel("Energy")

# ylim(-100, 100)

# title('The energy of the N body system over time')

# show()


#creating a for loop to sum the energies 
total = []

for i in range(len(PE_total)):

    total.append(PE_total[i] + KE_total[i])


#plotting the total energy 
figure(4)

plt.plot(arange(T), (total), label='Total')

xlabel("Timesteps (dt = 0.01 s)")

ylabel("Energy")

ylim(-100, 700)

title('The energy of the N body system over time')

show()

