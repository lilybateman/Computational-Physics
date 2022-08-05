# import modules
import numpy as np
from random import random, randrange
from matplotlib import pyplot as plt
import matplotlib
from matplotlib import animation


# definng constants
kB = 1.0
T = 1.0
J = 1.0
num_dipoles = 20
N = 100
time = 1000


def energyfunction(J_, dipoles): #define the function to calculate the energy of the system for given dipoles
    energy = -J_*np.sum(dipoles[0:-1]*dipoles[1:])
    return energy


def acceptance(E, Enew): #define the function that will say whether spin states should be flipped or not
    
    p = np.exp(-(Enew - E)/(kB*T)) #boltzmann probability
    
    if Enew - E <= 0 or Enew - E > 0 and p > random(): #defining the conditions for the probability to either accept or reject a flip
        
        return 'accepted' #return this string, this will be used later in the nested for loop
    
    else:
        return 'rejected' 
    


dipoles = np.zeros(((20),(20))) #defining an array of dimensions 20 by 20 as asked, filled with zeros so it can be added to later

for i in range(20): #creating a nested for loop because two dimensions 
    for j in range(20):
        dipoles[i][j] = (2*np.random.randint(0, 2)-1) #producing a random value that is confined to be either -1 or 1





energy = []  # empty list to fill with energies
E_total = [] # empty list to fill with energy total 
magnet = []  # empty list to add magnetization to for each step
M_total = [] #the total magnetization which is the sum of the magnetization at each i in N 
Dipoles = np.zeros((((time)), (20),(20))) #create a three dimensional array so that the spin states at each respective time can be animated


for t in range(time):
    for i in range(N):
        picked_1st = randrange(num_dipoles)  #pick a random value between 1 and 20 which is one dimension of the dipoles array
        picked_2nd = randrange(num_dipoles)  #pick a second random value between 1 and 20 which is the other dimension of the dipoles array
        dipoles[picked_1st, picked_2nd] *= 1 #have the original dipoles defined to store the original energy
        E = energyfunction(J, dipoles) #calculate original energy
        dipoles[picked_1st, picked_2nd] *= -1  # propose to flip 
        Enew1 = energyfunction(J, dipoles)      # calculate the energy of this new state post flip


        if acceptance(E, Enew1) == 'rejected': #if new flip is rejected
            dipoles[picked_1st, picked_2nd] *= -1 #then flip that state back 
            Enew2 = energyfunction(J, dipoles)   #calculate this new energy
            energy.append(Enew2) 
            magnet.append(np.sum(dipoles)) #add to magnet list 


        else:
            energy.append(Enew1) 
            magnet.append(np.sum(dipoles)) #if the proposed flip is not rejected, (it is accepted) then add the original magnet 

            
    M_total.append(np.sum(magnet)) #add values to total magnetization of each N iteration
    Dipoles[t] = dipoles #add the 20 by 20 array at the correct index of time for the animation
    E_total.append(np.average(energy)) #append energy of system averages for each N
    
        

        
plt.plot(np.arange(time), M_total) #plot total magnetization over time
plt.title('The Total Magnetization Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Magnetization')
plt.show()

plt.plot(np.arange(time), E_total) #plot total magnetization over time
plt.title('The Total Energy Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.show()





fig = plt.figure()
im = plt.imshow(Dipoles[0], interpolation="none", vmin = -1, vmax = 1) #setting inital state of animation


def update(t):
    im.set_array(Dipoles[t]) #now for the time set the array of the dipoles at that given time (which is a 20 by 20 array)


ani = matplotlib.animation.FuncAnimation(fig, func=update, frames=time, repeat = False, interval=time) #create animation
plt.title("Spin States")
plt.show()


    



