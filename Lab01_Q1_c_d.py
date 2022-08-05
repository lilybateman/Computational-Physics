#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

PHY407: Assignment 1, Question 1. c) and d)

"""



import numpy as np
from numpy import zeros
from matplotlib import pyplot as plt #to make the plots

#1.c)

# defining constants
t0 = 0 #unit is Earth year, initial time
tf = 1 #unit is Earth year, final time
vx0 = 0 #unit is AU/year, intial velocity x component
x0 = 0.47 #unit is AU, initial position x component
vy0 = 8.17 #unit is AU/year, initial velocity y component
y0 = 0 #unit is AU, initial position y component
r0 = (0.47)**(1/2) #unit is AU, this is what it is because r0 = sqrt(x0^2+y0^2) and y0 is 0
Ms = 1 #unit of mass of the sun because it is mass of the sun
G = 39.5 #unit of AU^2Ms^-1Earthyear^-2, gravitional constant

TimeStep = 0.0001 #unit is Earth year, this is delta t, the amount the time increases for each iteration through the for loop

n =10001 #the amount of iterations, found by 1/0.0001+1 the plus one because the iterations start at index 1 

t = np.linspace(t0,tf,num=n) #making time array, starting at 0 going to 1 earth year, with n subsections
x = zeros ([n]) #this and the next 4 rows are creating arrays of 0 the length of the amount of iterations that will be filled up in the for loop
vx = zeros ([n]) 
y = zeros ([n]) 
vy = zeros ([n]) 
r = zeros ([n]) 

vx[0] = vx0 #this and the next 4 rows are to tell the code the initial value of the arrays that were just created above
x[0] = x0
vy[0] = vy0
y[0] = y0
r[0] = r0 


for i in range(1,n): #creating the for loop from 1 up to n because initial values are already stated
    
    r[i-1] = ((x[i-1]**2+y[i-1]**2)**(1/2)) #r is indexed at i-1 so the first r0 can be passed into the vx[i] and so on, it could theorectially be at the end of the loop too with just r[i]
    
    vx[i] = TimeStep*(-G*Ms*(x[i-1]))/((r[i-1])**3) + vx[i-1] # filling in the x component of velocity, the indexes after the = are all[i-1] so the inital values can be used to build upon

    vy[i] = TimeStep*(-G*Ms*(y[i-1]))/((r[i-1]**3)) + vy[i-1] #same as the above line but for the velocity y component

    x[i] = TimeStep*vx[i]+x[i-1] # this is to fill the x component of position however take note the index of vx after the = is [i] so the velocity that was just created in this same iteration is used to conserve approximate energy

    y[i] = TimeStep*vy[i]+y[i-1] # same as above line but for y component of position

    #now the for loop will use the index values it just added to the zero arrays to build the next set of values and in turn, integrate


#plotting the x component of velocity vs time, y component of velocity vs time and the x and y position (orbit), respectively
plt.plot(t, vx)
plt.xlabel('Time (Earth Year)')
plt.ylabel('Velocity (AU/Earth Year)')
plt.title('Newtonian Velocity of Mercury (X Component) vs. Time') 
plt.show()

plt.plot(t, vy)
plt.xlabel('Time (Earth Year)')
plt.ylabel('Velocity (AU/Earth Year)')
plt.title('Newtonian Velocity of Mercury (Y Component) vs. Time') 
plt.show()

plt.plot(x, y)
plt.xlabel('X Distance between Mercury and the sun (AU)')
plt.ylabel('Y Distance between Mercury and the sun (AU) ')
plt.title('Newtonian Orbit of Mercury') 
plt.show()


Px = [] #making emptylist for x component of momentum
Py = [] #making emptylist for y component of momentum
for i in range(0, len(x)): #creating a for loop to multiply each element corresponding to the same index 
    Px.append( x[i]*r[i] ) #this is x component of momentum because position*radius for every value
    Py.append( y[i]*r[i] ) #same as line above but for y 
    

L = [] #angular momentum list
for i in range(0, len(x)):
    L.append( x[i]*Py[i] - y[i]*Px[i] ) #filling in the list with cross product expression L = rXP rXP = xP-yP

 
plt.ylim(-1, 1) #fixing the axis so it only goes from -1 to 1 (arbitrary, could have been bigger)
plt.plot(t, L)
plt.xlabel('Time (Earth year)')
plt.ylabel('Angular Momentum (Divided by M) (AU^2/Earth Year) of Mercury')
plt.title('Angular Momentum of Mercury vs. Time ') 
plt.show()





#1.d)

#Now do the for loop with general relativity gravitational force, the only difference is the function for vx and vy in the for loop now contains the term (1+alpha/r^2)
# Using values that are already listed above such as the same n, same initial velocities and positions, or function constants but did not want to list again and add clutter

x2 = zeros ([n])  
vx2 = zeros ([n]) 
y2 = zeros ([n]) 
vy2 = zeros ([n]) 
r2 = zeros ([n]) 


vx2[0] = vx0
x2[0] = x0
vy2[0] = vy0
y2[0] = y0
r2[0] = r0


alpha = 0.01 #Units of AU^2



for i in range(1,n):
    
    r2[i-1] = ((x2[i-1]**2+y2[i-1]**2)**(1/2))
      
    vx2[i] = TimeStep*(((-G*Ms*(x2[i-1]))/((r2[i-1])**3))*(1+alpha/(r2[i-1])**2)) + vx2[i-1] #same as previous code except extra term as mentioned in the first comment before this section

    vy2[i] = TimeStep*(((-G*Ms*(y2[i-1]))/((r2[i-1]**3)))*(1+alpha/(r2[i-1])**2)) + vy2[i-1]

    x2[i] = TimeStep*vx2[i]+x2[i-1]

    y2[i] = TimeStep*vy2[i]+y2[i-1]





plt.plot(t, vx2)
plt.xlabel('Time (Earth Year)')
plt.ylabel('Velocity (AU/Earth Year)')
plt.title('General Relativity Velocity of Mercury (X Component) vs. Time') 
plt.show()

plt.plot(t, vy2)
plt.xlabel('Time (Earth Year)')
plt.ylabel('Velocity (AU/Earth Year)')
plt.title('General Relativity Velocity of Mercury (Y Component) vs. Time') 
plt.show()

plt.plot(x2, y2)
plt.plot(-0.428, 0.09,'ro', label = 'Perihelion orbit 1') #approximated perihelion point for first orbit
plt.plot(-0.29, 0.065,'go', label = 'Perihelion orbit 4') #approximated perihelion point for fourth orbit
plt.xlabel('X Distance between Mercury and the sun (AU)')
plt.ylabel('Y Distance between Mercury and the sun (AU) ')
plt.title('General Relativity Orbit of Mercury') 
plt.legend()
plt.show()
  
  
  
  
  
