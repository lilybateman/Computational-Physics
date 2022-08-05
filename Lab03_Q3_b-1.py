#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Question 3. b)

"""

import struct # for reading binary files
import numpy as np
from matplotlib import pyplot as plt # for plotting



f = open("N46E006.hgt", 'rb') #opening the file


W = np.zeros([1201],[1201]) #creating w array of zeros to fill, 1201 by 1201 shape

for i in range(1201): #creating a nested for loop to go through vertical and horizontal 
    for j in range(1201):
        buf = f.read(2) #read 2 bytes of data each iteration for each index above
        value = struct.unpack('>h', buf)[0] #gives value of height when unpacked
        W[i, j] = value #fill up the W array with these heights corresponding to veritcal and horizontal coordinates


#forward for the right side and backward for the left side but I am not quite sure how to do this. 
#so I am doing a central derivative estimation for every point
 
#for x 

#Forward differnce scheme

h = np.array([0, 1201*1201, 83]) #step size array


dwdx_forward = np.zeros([len(h)]) #creating empty array to fill up in for loop 
for i in range(len(h)):
    dwdx_forward[i] = ((W(0.5+h[i], 0) - W(0.5, 0))/h[i])
    

#Central difference scheme

dwdx = np.zeros([len(h)])

for i in range(len(h)):
    dwdx[i] = (W(0.5+(h[i]/2), 0) - W(0.5 -(h[i]/2)), 0)/(h[i])
    
    
#Now for y 

#Forward differnce scheme
    
dwdy_forward = np.zeros([len(h)])
for i in range(len(h)):
    dwdy_forward[i] = ((W(0, 0.5+h[i]) - W(0.5, 0))/h[i])
    

    

#Central difference scheme

dwdy = np.zeros([len(h)])

for i in range(len(h)):
    dwdy[i] = (W(0, 0.5+(h[i]/2)) - W(0, 0.5 -(h[i]/2)))/(h[i])




I = np.zeros([1201], [1201]) #creating 2d intensity array

for i in range(len(dwdx)): #nested for loop to fill both dimensions
    for j in range(len(dwdy)):
        phi = np.pi/6
        I[i, j] = (np.cos(phi)*(dwdx[i])+np.sin(phi)*(dwdy[j]))/(((dwdx[i])**2+(dwdy[j])**2+1)**(1/2))





plt.imshow(W, extent = [46, 47, 6, 7]) #plotting the height density data
    
plt.imshow(I, extent = [0, 1201, 1201, 0]) #plotting the light intensity density data