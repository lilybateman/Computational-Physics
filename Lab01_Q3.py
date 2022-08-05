#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

PHY407: Assignment 1, Question 3

"""

import numpy as np 
from numpy import ones
from numpy import zeros
from matplotlib import pyplot as plt #to make the plots
import time # to time how long the code takes 




N = 2 #Setting the initial size of the matrix, it will be a 2x2, as desired in the question

times = [] # setting an empty list to load the times from each iteration of the for loop into
while N < 100: #this means the dimensions of the matrices being multplied together will be less than 100
    TimeStart = time.time() #starting the timer

    A = ones([N, N], float)*4 #making the first matrix to be multplied, it is comprised of only 4s, arbitarily
    B = ones([N, N], float)*5 #making the first matrix to be multplied, it is comprised of only 5s, arbitarily
    C = zeros([N, N], float) #making C the resultant matrix out of 0s so that content can be put into it via the for loop

    for i in range(N): 
        for j in range(N):
            for k in range(N): #there are this many for loops so all indexes are covered and specific elements corresponding to specific rows or coloumns can be multiplied
                C[i,j] += A[i,k]*B[k,j] #the actual multiplication of elements to put into C the zero array
                
    N = N + 1 #to make sure N increases with each iteration through the for loop
    TimeEnd = time.time() #ending the timer
    times.append(TimeEnd - TimeStart) #appending the difference of time, so the total time for each iteration into the time list 
 
    


N = 2

times2 = [] 

while N < 100:
    TimeStart = time.time()

    A = ones([N, N], float)*4
    B = ones([N, N], float)*5
    C = np.dot(A, B) #for this portion, a lot is similar to above except there is no for loop needed to use the numpy dot operation, it already knows which indexes to multiply
        
    N = N + 1
    TimeEnd = time.time()
    times2.append(TimeEnd - TimeStart)
    
    


N_arr = np.arange(2, 100) #creating the array of N values inorder to plot
N3_arr = (np.arange(2, 100))**3 #creating array of N^3 values inorder to plot



#will now be plotting all the graphs
plt.plot(N_arr, times)
plt.xlabel('Matrix Size (N)')
plt.ylabel('Time (s)')
plt.title('Manual Code: N vs. Time') 
plt.show()

plt.plot(N3_arr, times)
plt.xlabel('Matrix Size (N^3)')
plt.ylabel('Time (s)')
plt.title('Manual Code: N^3 vs. Time') 
plt.show()

plt.plot(N_arr, times2)
plt.xlabel('Matrix Size (N)')
plt.ylabel('Time (s)')
plt.title('Numpy.dot Operation: N vs. Time') 
plt.show()

plt.plot(N3_arr, times2)
plt.xlabel('Matrix Size (N^3)')
plt.ylabel('Time (s)')
plt.title('Numpy.dot Operation: N^3 vs. Time') 
plt.show()



