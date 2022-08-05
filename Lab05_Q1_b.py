#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Lab 5, Question 1, b

"""


import numpy as np
from matplotlib import pyplot as plt #for graphing



#a)


dow = np.loadtxt("dow.txt", float) #loading the data into array


#plotting data
plt.plot(dow, )
plt.xlabel('Time (Days)')
plt.ylabel('Closing Values ')
plt.title('Store closing values over time') 
plt.show()



#b)


dow_fourier = np.fft.rfft(dow) #do fourier on the dow data



#c)


N = len(dow_fourier) #this is for creating the 10 percent data array


dow_fourier10 = np.zeros(N, complex) #creating empty 0 array to load the 10 percent of fourier coefficiants into

for i in range(len(dow_fourier)): 
    if i < N*0.1: #this is ten percent of the length of the fourier of the dow so if i is less than this number, it will be for the first 10 percent of the data
        dow_fourier10[i] = dow_fourier[i] #keep the values
    else: dow_fourier10[i] = 0 #make them 0 if they are beyond the 10 percent

    
    
#d)

        
dow_fourier10i = np.fft.irfft(dow_fourier10) #now do the inverse of the fourier


#plotting data
plt.plot(dow,  label = 'Closing Values')
plt.plot(dow_fourier10i,  label ='Closing Values Fourier 10%')
plt.xlabel('Time (Days)')
plt.ylabel('Closing Values ')
plt.title('Store closing values over time') 
plt.legend()
plt.show()



#e)


#below is the same process as in d) above except it is 2 percent instead of 10 percent
dow_fourier2 = np.zeros(N, complex)

for i in range(len(dow_fourier)):
    if i < N*0.02:
        dow_fourier2[i] = dow_fourier[i]
    else: dow_fourier2[i] = 0
        

dow_fourier2i = np.fft.irfft(dow_fourier2)


#plotting data
plt.plot(dow,  label = 'Closing Values')
plt.plot(dow_fourier10i,  label ='Closing Values Fourier 10%')
plt.plot(dow_fourier2i,  label ='Closing Values Fourier 2%')
plt.xlabel('Time (Days)')
plt.ylabel('Closing Values ')
plt.title('Store closing values over time') 
plt.legend()
plt.show()