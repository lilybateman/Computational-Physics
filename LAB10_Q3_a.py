import matplotlib.pyplot as plt
import numpy as np



#a)


#importance sampling method

N = 10000 #amount of values to sum through

I_imp = [] #empty list to load I values into 

def f(x): #define main function (under integral )
    return (x**(-1/2))*(1/(np.e**(x)+1))

def w(x): #define given importance function w(x) 
    return (x**(-1/2))

for i in range(100): # want to redo the calculation 100 times
    
    z1 = np.random.rand(N) #random numbers between 0 and 1, with a length of N 
    
    x1 = (z1)**2 #this was found by integrating the probability density function 

    Sum1 = np.sum(f(x1)/w(x1)) #for each x at length N, it will be the input of the function and all of those are summed for f(x) over w(x) by importance formula

    I1 = (2/N)*(Sum1) #scale the sum by 2 (because that is value of the integral of w(x) from 0 to 1) and then N because that is how much random numbers we summed over
    
    I_imp.append(I1) #append each caluclated I to a list 


plt.hist(I_imp, 10, range=[0.80, 0.88]) #plot the histogram with criteria given in lab handout
plt.title('Histogram for Importance Method (a)')
plt.xlabel('Integral Value (I)')
plt.ylabel('Count')
plt.show()



#regular mean value method 

a = 0 #bounds of integration 
b = 1

I_reg = [] #create an empty list to later append values to 
for i in range(100): #redo the calculation 100 times

       x2 = np.random.random_sample(N) #for this method it samples just values between a and b and there is not the probability density process 
       
       Sum2 = np.sum(f(x2)) #sum f of all x in the length N of random numbers 

       I2 = (((b-a)/N)*(Sum2)) #scale by the difference of bounds of integration and N 
    
       I_reg.append(I2) #append to list 
    
        
plt.hist(I_reg, 10, range=[0.80, 0.88])  #plot the histogram with the criteria given in lab handout
plt.title('Histogram for Regular Mean Value Method (a)')
plt.xlabel('Integral Value (I)')
plt.ylabel('Count')
plt.show()
    
    
    