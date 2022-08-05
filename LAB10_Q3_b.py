import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
from scipy import special



#b)


#importance sampling method

N = 10000 #amount of values to sum through

I_imp = [] #empty list to load I values into 

def f(x): #define main function under integral 
    return (np.exp(-2*(np.abs(x-5))))

def w(x): #define importance function given 
    return (1/(2*np.pi)**(1/2))*np.exp((-(x-5)**2)/2)

wconst = sc.special.erf(5/(2**(1/2))) #this is the integral of w(x) from 0 to 10

x1 = np.zeros((N)) #empty array to load values into 
for i in range(100):
    for j in range(N): #create a for loop to have a unique value around 5 for each index of the N length array
        x1[j] = np.array(np.random.normal(5)) #because gaussian distribution want values centered around 5 

    Sum1 = np.sum(f(x1)/w(x1)) #take the sum of f(x) over w(x) inputting each of the x values generated above

    I1 = (wconst/N)*(Sum1) #scale by the integral of w(x) from 0 to 10 over N 

    I_imp.append(I1)


plt.hist(I_imp, 10, range=[0.95, 1.05])  #take a reasonable range 
plt.title('Histogram for Importance Method (b)')
plt.xlabel('Integral Value (I)')
plt.ylabel('Count')
plt.show()



#regular mean value method 

a = 0 #integration bounds
b = 10

I_reg = [] #empty list to put values into 

for i in range(100):
    
    x2 = (np.random.random_sample(N))*10 #multplying by 10 to have the values between a and b 

    Sum2 = np.sum(f(x2)) #summing over the values of N for f(x)

    I2 = ((b-a)/N)*(Sum2) #scale sum by the difference of integration bounds over N 

    I_reg.append(I2) #append values into list 
    
    
    
plt.hist(I_reg, 10, range=[0.95, 1.05])  #take the same range as importance method to better compare them 
plt.title('Histogram for Regular Mean Value Method (b)')
plt.xlabel('Integral Value (I)')
plt.ylabel('Count')
plt.show()
    
    