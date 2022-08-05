import numpy as np



N = 1000000 #number to sum through 


def f(r): #define the function of the radius 
    
    rsquare = 0 #start with 0 to accumulate data
    for i in range(10): #range of 10 because 10 dimensions 
        rsquare += r[i]**2 #doing this so r squared is each of the values in r squared not the sum of all the values then squared
    
    if rsquare <= 1: #if the radius squared is less than or equal to 1 
        return 1 #return 1 
    
    elif rsquare > 1: #if the radius squared is greater than 1 it is out of bounds so return 0
        return 0  

    
Sum = 0 #set to 0 to then add values 
Sumsq = 0 #set to 0 to then add values 
for i in range(N): #create a for loop to go through N times 
    r = np.random.rand(10) #get 10 random values as r
    Sum += (f(r)) # accumulating f(r) values for the sum, for I and for error, either gives 0 or 1
    Sumsq += (f(r)**2) #this is for the error


I = (2**10)/N*Sum #scale the sum by 2^10 because 10 dimensions (volume of the unit square) and then over N the amount of values summed over


print(I)

#Error

f1 = (1/N)*(Sum) #normal sum scaled by 1/N 

f2 = (1/N)*(Sumsq) #sum of the f(r) values squared, scaled by 1/N

varf = f2 - f1**2 #the varience, which is the difference of the above f where the normal sum f1 is squared

sigma = ((2**10)/N)*((N*varf)**(1/2)) #full error sigma which is the varience multiplied by the amount of values summed through to the square root then scaled by again volume of unit square over N


print(sigma)




    




