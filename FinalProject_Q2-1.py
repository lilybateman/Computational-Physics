

#Gauss-Seidel

#import module
import numpy as np


#define constants
N = 1000 #iteration value

dim = 3 #dimension of the matrix (3x3)


Arand = np.random.rand(dim,dim) #make a random matrix of desired size

A = np.copy(Arand) #copy it to get the same one


brand = np.random.rand(dim, 1) #make a random coloumn vector

b = np.copy(brand) #copy to get the same one 





x1 = np.zeros(((dim), (1))) #create a coloumn vector of zeros 
x = np.array(((x1),)*N) #turn this into N coloumn vectors 



#creating the diagonal matrix 
AD = np.zeros(((dim), (dim))) #create array of zeros the size of A
for i in range (dim): #coloumns 
    for j in range (dim): #rows
        if i == j: #if it is the diagonal element 
            AD[i][j] = np.diag(A)[i] #make it the diagonal of the matrix A



L = np.tril(A, k = 0) #set the top triangle to 0 
U = np.triu(A, k = 0) #set the bottom triangle to 0 


#set the diagonal elements to 0 for U
for i in range(dim): #coloumns
    for j in range(dim): #rows
        if i == j: #if it is a diagonal element
            U[i][i] = 0 #make it 0 
            



for i in range(N): #iteration loop 

       
        x[i] = np.dot((np.linalg.inv(L)), (b - np.dot(U, x[i-1]))) #calculate the new x based upon the previous one


        error = np.max(np.abs((x[i] - x[i-1])/(x[i])))
        
        if error < (0.001):

            print(x[i]) #see the chosen x 
        
            gsol = (x[i]) #keep it so do standard deviation

            break #stop the loop 
        else:
             continue #continue the loop to look for x 

    
nsol = np.linalg.solve(A, b) #use the numpy function to get a 'theorectial' value of x

print(nsol) #see this x 


diff = np.abs(gsol - nsol) #calculate the difference of the jacobian x and numpy x

print(diff) # see difference

print(np.std(diff)) #compute the standard deviation and see it 


    
    
    
    
    
    
    



    

