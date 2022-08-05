

#Jacobian

#import module
import numpy as np

#define constants
N = 1000 #iteration value

dim = 3 #dimension of the matrix (3x3)

Arand = np.random.rand(dim,dim) #make a random matrix of desired size

A = np.copy(Arand) #copy it to get the same one


brand = np.random.rand(dim, 1) #make a random coloumn vector

b = np.copy(brand) #copy to get the same one 



#Jacobian Method

#make a diagonally dominant matrix


mat = False
while mat is False:
    list1 = []
    for i in range(dim): #for the range of the coloumns
        Sum = 0 #make a sum to compare it for each row
        for y in range(dim): #row the range of the rows
            if y != i: #if it is not a diagonal element
                Sum += (np.abs(A[i][y])) #sum the values of that row

        if np.abs(A[i][i]) < Sum: #is the sum of the diagonal is less than the sum of the entries in the same row
            list1.append(1) #note this 
            mat = False #condition is not met 
            Arand = np.random.rand(dim,dim) #create a new random matrix and try it again
            A = np.copy(Arand)

        else: #if the condition is met for that row
            continue #go to the next one
    if len(list1) < 1: #if values were in the list it means one of the rows had a larger sum 
        mat = True #then the condition is true and we can use this matrix




x1 = np.zeros(((dim), (1))) #create a coloumn vector of zeros 
x = np.array(((x1),)*N) #turn this into N coloumn vectors 
deltax1 = np.zeros(((dim), (1))) #same process for delta x 
deltax = np.array(((x1),)*N)



#creating the diagonal matrix 
AD = np.zeros(((dim), (dim))) #create array of zeros the size of A
for i in range (dim): #coloumns 
    for j in range (dim): #rows
        if i == j: #if it is the diagonal element 
            AD[i][j] = np.diag(A)[i] #make it the diagonal of the matrix A




for i in range(N): #iteration loop 

        deltax[i-1] = np.dot((np.linalg.inv(AD)), (b - np.dot(A, x[i-1]))) #the change in x component

        x[i] = x[i-1] + deltax[i-1] #the new computed x with the addition of the change in x 

        error = np.max(np.abs(((deltax[i-1])/(x[i]+deltax[i-1])))) #computing the error value
    
        if error < (0.001): #imposing a threshold 
            
            print(x[i]) #see what the chosen x is 
            
            jsol = (x[i]) #keep it so do standard deviation
        
            break #stop the loop 
        else:
            continue #continue the loop to look for x 
    
    


nsol = np.linalg.solve(A, b) #use the numpy function to get a 'theorectial' value of x

print(nsol) #see this x 


diff = np.abs(jsol - nsol) #calculate the difference of the jacobian x and numpy x

print(diff) # see difference

print(np.std(diff)) #compute the standard deviation and see it 





    
    
    
    

    

