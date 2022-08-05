#Pseudocode:
#Define the Hermite Polynomial for given x and any n > 0 and n = 0 using Linear Iterative Process (LIP)
#Plot Hermite Polynomial for n = 0, 1, 2, 3, for x between -4 to 4

import numpy as np
import matplotlib.pyplot as plt

def Hermite(n,x):
    
    def Hermite_LIP(X,Y,Z):
        if Z==0 :
            return Y
        else:
            return Hermite_LIP(2 * x * X - 2 * (Z - 1) * Y, X, Z - 1)
    
    return Hermite_LIP(x*2, 1, n)


position = np.linspace(-4,4, 1000)
Energy0= []
Energy1= []
Energy2= [] 
Energy3= []
for k in range(1000):
    Energy0.append(Hermite(0,position[k]))
    Energy1.append(Hermite(1,position[k]))
    Energy2.append(Hermite(2,position[k]))
    Energy3.append(Hermite(3,position[k]))
    
plt.plot(position, Energy0)
plt.plot(position, Energy1)
plt.plot(position, Energy2)
plt.plot(position, Energy3)

plt.xlabel("Position (x)")
plt.ylabel("Energy Levels (H)")
plt.legend("0123")

