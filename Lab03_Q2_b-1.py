import numpy as np
import math as m
import matplotlib.pyplot as plt

def Hermite(n,x):
    
    def Hermite_LIP(X,Y,Z):
        if Z==0 :
            return Y
        else:
            return Hermite_LIP(2 * x * X - 2 * (Z - 1) * Y, X, Z - 1)
    
    return Hermite_LIP(x*2, 1, n)

def wavefunction(n, x):
    return 1 / m.sqrt(2 ** n * m.factorial(n) * m.sqrt(m.pi)) * m.exp(- x ** 2 / 2) * Hermite(n, x)

position = np.linspace(-10, 10, 1000)
Energy30 = []
for k in range(1000):
    Energy30.append(wavefunction(30, position[k]))

plt.plot(position, Energy30)
plt.xlabel("position (x)")
plt.ylabel("Energy Levels (H)")
plt.title("H30")


    