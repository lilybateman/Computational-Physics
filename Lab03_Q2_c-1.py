import math as m
from gaussxw import gaussxwab

def Hermite(n,x):
    
    def Hermite_LIP(X,Y,Z):
        if Z==0 :
            return Y
        else:
            return Hermite_LIP(2 * x * X - 2 * (Z - 1) * Y, X, Z - 1)
    
    return Hermite_LIP(x*2, 1, n)


def wavefunction(n, x):
    return 1 / m.sqrt(2 ** n * m.factorial(n) * m.sqrt(m.pi)) * m.exp(- x ** 2 / 2) * Hermite(n, x)

def position(x):
    def X(x):
        return x / (1 - x)

    return X(x) ** 2 * abs(wavefunction(5, X(x))) ** 2 * (1 / (1 - x) ** 2)


s = 0
N = 100
x, w = gaussxwab(N, 0, 1)
for k in range(N):
    s += w[k] * position(x[k])
print("The Position x is: ", m.sqrt(2 * s))

def momentum(p):
    def P(p):
        return p / (1 - p)
    
    return P(p) ** 2 * abs(wavefunction(5, P(p))) ** 2 * (1 / (1 - p) ** 2)

s = 0
N = 100
p, w = gaussxwab(N, 0, 1)
for k in range(N):
    s += w[k] * momentum(p[k])
print("The Momentum p is: ", m.sqrt(2 * s))

print("The Energy is: ", 0.5*(m.sqrt(2 * s)**2 + m.sqrt(2 * s)**2))


