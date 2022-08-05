from pylab import plot, show, legend, subplot, xlabel, ylabel, tight_layout, figure, xlim, ylim, title
from numpy import zeros, empty, linspace, exp, arange, minimum, pi, sin, cos, array
from dcst import dst, idst, dct, idct
from numpy.linalg import solve
import cmath
from numpy import exp, linspace, matmul, dot, arange, diag
from scipy import integrate
from intfuncs import trapz_int
import numpy as np

###############################################################################
# 1a/b
# the infinite square well case

P = 1024 # number of segments over the x interval
L = 10**-8 # the width of the potential
m = 9.109*10**-31 # mass of the electron
sigma = L / 25 # a constant to make sure the initial wave function is correct 
k = 500 / L # a constant that determines the periodicity
N = 3000 # the number of time steps
a = L / P # the setp size in the xdirection
tou = 10**-18 # the length of a time step
hbar = 1.0545718 * 10 ** - 34  # plank's constant in units m2 kg / s
x0 = L / 5 # the initial x value

A = - (hbar)**2 / (2 * m * (a**2)) # part of the H matrix
B = -2 * A # part of the H matrix

# the potnetial energy function of the infinite square well
def V(x): 
    return 0

# the normalization constant of the wave function
def f(x_value):
    return np.exp(-(x_value - x0)**2 / 2 * (sigma**2))

w_int = trapz_int(-L/2, L/2, 1024, f)
print('w int', w_int)
# w0 = 1/((w_int)**0.5) 
w0 = 1/((1e-08)**0.5) 

# the list of x value
x = linspace(-L / 2, L / 2, P) 
xpoints = linspace(-L / 2, L / 2, P) 

# the wave function
w = w0 * exp(-(x - x0)**2 / (4 * (sigma)**2) + 1j * k * x)  
# the complex conjugate of the wave function
wstar = w0 * exp(-(x - x0)**2 / (4 * (sigma)**2) - 1j * k * x)

# define the discretized Hamiltonian matrix
H = zeros((P, P), complex)
for l in range(len(H)):
    for k in range(len(H)):
        if l == k:
            H[l][k] = B # on the diagnal
            
        elif (l == k + 1) or (l == k - 1): 
            H[l][k] = A # just above and below the diagnal
            
# define the identity matrix
I = zeros((P, P), complex)
for l in range(len(I)):
    for k in range(len(I)):
        if l == k:
            I[l][k] = 1

# define the L matrix
L = I + (1j * H*(tou / (2 * hbar)))
# define the R matrix
R = I - (1j * H*(tou / (2 * hbar)))

# define the v vector
v = matmul(R, w)
# solve for x
x = solve(L, v)

# intialize lists of Energy and time averaged x position points
Epoints = []
Ppoints = []

# define starting time
t = 0
# initialize list of time points
tpoints = []
# define maximum time
T = N*tou

# initialize a list for keeping track of the normalization over time
w0_points = []
# initialize a list for keeping track of the probability density over time
Prob_points = []

# loop through time
while t <= (N * tou):
    
    # solve for the updated wave function
    # Lw = v where v = Rw
    v = matmul(R, w)
    x = solve(L, v)
    for i in range(len(w)):
        w[i] = x[i]

    tpoints.append(t)
    
    #energy 
    intermediate = matmul(H,w)
    energy = dot(wstar, intermediate)
    Epoints.append((energy))
    
    #expected position
    intermediate_p = dot(x, w)
    position = wstar * intermediate_p
    Ppoints.append(integrate.simps(position))
    
    #expected probability density
    probability = dot(wstar, w)
    Prob_points.append((probability))
    
    # for plotting the normalization of the wavefunction
    vals = f(x)
    w_int = integrate.simps(vals)
    w0_points.append(w_int)

    # plots for part b 
    
    if t == 0:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(1)
        title('The Real Part of the Wave Function at t=0 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
    
    elif t == 7.5000000000000775e-16:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(2)
        title('The Real Part of the Wave Function at t=T/4 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
        
    elif t == 1.5000000000000198e-15:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(3)
        title('The Real Part of the Wave Function at t = T/2 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
        
    elif t == 2.999999999999803e-15:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(4)
        title('The Real Part of the Wave Function at t = T')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
            
    t += tou

# plots for part a 

# plot of the normalized wave function with time
figure(5)
title('The Normalization of the Wave Function with Time')
xlabel('Time')
ylabel('The normalization')
plot(tpoints, w0_points)
show()

# plot of the energy with time
figure(6)
title('Energy as a function of Time')
xlabel('Time')
ylabel('Energy')
plot(tpoints, Epoints)
show()

# for part b 

# the trajectory
figure(7)
title('The Trajectory: <x>(t)')
xlabel('Time')
ylabel('Position')
plot(tpoints, Ppoints)
show()

# the probability density
figure(8)
title('The probability density')
xlabel('Time')
ylabel('The probabilty density')
plot(tpoints, Prob_points)
show()

###############################################################################
# part c
# Harmonic osscilator

P = 1024 # number of segments over the x interval
L = 10**-8 # the width of the potential
m = 9.109*10**-31 # mass of the electron
sigma = L / 25 # a constant to make sure the initial wave function is correct 
k = 500 / L # a constant that determines the periodicity
N = 4000 # the number of time steps
a = L / P # the setp size in the xdirection
tou = 10**-18 # the length of a time step
hbar = 1.0545718 * 10 ** - 34  # plank's constant in units m2 kg / s
x0 = L / 5 # the initial x value

p = linspace(0, P, P) # initialize the p values 1, ... P-1


omega = 3e15

# initialize the potential energy function for the Harmonic osscilator
def V(arg):
    return 0.5*m*(omega**2)*(arg**2) 


A = - (hbar)**2 / (2 * m * (a**2))  # part of the H matrix

# part of the H matrix
B = []
for i in range(len(p)):
    arg = p[i]*a -L/2 
    B.append(V(arg) -(2 * A))
    
vec_diag = B
B_diag = diag(vec_diag, k=0)

def f(x_value):
    return np.exp(-(x_value - x0)**2 / 2 * (sigma**2))

w0 = 1/((1e-08)**0.5) # normalization constant

# the list of x value
x = linspace(-L / 2, L / 2, P) 
xpoints = linspace(-L / 2, L / 2, P) 

# the wave function
w = w0 * exp(-(x - x0)**2 / (4 * (sigma)**2) + 1j * k * x)  
# the complex conjugate of the wave function
wstar = w0 * exp(-(x - x0)**2 / (4 * (sigma)**2) - 1j * k * x)

# define the discretized Hamiltonian matrix
H = zeros((P, P), complex)
for l in range(len(H)):
    for k in range(len(H)):
        if l == k:
            H[l][k] = B_diag[l][k]
            
        elif (l == k + 1) or (l == k - 1):
            H[l][k] = A
            
# define the identity matrix
I = zeros((P, P), complex)
for l in range(len(I)):
    for k in range(len(I)):
        if l == k:
            I[l][k] = 1

# define the L matrix
L = I + (1j * H*(tou / (2 * hbar)))
# define the R matrix
R = I - (1j * H*(tou / (2 * hbar)))

# define the v vector
v = matmul(R, w)
# solve for x
x = solve(L, v)

# intialize lists of Energy and time averaged x position points
Epoints = []
Ppoints = []

# define starting time
t = 0
# initialize list of time points
tpoints = []
# define maximum time
T = N*tou

# initialize a list for keeping track of the normalization over time
w0_points = []
# initialize a list for keeping track of the probability density over time
Prob_points = []

# loop through time
while t <= (N * tou):
    
    # solve for the updated wave function
    # Lw = v where v = Rw
    v = matmul(R, w)
    x = solve(L, v)
    for i in range(len(w)):
        w[i] = x[i]

    tpoints.append(t)
    
    #energy 
    intermediate = matmul(H,w)
    energy = dot(wstar, intermediate)
    Epoints.append((energy))
    
    #expected position
    intermediate_p = dot(x, w)
    position = wstar * intermediate_p
    Ppoints.append(integrate.simps(position))
    
    #expected probability density
    probability = dot(wstar, w)
    Prob_points.append((probability))
    
    # for plotting the normalization of the wavefunction
    vals = f(x)
    w_int = integrate.simps(vals)
    w0_points.append(w_int)

    # plots for part b 
    
    if t == 0:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(1)
        title('The Real Part of the Wave Function at t=0 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
    
    elif t == 7.5000000000000775e-16:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(2)
        title('The Real Part of the Wave Function at t=T/4 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
        
    elif t == 1.5000000000000198e-15:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(3)
        title('The Real Part of the Wave Function at t = T/2 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
        
    elif t == 2.999999999999803e-15:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(4)
        title('The Real Part of the Wave Function at t = T')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
            
    t += tou

# plots for part a 

# plot of the normalized wave function with time
figure(5)
title('The Normalization of the Wave Function with Time')
xlabel('Time')
ylabel('The normalization')
plot(tpoints, w0_points)
show()

# plot of the energy with time
figure(6)
title('Energy as a function of Time')
xlabel('Time')
ylabel('Energy')
plot(tpoints, Epoints)
show()

# for part b 

# the trajectory
figure(7)
title('The Trajectory: <x>(t)')
xlabel('Time')
ylabel('Position')
plot(tpoints, Ppoints)
show()

# the probability density
figure(8)
title('The probability density')
xlabel('Time')
ylabel('The probabilty density')
plot(tpoints, Prob_points)
show()

################################################################################
# part d

P = 1024 # number of segments over the x interval
L = 10**-8 # the width of the potential
m = 9.109*10**-31 # mass of the electron
sigma = L / 25 # a constant to make sure the initial wave function is correct 
k = 500 / L # a constant that determines the periodicity
N = 6000 # the number of time steps
a = L / P # the setp size in the xdirection
tou = 10**-18 # the length of a time step
hbar = 1.0545718 * 10 ** - 34  # plank's constant in units m2 kg / s
x0 = L / 3 # initial x
x1 = L/4 # constant used in defining the potential
V0 = 6e-17 # the initial potential energy

p = linspace(0, P, P) # initialize the p values 1, ... P-1


omega = 3e15
def V(arg):
    return V0*((((arg**2/x1**2) - 1))**2)


A = - (hbar)**2 / (2 * m * (a**2))
B = []
for i in range(len(p)):
    arg = p[i]*a -L/2 
    B.append(V(arg) -(2 * A))
    
vec_diag = B
B_diag = diag(vec_diag, k=0)

def f(x_value):
    return np.exp(-(x_value - x0)**2 / 2 * (sigma**2))

w0 = 1/((1e-08)**0.5) # normalization constant

# the list of x value
x = linspace(-L / 2, L / 2, P) 
xpoints = linspace(-L / 2, L / 2, P) 

# the wave function
w = w0 * exp(-(x - x0)**2 / (4 * (sigma)**2) + 1j * k * x)  
# the complex conjugate of the wave function
wstar = w0 * exp(-(x - x0)**2 / (4 * (sigma)**2) - 1j * k * x)

# define the discretized Hamiltonian matrix
H = zeros((P, P), complex)
for l in range(len(H)):
    for k in range(len(H)):
        if l == k:
            H[l][k] = B_diag[l][k]
            
        elif (l == k + 1) or (l == k - 1):
            H[l][k] = A
            
# define the identity matrix
I = zeros((P, P), complex)
for l in range(len(I)):
    for k in range(len(I)):
        if l == k:
            I[l][k] = 1

# define the L matrix
L = I + (1j * H*(tou / (2 * hbar)))
# define the R matrix
R = I - (1j * H*(tou / (2 * hbar)))

# define the v vector
v = matmul(R, w)
# solve for x
x = solve(L, v)

# intialize lists of Energy and time averaged x position points
Epoints = []
Ppoints = []

# define starting time
t = 0
# initialize list of time points
tpoints = []
# define maximum time
T = N*tou

# initialize a list for keeping track of the normalization over time
w0_points = []
# initialize a list for keeping track of the probability density over time
Prob_points = []

# loop through time
while t <= (N * tou):
    
    # solve for the updated wave function
    # Lw = v where v = Rw
    v = matmul(R, w)
    x = solve(L, v)
    for i in range(len(w)):
        w[i] = x[i]

    tpoints.append(t)
    
    #energy 
    intermediate = matmul(H,w)
    energy = dot(wstar, intermediate)
    Epoints.append((energy))
    
    #expected position
    intermediate_p = dot(x, w)
    position = wstar * intermediate_p
    Ppoints.append(integrate.simps(position))
    
    #expected probability density
    probability = dot(wstar, w)
    Prob_points.append((probability))
    
    # for plotting the normalization of the wavefunction
    vals = f(x)
    w_int = integrate.simps(vals)
    w0_points.append(w_int)

    # plots for part b 
    
    if t == 0:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(1)
        title('The Real Part of the Wave Function at t=0 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
    
    elif t == 7.5000000000000775e-16:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(2)
        title('The Real Part of the Wave Function at t=T/4 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
        
    elif t == 1.5000000000000198e-15:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(3)
        title('The Real Part of the Wave Function at t = T/2 s')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
        
    elif t == 2.999999999999803e-15:
        wplot = []
        for i in range(len(w)):
            wplot.append(w[i].real)
        figure(4)
        title('The Real Part of the Wave Function at t = T')
        xlabel('Position')
        ylabel('The Real Part of the Wave Function')
        plot(xpoints, wplot)
        show()
            
    t += tou

# plots for part a 

# plot of the normalized wave function with time
figure(5)
title('The Normalization of the Wave Function with Time')
xlabel('Time')
ylabel('The normalization')
plot(tpoints, w0_points)
show()

# plot of the energy with time
figure(6)
title('Energy as a function of Time')
xlabel('Time')
ylabel('Energy')
plot(tpoints, Epoints)
show()

# for part b 

# the trajectory
figure(7)
title('The Trajectory: <x>(t)')
xlabel('Time')
ylabel('Position')
plot(tpoints, Ppoints)
show()

# the probability density
figure(8)
title('The probability density')
xlabel('Time')
ylabel('The probabilty density')
plot(tpoints, Prob_points)
show()




