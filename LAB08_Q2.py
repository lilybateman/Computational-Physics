# QUESTION 2

# code adapted form Newman's heat.py from the textbook

# import statements
from numpy import empty, exp, arange, copy, zeros, mean, linspace
from pylab import plot,xlabel,ylabel,show, figure, title

# constants
L = 1.0      # maximum distance in meters
J = 50       # number of grid spaces
a = L / J    # width of each grid space
h = 0.01     # time increment
epsilon = h/1000
g = 9.81     # gravity
H = 0.01     # in meters

# define different times we want to plot the graphs
t1 = 0.0  # in seconds
t2 = 1.0
t3 = 4.0

# define final time iteration
tfinal = t3 + epsilon

# Boundary and initial conditions

# initial u array
u = zeros(J + 1, float)
# array to hold new u values
up = zeros(J + 1, float)

# initial eta array
eta = zeros(J + 1, float)
# array to hold new eta values
etap = zeros(J+1,float)

# constants used to define initial conditions of eta
A = 0.002
mew = 0.5
sigma = 0.05
xvals = linspace(0, L, J+1)

constant = mean(A * exp(-(xvals - mew)**2 / sigma**2))                    
for i in range(len(xvals)):
    eta[i] = H + A* exp(-(xvals[i] - mew)**2 / sigma**2) - (constant) 

# define initial time
t = 0.0
c = -h/(2*a)

# iterate through time
while t < tfinal:
    
    # plot at t=0
    if t == 0.0:
        figure(1)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of the surface")
        title('FTCS Scheme: Shallow Water System at t = 0 s')
        show()

    # iterate through space
    for j in range(0, J + 1):  # 0, ..., 50


        # forward difference scheme for eta at the beginning
        if j == 0:
            etap[0] = eta[0] + c * (u[1] * eta[1] - u[0] * eta[0])
            up[0] = 0
        
        # backwards difference scheme for eta at the end
        elif j == (J):
            etap[j] = eta[j] + c * (u[j] * eta[j] - u[j - 1] * eta[j - 1])
            up[j] = 0
        
        # for interior points
        else:
            up[j] = u[j] + c*(g*eta[j+1] + u[j+1]**2 /2 - g*eta[j-1] - u[j-1]**2 /2)
            etap[j] = eta[j] + c*(u[j+1]*eta[j+1]-u[j-1]*eta[j-1])

    # use new values of u and eta for the next iteration
    eta = copy(etap)
    u = copy(up)
    
    # update time
    t += h

    # Make plots at the given times
    # plot at t=1
    if abs(t-t2)<epsilon:
        figure(2)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of the surface")
        title('FTCS Scheme: Shallow Water System at t = 1 s')
        show()

    # plot at t=4
    if abs(t - t3) < epsilon:
        figure(3)
        plot(xvals, eta)
        xlabel("x (m)")
        ylabel("Altitude of the surface")
        title('FTCS Scheme: Shallow Water System at t = 4 s')
        show()




