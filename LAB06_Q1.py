
# import statements
from math import sin
from numpy import arange, array
from pylab import plot, xlabel, ylabel, show, xlim, figure, title

# 8.8
# part b : use RK4 to solve the equations

# want to solve dx/dt = a
# da/dt = -G * M * x / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))
# dx/dt = a
# dy/dt = b
# db/dt = -G * M * y / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))
# from t = 0 to t = 10, M = 10, G = 1, L = 2, IC : (x,y) = (1,0) @ t=0

# define constants
G = 1
L = 2
M = 10

# define the four fisrt order ODEs
def f(r, t):
    a = r[0]  # update x component velocity
    b = r[1]  # update y component velocity
    x = r[2]  # update x position
    y = r[3]  # update y position
    fx = a  # dx/dt = a
    fy = b  # dy/dt = b
    fb = -G * M * y / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))  # db/dt
    fa = -G * M * x / ((x**2 + y**2) * (((x**2 + y**2) + L**2 / 4)**0.5))  # da/dt
    return array([fa, fb, fx, fy], float)

# define the interval
start = 0.0
stop = 10.0
N = 1000
h = (stop - start) / N

# initialize lists
tpoints = arange(start, stop, h)
apoints = []  # this is a list of x component velcoities at each time step
bpoints = []  # this is a list of y component velcoities at each time step
xpoints = []  # this is a list of x positions at each time step
ypoints = []  # this is a list of y positions at each time step

# this array starts with the initial x component velocity, y component velocity, x position, yposition, in that order and updates it at each time iteration
r = array([0.0, 1.0, 1.0, 0.0], float)

# use RK4 to solve the ODEs
for t in tpoints:
    apoints.append(r[0])  # update list of x component velcoities at each time step
    bpoints.append(r[1])  # update list of y component velcoities at each time step
    xpoints.append(r[2])  # update list of x positions at each time step
    ypoints.append(r[3])  # update list of y positions at each time step

    # RK4 equations
    k1 = h * f(r, t)
    k2 = h * f(r + (0.5 * k1), t + (0.5 * h))
    k3 = h * f(r + (0.5 * k2), t + (0.5 * h))
    k4 = h * f(r + k3, t + h)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

# plot the orbit of the ball bearing
plot(xpoints, ypoints)
xlabel("x")
ylabel("y")
title('The orbit of the ball bearing')
show()

