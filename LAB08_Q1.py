# QUESTION 1
# code adapted from Newman's laplace.py

# import statements
from numpy import empty, zeros, max, meshgrid, arange, stack, linspace, gradient
from pylab import imshow, gray, show, title, figure, ylabel, xlabel
from matplotlib.pyplot import streamplot, figure, colorbar, tight_layout, axis
from time import time

# constants
M = 100         # number of grid squares per side
target = 1e-6   # desired accuracy
Vp = 1.0        # voltage at first plate
Vn = -1.0       # voltage at second plate

# Q1 PART a ###################################################################
# Gauss-Seidel method without relaxation

#start_firstmethod = time()

phi = zeros([M + 1, M + 1], float)  # initialize matrix of phi values to store potentials
phi[20:80, 20] = Vp  # define voltage of first plate
phi[20:80, 80] = Vn  # define voltage of second plate

delta = 1.0
while delta > target:  # calculate each grid point again untill accurate enough
    error_list = []  # initialize list of errors

    # calculate new potential at each point in the grid
    for i in range(M + 1):
        for j in range(M + 1):
            if i == 0 or i == M or j == 0 or j == M:  # if on boundary, remain the same
                phiprime = phi[i, j]
            elif 20 <= i <= 80 and (j == 20 or j == 80):  # if on one of the plates remain the same
                phiprime = phi[i, j]
            else:
                phiprime = (phi[i + 1, j] + phi[i - 1, j]
                                 + phi[i, j + 1] + phi[i, j - 1]) / 4
            each_delta = abs(phi[i, j] - phiprime)  # calculate error at grid point
            error_list.append(each_delta)  # add this error to list of errors
            phi[i, j] = phiprime  # replace with new value of potential

    # find maximum error from old values of potential
    delta = max(error_list)

# plot
figure(1)
imshow(phi)
gray()
xlabel('x')
ylabel('y')
title('Electrostatic Potential of an Electronic Capacitor')
show()

#end_firstmethod = time()
#difference_firstmethod = end_firstmethod - start_firstmethod
#print('difference_firstmethod=', difference_firstmethod)

###########
# plot the electric feild lines, each line a different colour representing a voltage
# E = - gradient(V)

x = linspace(0, 1, 101)
y = linspace(0, 1, 101)
X, Y = meshgrid(x, y)

V = phi
Ey, Ex = gradient(-V, y, x)

figure(2)
strm = streamplot(X, Y, Ex, Ey, color=V, linewidth=2, cmap='autumn')
cbar = colorbar(strm.lines)
cbar.set_label('Potential $V$')
title('Electric field lines')
xlabel('$x$')
ylabel('$y$')
axis('equal')
tight_layout()
show()

###############################################################################

#part b #######################################################################
# gauss seidel with relaxation # w = 0.1
#start_smethod = time()

phi = zeros([M + 1, M + 1], float)  # initialize matrix of phi values to store potentials
phi[20:80, 20] = Vp  # define voltage of first plate
phi[20:80, 80] = Vn  # define voltage of second plate
w = 0.1

delta = 1.0
while delta > target:  # calculate each grid point again untill accurate enough
    error_list = []  # initialize list of errors

    # calculate new potential at each point in the grid
    for i in range(M + 1):
        for j in range(M + 1):
            if i == 0 or i == M or j == 0 or j == M:  # if on boundary, remain the same
                phiprime = phi[i, j]
            elif 20 <= i <= 80 and (j == 20 or j == 80):  # if on one of the plates remain the same
                phiprime = phi[i, j]
            else:
                phiprime = ((phi[i + 1, j] + phi[i - 1, j]
                                 + phi[i, j + 1] + phi[i, j - 1]) * (1 + w) / 4) - w * phi[i, j]
            each_delta = abs(phi[i, j] - phiprime)  # calculate error at grid point
            error_list.append(each_delta)  # add this error to list of errors
            phi[i, j] = phiprime  # replace with new value of potential

    # find maximum error from old values of potential
    delta = max(error_list)

# plot
figure(3)
xlabel('x')
ylabel('y')
imshow(phi)
gray()
title('Electrostatic Potential of an Electronic Capacitor')
show()

#end_smethod = time()
#difference_smethod = end_smethod - start_smethod
#print('difference_smethod=', difference_smethod)

# gauss seidel with relaxation # w = 0.5 ######################################
#start_tmethod = time()

phi = zeros([M + 1, M + 1], float)  # initialize matrix of phi values to store potentials
phi[20:80, 20] = Vp  # define voltage of first plate
phi[20:80, 80] = Vn  # define voltage of second plate
w = 0.5

delta = 1.0
while delta > target:  # calculate each grid point again untill accurate enough
    error_list = []  # initialize list of errors

    # calculate new potential at each point in the grid
    for i in range(M + 1):
        for j in range(M + 1):
            if i == 0 or i == M or j == 0 or j == M:  # if on boundary, remain the same
                phiprime = phi[i, j]
            elif 20 <= i <= 80 and (j == 20 or j == 80):  # if on one of the plates remain the same
                phiprime = phi[i, j]
            else:
                phiprime = ((phi[i + 1, j] + phi[i - 1, j]
                                 + phi[i, j + 1] + phi[i, j - 1]) * (1 + w) / 4) - w * phi[i, j]
            each_delta = abs(phi[i, j] - phiprime)  # calculate error at grid point
            error_list.append(each_delta)  # add this error to list of errors
            phi[i, j] = phiprime  # replace with new value of potential

    # find maximum error from old values of potential
    delta = max(error_list)

# plot
figure(4)
imshow(phi)
xlabel('x')
ylabel('y')
gray()
title('Electrostatic Potential of an Electronic Capacitor')
show()

#end_tmethod = time()
#difference_tmethod = end_tmethod - start_tmethod
#print('difference_tmethod=', difference_tmethod)
###############################################################################



