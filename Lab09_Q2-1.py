from dcst import dct, idct
import numpy as np
from numpy import zeros
from numpy import empty,arange,exp,real,imag,pi
from numpy.fft import rfft, irfft
from matplotlib import pyplot as plt
from pylab import imshow, gray, show, title, figure, ylabel, xlabel

def dXXt2(f):
    """ Takes DXT along x, then DXT along y (X = C/S) IN: f, the input 2D numpy array
    OUT: b, the 2D transformed array """
    M = f.shape[0]  # Number of rows
    N = f.shape[1]  # Number of columns
    a = np.zeros((M, N))  # Intermediate array
    b = np.zeros((M, N))  # Final array

    # Take transform along x
    for j in range(N):
        # DXT f[:, j] and set as a[:, j]
        a[:, j] = dct(f[:, j])

    # Take transform along y
    for i in range(M):
        # DXT a[i, :] and set as b[i, :]
        b[i, :] = dct(a[i, :])

    return b


def idXXt2(b, f):
    """ Takes iDXT along y, then iDXT along x (X = C/S) IN: b, the input 2D numpy array
    OUT: f, the 2D inverse-transformed array """
    M = f.shape[0]  # Number of rows
    N = f.shape[1]  # Number of columns
    a = np.zeros((M, N), dtype=complex)  # Intermediate array
    f = np.zeros((M, N), dtype=complex)  # Final array
    
    # Take inverse transform along y
    for i in range(M):
        # iDXT b[i,:] and set as a[i,:]
        a[i, :] = idct(b[i, :])
    # Take inverse transform along x
    for j in range(N):
        # iDXT a[:,j] and set as f[:,j]
        f[:, j] = idct(a[:, j])
        
    return f

# the above code was adapted from the lab manual and the textbook's dcst.py

################################################################################
# part a 

N = 200#15 # should be 2000 # number of time steps
P = 32 # should be 32 # number of y and x steps
tou = 0.01 # duration of each time step

Lx = 1#P * ax the wall in the x direction
Ly = 1#P * ay the wall in the y direction
omega = 3.75 # driving frequency
J0 = 1 # constant
ax = Lx/P  # size of the cells in the x direction
ay = Ly/P# size of the cells in the y direction


m = 1
n = 1
c = 1

# initialize lists of points for space and time
nlist = (np.arange(0, N))
p = (arange(0, P))
q = (np.arange(0, P))

# disrcitize t, x, and y
t = tou * nlist
x = (ax * p) 
y = (ay * q)

# initialize matrices
Hx = zeros((P, P), complex)
Hy = zeros((P, P), complex)
Ez = zeros((P, P), complex)
Jz = zeros((P, P), complex)

# initialize matrices to hold fourier coefficients
E = zeros((P, P), complex)
X = zeros((P, P), complex)
Y = zeros((P, P), complex)
J = zeros((P, P), complex)

Dx = np.pi * c * tou / (2 * Lx)
Dy = np.pi * c * tou / (2 * Ly)

p_prime = 0
q_prime = 0

# initialize points for plotting
Hxplot=[]
Hyplot=[]
Ezplot=[]

# define final time
T = N * tou
for t_val in range(len(t)):
    
    for x_val in range(len(x)):
       
        for y_val in range(len(y)):
            # use equation (8) to find Jz(x,y,t)
            if (x_val == 0) or (x_val == Lx):
                Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega * t[t_val])

            elif (y_val == 0) or (x_val == Ly):
                Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega * t[t_val])

            else:
                Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega * t[t_val])

                # use equations 10 to define the matrices
                pval = x_val / ax
                qval = y_val / ax

                for pprime in range(0,P):
                    
                    Ez_sum = 0
                    Hx_sum = 0
                    Hy_sum = 0
                    Jz_sum = 0

                    for qprime in range(0,P):
                        
                        xprime = pprime * ax
                        yprime = qprime * ay

                        Ez_sum += E[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)
                        Hx_sum += X[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.cos(qval * qprime * np.pi / P)
                        Hy_sum += Y[pprime][qprime] * np.cos(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)
                        Jz_sum += J[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)

                    Ez[x_val][y_val] += Ez_sum  
                    Hx[x_val][y_val] += Hx_sum
                    Hy[x_val][y_val] += Hy_sum
                    Jz[x_val][y_val] += Jz_sum

            # take the fourier transforms of Hx, Hy, Jz, Ez 
            fourier_transform_Hx = dXXt2(Hx)
            fourier_transform_Hy = dXXt2(Hy)
            fourier_transform_Ez = dXXt2(Ez)
            fourier_transform_Jz = dXXt2(Jz)
            
            pval = x_val / ax
            qval = y_val / ax  
            
            # evolve the fourier coefficients using equations (11) which follow the Crank-Nicolson method
            E[x_val][y_val] = (((1 - (pval**2 * Dx**2) - (qval**2 * Dy**2)) * E[x_val - 1][y_val - 1]) + (2 * qval * Dy * X[x_val - 1][y_val - 1]) + ((2 * pval * Dx * Y[x_val - 1][y_val - 1]) + (tou * fourier_transform_Jz[x_val - 1][y_val - 1]))) / (1 + (pval**2 * Dx**2) + (qval**2 * Dy**2))
        
            X[x_val][y_val] = X[x_val - 1][y_val - 1] - (qval * Dy * (E[x_val][y_val] + E[x_val - 1][y_val - 1]))
            Y[x_val][y_val] = Y[x_val - 1][y_val - 1] - (pval * Dx * (E[x_val][y_val] + E[x_val - 1][y_val - 1]))
            
            # reconstruct Hx, Hy, Ez with the inverse Fourier transform
            fourier_transform_Hx = idXXt2(X, Hx)
            fourier_transform_Hy = idXXt2(Y, Hy)
            fourier_transform_Ez = idXXt2(E, Ez)
            fourier_transform_Jz = idXXt2(fourier_transform_Jz, Jz)
            
            # update for next time step
            Hx = fourier_transform_Hx
            Hy = fourier_transform_Hy
            Ez = fourier_transform_Ez
            Jz = fourier_transform_Jz
            
            # for plotting traces
            if x[x_val] == 0.5 and y[y_val] == 0:
                Hxplot.append(Hx[x_val][y_val])
                
            if y[y_val] == 0.5 and x[x_val] == 0:
                Hyplot.append(Hy[x_val][y_val])
            
            if y[y_val] == 0.5 and x[x_val] == 0.5:
                Ezplot.append(Ez[x_val][y_val]) 

# plots                
plt.plot(t, Hxplot)
plt.xlabel('Time (s)')
plt.ylabel('H$_x$(x = 0.5, y = 0)')
plt.title('H$_x$(x = 0.5, y = 0) Trace Over Time')
plt.show()

plt.plot(t, Hyplot)
plt.xlabel('Time (s)')
plt.ylabel('H$_y$(x = 0, y = 0.5)')
plt.title('H$_y$(x = 0, y = 0.5) Trace Over Time')
plt.show()

plt.plot(t, Ezplot)
plt.xlabel('Time (s)')
plt.ylabel('E$_z$(x = 0.5, y = 0.5)')
plt.title('E$_z$(x = 0.5, y = 0.5) Trace Over Time')
plt.show()

###############################################################################
# part d

N = 20#15 # should be 2000 # number of time steps
P = 32 # should be 32 # number of y and x steps
tou = 0.01 # duration of each time step

Lx = 1#P * ax the wall in the x direction
Ly = 1#P * ay the wall in the y direction
J0 = 1 # constant
ax = Lx/P  # size of the cells in the x direction
ay = Ly/P# size of the cells in the y direction


m = 1
n = 1
c = 1

# initialize lists of points for space and time
nlist = (np.arange(0, N))
p = (arange(0, P))
q = (np.arange(0, P))

# disrcitize t, x, and y
t = tou * nlist
x = (ax * p) 
y = (ay * q)

# initialize matrices
Hx = zeros((P, P), complex)
Hy = zeros((P, P), complex)
Ez = zeros((P, P), complex)
Jz = zeros((P, P), complex)

# initialize matrices to hold fourier coefficients
E = zeros((P, P), complex)
X = zeros((P, P), complex)
Y = zeros((P, P), complex)
J = zeros((P, P), complex)

Dx = np.pi * c * tou / (2 * Lx)
Dy = np.pi * c * tou / (2 * Ly)

p_prime = 0
q_prime = 0

# initialize points for plotting
Ezplot=[]
Ezmax = []

# define final time
T = N * tou
omega = np.arange(0, 9) # driving frequency

for k in range(len(omega)):
    for t_val in range(len(t)):
    
        for x_val in range(len(x)):
       
            for y_val in range(len(y)):
                # use equation (8) to find Jz(x,y,t)
                if (x_val == 0) or (x_val == Lx):
                    Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega[k] * t[t_val])

                elif (y_val == 0) or (x_val == Ly):
                    Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega[k] * t[t_val])

                else:
                    Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega[k] * t[t_val])

                    # use equations 10 to define the matrices
                    pval = x_val / ax
                    qval = y_val / ax

                    for pprime in range(0,P): 
                    
                        Ez_sum = 0
                        Hx_sum = 0
                        Hy_sum = 0
                        Jz_sum = 0

                        for qprime in range(0,P):
                        
                            xprime = pprime * ax
                            yprime = qprime * ay

                            Ez_sum += E[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)
                            Hx_sum += X[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.cos(qval * qprime * np.pi / P)
                            Hy_sum += Y[pprime][qprime] * np.cos(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)
                            Jz_sum += J[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)

                        Ez[x_val][y_val] += Ez_sum  
                        Hx[x_val][y_val] += Hx_sum
                        Hy[x_val][y_val] += Hy_sum
                        Jz[x_val][y_val] += Jz_sum

                # take the fourier transforms of Hx, Hy, Jz, Ez (use equation 10)
                fourier_transform_Hx = dXXt2(Hx)
                fourier_transform_Hy = dXXt2(Hy)
                fourier_transform_Ez = dXXt2(Ez)
                fourier_transform_Jz = dXXt2(Jz)
            
                pval = x_val / ax
                qval = y_val / ax  
            
                # evolve the fourier coefficients using equations (11) which follow the Crank-Nicolson method
                E[x_val][y_val] = (((1 - (pval**2 * Dx**2) - (qval**2 * Dy**2)) * E[x_val - 1][y_val - 1]) + (2 * qval * Dy * X[x_val - 1][y_val - 1]) + ((2 * pval * Dx * Y[x_val - 1][y_val - 1]) + (tou * fourier_transform_Jz[x_val - 1][y_val - 1]))) / (1 + (pval**2 * Dx**2) + (qval**2 * Dy**2))
        
                X[x_val][y_val] = X[x_val - 1][y_val - 1] - (qval * Dy * (E[x_val][y_val] + E[x_val - 1][y_val - 1]))
                Y[x_val][y_val] = Y[x_val - 1][y_val - 1] - (pval * Dx * (E[x_val][y_val] + E[x_val - 1][y_val - 1]))
            
                # reconstruct Hx, Hy, Ez with the inverse Fourier transform
                fourier_transform_Hx = idXXt2(X, Hx)
                fourier_transform_Hy = idXXt2(Y, Hy)
                fourier_transform_Ez = idXXt2(E, Ez)
                fourier_transform_Jz = idXXt2(fourier_transform_Jz, Jz)
            
                # update for next time step
                Hx = fourier_transform_Hx
                Hy = fourier_transform_Hy
                Ez = fourier_transform_Ez
                Jz = fourier_transform_Jz

                # for plotting Ez as a function of frequency
                if y[y_val] == 0.5 and x[x_val] == 0.5:
                    Ezplot.append(Ez[x_val][y_val]) 

        
    Ezmax.append(np.max(Ezplot))    
                
                
# plot                
plt.plot(omega, Ezmax)
plt.xlabel('omega')
plt.ylabel('E$_z$(x = 0.5, y = 0.5)')
plt.title('E$_z$(x = 0.5, y = 0.5) as a function of omega')
plt.show()

###############################################################################

# part e

N = 20#15 # should be 2000 # number of time steps
P = 32 # should be 32 # number of y and x steps
tou = 0.01 # duration of each time step

Lx = 1#P * ax the wall in the x direction
Ly = 1#P * ay the wall in the y direction
J0 = 1 # constant
ax = Lx/P  # size of the cells in the x direction
ay = Ly/P# size of the cells in the y direction


m = 1
n = 1
c = 1
omega = np.pi*c*np.sqrt((n*Lx)**-2 + (m*Ly)**-2) # driving frequency

# initialize lists of points for space and time
nlist = (np.arange(0, N))
p = (arange(0, P))
q = (np.arange(0, P))

# disrcitize t, x, and y
t = tou * nlist
x = (ax * p) 
y = (ay * q)

# initialize matrices
Hx = zeros((P, P), complex)
Hy = zeros((P, P), complex)
Ez = zeros((P, P), complex)
Jz = zeros((P, P), complex)

# initialize matrices to hold fourier coefficients
E = zeros((P, P), complex)
X = zeros((P, P), complex)
Y = zeros((P, P), complex)
J = zeros((P, P), complex)

Dx = np.pi * c * tou / (2 * Lx)
Dy = np.pi * c * tou / (2 * Ly)

p_prime = 0
q_prime = 0

# initialize points for plotting
Hxplot=[]
Hyplot=[]
Ezplot=[]

# define final time
T = N * tou
for t_val in range(len(t)):
    
    for x_val in range(len(x)):
       
        for y_val in range(len(y)):
            # use equation (8) to find Jz(x,y,t)
            if (x_val == 0) or (x_val == Lx):
                Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega * t[t_val])

            elif (y_val == 0) or (x_val == Ly):
                Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega * t[t_val])

            else:
                Jz[x_val][y_val] = J0 * np.sin((m * np.pi * x[x_val]) / Lx) * np.sin((n * np.pi * y[y_val]) / Ly) * np.sin(omega * t[t_val])

                # use equations 10 to define the matrices
                pval = x_val / ax
                qval = y_val / ax

                for pprime in range(1,P): # 
                    
                    Ez_sum = 0
                    Hx_sum = 0
                    Hy_sum = 0
                    Jz_sum = 0

                    for qprime in range(1,P):
                        
                        xprime = pprime * ax
                        yprime = qprime * ay

                        Ez_sum += E[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)
                        Hx_sum += X[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.cos(qval * qprime * np.pi / P)
                        Hy_sum += Y[pprime][qprime] * np.cos(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)
                        Jz_sum += J[pprime][qprime] * np.sin(pval * pprime * np.pi / P) * np.sin(qval * qprime * np.pi / P)

                    Ez[x_val][y_val] += Ez_sum  
                    Hx[x_val][y_val] += Hx_sum
                    Hy[x_val][y_val] += Hy_sum
                    Jz[x_val][y_val] += Jz_sum

            # take the fourier transforms of Hx, Hy, Jz, Ez (use equation 10)
            fourier_transform_Hx = dXXt2(Hx)
            fourier_transform_Hy = dXXt2(Hy)
            fourier_transform_Ez = dXXt2(Ez)
            fourier_transform_Jz = dXXt2(Jz)
            
            pval = x_val / ax
            qval = y_val / ax  
            
            # evolve the fourier coefficients using equations (11) which follow the Crank-Nicolson method
            E[x_val][y_val] = (((1 - (pval**2 * Dx**2) - (qval**2 * Dy**2)) * E[x_val - 1][y_val - 1]) + (2 * qval * Dy * X[x_val - 1][y_val - 1]) + ((2 * pval * Dx * Y[x_val - 1][y_val - 1]) + (tou * fourier_transform_Jz[x_val - 1][y_val - 1]))) / (1 + (pval**2 * Dx**2) + (qval**2 * Dy**2))
        
            X[x_val][y_val] = X[x_val - 1][y_val - 1] - (qval * Dy * (E[x_val][y_val] + E[x_val - 1][y_val - 1]))
            Y[x_val][y_val] = Y[x_val - 1][y_val - 1] - (pval * Dx * (E[x_val][y_val] + E[x_val - 1][y_val - 1]))
            
            # reconstruct Hx, Hy, Ez with the inverse Fourier transform
            fourier_transform_Hx = idXXt2(X, Hx)
            fourier_transform_Hy = idXXt2(Y, Hy)
            fourier_transform_Ez = idXXt2(E, Ez)
            fourier_transform_Jz = idXXt2(fourier_transform_Jz, Jz)
            
            # update for next time step
            Hx = fourier_transform_Hx
            Hy = fourier_transform_Hy
            Ez = fourier_transform_Ez
            Jz = fourier_transform_Jz
            
            # for plotting traces
            if x[x_val] == 0.5 and y[y_val] == 0:
                Hxplot.append(Hx[x_val][y_val])
                
            if y[y_val] == 0.5 and x[x_val] == 0:
                Hyplot.append(Hy[x_val][y_val])
            
            if y[y_val] == 0.5 and x[x_val] == 0.5:
                Ezplot.append(Ez[x_val][y_val]) 
  
# plot              
plt.plot(t, Hxplot)
plt.xlabel('Time (s)')
plt.ylabel('H$_x$(x = 0.5, y = 0)')
plt.title('H$_x$(x = 0.5, y = 0) Trace Over Time')
plt.show()

plt.plot(t, Hyplot)
plt.xlabel('Time (s)')
plt.ylabel('H$_y$(x = 0, y = 0.5)')
plt.title('H$_y$(x = 0, y = 0.5) Trace Over Time')
plt.show()

plt.plot(t, Ezplot)
plt.xlabel('Time (s)')
plt.ylabel('E$_z$(x = 0.5, y = 0.5)')
plt.title('E$_z$(x = 0.5, y = 0.5) Trace Over Time')
plt.show()
