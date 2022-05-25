import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

import InitialiseSkyrmion

# if not os.path.exists('results'):
#     os.makedirs('results')
#
# os.chdir('results')


Lx, Ly = 101, 101
R = 40

def sphere(a,b,c):
    dx = 1
    dy = 1
    E=0.0
    phi = np.loadtxt('phi.dat')
    theta = np.loadtxt('theta.dat')

    for x in range(0,Lx):

        if x == 0:
                xl = Lx-1
                xr = 1
        elif x == Lx-1:
            xl = Lx-2
            xr = 0
        else:
            xl = x - 1
            xr = x + 1

        for y in range(0,Ly):
            if y == 0:
                ya = y+1
                yb = Ly-1
            elif y == Ly-1:
                ya = 0
                yb = Ly-2
            else:
                ya = y+1
                yb = y-1

            E += c*(-1 + math.cos(theta[x,y])**2) + b*(-(math.cos(phi[x,y])*(-theta[x,yb] + theta[x,ya]))/(2.*dy) + (math.sin(phi[x,y])*(-theta[xl,y] + theta[xr,y]))/(2.*dx) +\
                math.cos(theta[x,y])*math.sin(theta[x,y])*((math.sin(phi[x,y])*(-phi[x,yb] + phi[x,ya]))/(2.*dy) + (math.cos(phi[x,y])*(-phi[xl,y] + phi[xr,y]))/(2.*dx))) + \
                a*(math.cos(phi[x,y])*((math.cos(theta[x,y])*(-theta[xl,y] + theta[xr,y]))/(2.*dx) + (math.sin(theta[x,y])*(-phi[x,yb] + phi[x,ya]))/(2.*dy)) + \
                math.sin(phi[x,y])*((math.cos(theta[x,y])*(-theta[x,yb] + theta[x,ya]))/(2.*dy) - (math.sin(theta[x,y])*(-phi[xl,y] + phi[xr,y]))/(2.*dx)))**2
    print(E)
    return(E)








def projection(u,v):

    E = (c*(1 - u[x,y]**2 - v[x,y]**2)**2)/(1 + u[x,y]**2 + v[x,y]**2)**2 +\
        (b*(-(((-u[x,yb] + u[x,ya])*(1 + u[x,y]**2 - v[x,y]**2))/dy) + 4*u[x,y]*v[x,y]*((-u[xl,y] + u[xr,y])/(2.*dx) - (-v[x,yb] + v[x,ya])/(2.*dy)) +\
        ((1 - u[x,y]**2 + v[x,y]**2)*(-v[xl,y] + v[xr,y]))/dx))/(1 + u[x,y]**2 + v[x,y]**2)**2 +\
        (4*a*(((-u[xl,y] + u[xr,y])*(1 + v[x,y]**2))/(2.*dx) - ((-1 + v[x,y]**2)*(-v[x,yb] + v[x,ya]))/(2.*dy) +\
        u[x,y]**2*(-(-u[xl,y] + u[xr,y])/(2.*dx) + (-v[x,yb] + v[x,ya])/(2.*dy)) -\
        2*u[x,y]*v[x,y]*((-u[x,yb] + u[x,ya])/(2.*dy) + (-v[xl,y] + v[xr,y])/(2.*dx)))**2)/(1 + u[x,y]**2 + v[x,y]**2)**4

    return(E)

sphere(1,1,1)
