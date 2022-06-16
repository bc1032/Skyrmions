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

def sphere(a,b,c,phi,theta,Lx,Ly):
    dx = 1
    dy = 1
    E=0.0
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

            dphidx1 = (phi[xr,y]-phi[xl,y])/(2.*dx)
            dphidxm = (phi[xr,y]-phi[xl,y] - 2*math.pi)/(2.*dx)
            dphidxp = (phi[xr,y]-phi[xl,y] + 2*math.pi)/(2.*dx)

            dphidx = min(dphidx1,dphidxm,dphidxp)

            dphidy1 = (phi[x,ya]-phi[x,yb])/(2.*dx)
            dphidym = (phi[x,ya]-phi[x,yb] - 2*math.pi)/(2.*dx)
            dphidyp = (phi[x,ya]-phi[x,yb] + 2*math.pi)/(2.*dx)

            dphidy = min(abs(dphidy1),abs(dphidym),abs(dphidyp))

            E += c*(-1 + math.cos(theta[x,y]))**2 + b*(math.cos(theta[x,y])*math.sin(theta[x,y])*(dphidx*math.cos(phi[x,y]) + dphidy*math.sin(phi[x,y])) - (math.cos(phi[x,y])*(theta[x,ya] - theta[x,yb]))/(2.*dy) + \
                (math.sin(phi[x,y])*(-theta[xl,y] + theta[xr,y]))/(2.*dx)) + a*(math.sin(phi[x,y])*(-(dphidx*math.sin(theta[x,y])) + (math.cos(theta[x,y])*(theta[x,ya] - theta[x,yb]))/(2.*dy)) + \
                math.cos(phi[x,y])*(dphidy*math.sin(theta[x,y]) + (math.cos(theta[x,y])*(-theta[xl,y] + theta[xr,y]))/(2.*dx)))**2
            #print(E)
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
# phi = np.loadtxt('phi.dat')
# theta = np.loadtxt('theta.dat')
# sphere(1,1,1, phi, theta)
