import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

import InitialiseSkyrmion
import SkyrmionEnergy

def derivativesphere(gamma):
    dx, dy = 1, 1
    tolerance = 1e5
    gamma = 1
    phi = np.loadtxt('phi.dat')
    theta = np.loadtxt('theta.dat')
    a,b,c = 1,1,1

    phi,theta,dphi, dtheta = minimise(gamma,phi,theta)
    SkyrmionEnergy.sphere(a,b,c)

    for i in range(0,1000):
        print(np.sum(dphi),np.sum(dtheta))
        print(i)
        SkyrmionEnergy.sphere(a,b,c)

        if (np.sum(dtheta)) < tolerance and np.sum(dtheta) < tolerance:
            np.savetxt('phifinal',phi)
            np.savetxt('thetafinal',theta)

            return(phi,theta)
        else:
            phi,theta,dphi, dtheta = minimise(gamma,phi,theta)

def minimise(gamma,phi,theta):
    a,b,c = 1,1,1
    Lx,Ly = 101, 101
    dphi = np.zeros((Lx,Ly))
    dtheta = np.zeros((Lx,Ly))
    dx,dy = 1,1
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

            dphi[x,y] = gamma*math.sin(theta[x,y])*((b*math.sin(theta[x,y])*math.sin(phi[x,y])*(-theta[x,yb] + theta[x,ya]))/dy \
                        + (b*math.cos(phi[x,y])*math.sin(theta[x,y])*(-theta[xl,y] + theta[xr,y]))/dx -\
                        2*a*math.cos(theta[x,y])*(((-theta[x,yb] + theta[x,ya])*(-phi[x,yb] + phi[x,ya]))/(4.*dy**2) + ((-theta[xl,y] + theta[xr,y])*(-phi[xl,y] + phi[xr,y]))/(4.*dx**2)) - \
                        a*math.sin(theta[x,y])*((phi[x,yb] - 2*phi[x,y] + phi[x,ya])/dx + (phi[xl,y] - 2*phi[x,y] + phi[xr,y])/dx) + \
                        a*math.cos(2*phi[x,y])*(-2*math.cos(theta[x,y])*((-theta[xl,yb] + theta[xr,ya])/(4.*dx*dy) + ((-theta[x,yb] + theta[x,ya])*(-phi[x,yb] + phi[x,ya]))/(4.*dy**2) - \
                        ((-theta[xl,y] + theta[xr,y])*(-phi[xl,y] + phi[xr,y]))/(4.*dx**2)) + \
                        math.sin(theta[x,y])*(((-theta[x,yb] + theta[x,ya])*(-theta[xl,y] + theta[xr,y]))/(2.*dx*dy) - (phi[x,yb] - 2*phi[x,y] + phi[x,ya])/dx + \
                        ((-phi[x,yb] + phi[x,ya])*(-phi[xl,y] + phi[xr,y]))/(2.*dx*dy) + (phi[xl,y] - 2*phi[x,y] + phi[xr,y])/dx)) + \
                        a*math.sin(2*phi[x,y])*(math.cos(theta[x,y])*(-((theta[x,yb] - 2*theta[x,y] + theta[x,ya])/dx) + (theta[xl,y] - 2*theta[x,y] + theta[xr,y])/dx + \
                        ((-theta[xl,y] + theta[xr,y])*(-phi[x,yb] + phi[x,ya]))/(2.*dx*dy) + ((-theta[x,yb] + theta[x,ya])*(-phi[xl,y] + phi[xr,y]))/(2.*dx*dy)) + \
                        math.sin(theta[x,y])*((-theta[x,yb] + theta[x,ya])**2/(4.*dy**2) - (-theta[xl,y] + theta[xr,y])**2/(4.*dx**2) + (-phi[x,yb] + phi[x,ya])**2/(4.*dy**2) - \
                        (-phi[xl,y] + phi[xr,y])**2/(4.*dx**2) + (-phi[xl,yb] + phi[xr,ya])/(2.*dx*dy))))

            dtheta[x,y] = -gamma*((b*math.sin(theta[x,y])**2*math.sin(phi[x,y])*(-phi[x,yb] + phi[x,ya]))/dy)\
                        - (b*math.cos(phi[x,y])*math.sin(theta[x,y])**2*(-phi[xl,y] + phi[xr,y]))/dx + \
                        a*math.cos(theta[x,y])*math.sin(2*phi[x,y])*(-2*math.cos(theta[x,y])*((-theta[xl,yb] + theta[xr,ya])/(4.*dx*dy) + \
                        ((-theta[x,yb] + theta[x,ya])*(-phi[x,yb] + phi[x,ya]))/(4.*dy**2) - ((-theta[xl,y] + theta[xr,y])*(-phi[xl,y] + phi[xr,y]))/(4.*dx**2)) + \
                        math.sin(theta[x,y])*(((-theta[x,yb] + theta[x,ya])*(-theta[xl,y] + theta[xr,y]))/(2.*dx*dy) - (phi[x,yb] - 2*phi[x,y] + phi[x,ya])/dx + \
                        ((-phi[x,yb] + phi[x,ya])*(-phi[xl,y] + phi[xr,y]))/(2.*dx*dy) + (phi[xl,y] - 2*phi[x,y] + phi[xr,y])/dx)) + \
                        math.cos(theta[x,y])*(-(a*math.cos(theta[x,y])*((theta[x,yb] - 2*theta[x,y] + theta[x,ya])/dx + (theta[xl,y] - 2*theta[x,y] + theta[xr,y])/dx)) + \
                        math.sin(theta[x,y])*(-2*c + a*((-theta[x,yb] + theta[x,ya])**2/(4.*dy**2) + (-theta[xl,y] + theta[xr,y])**2/(4.*dx**2) + (-phi[x,yb] + phi[x,ya])**2/(4.*dy**2) + \
                        (-phi[xl,y] + phi[xr,y])**2/(4.*dx**2)))) - a*math.cos(theta[x,y])*math.cos(2*phi[x,y])*\
                        (math.cos(theta[x,y])*(-((theta[x,yb] - 2*theta[x,y] + theta[x,ya])/dx) + (theta[xl,y] - 2*theta[x,y] + theta[xr,y])/dx + \
                        ((-theta[xl,y] + theta[xr,y])*(-phi[x,yb] + phi[x,ya]))/(2.*dx*dy) + ((-theta[x,yb] + theta[x,ya])*(-phi[xl,y] + phi[xr,y]))/(2.*dx*dy)) + \
                        math.sin(theta[x,y])*((-theta[x,yb] + theta[x,ya])**2/(4.*dy**2) - (-theta[xl,y] + theta[xr,y])**2/(4.*dx**2) + (-phi[x,yb] + phi[x,ya])**2/(4.*dy**2) - \
                        (-phi[xl,y] + phi[xr,y])**2/(4.*dx**2) + (-phi[xl,yb] + phi[xr,ya])/(2.*dx*dy)))
    phi -= dphi
    theta -= dtheta
    np.savetxt('dtheta.dat', dtheta)
    np.savetxt('dphi.dat', dphi)

    print(dtheta)
    #phi = phi - dphi
    #theta = theta - dtheta
    return(phi,theta, dphi,dtheta)
phi,theta = derivativesphere(1)
file = open('finaldirectorfield.dat', 'w')

Lx,Ly = 101, 101

for x in range(0,Lx):
    for y in range(0,Ly):
        i = math.sin(theta[x,y])*math.cos(phi[x,y])
        j = math.sin(theta[x,y])*math.sin(phi[x,y])
        k = math.cos(theta[x,y])
        file.write("%f  %f  %f\n" % (i,j,k))
file.close()
