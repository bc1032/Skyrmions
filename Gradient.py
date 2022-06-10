import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

import InitialiseSkyrmion
import SkyrmionEnergy
Lx, Ly = 151, 151
R = 40
numtwists = 5.0
twist = 1.0/numtwists
dx,dy = 1,1
numits = 10000

tolerance = 5e-3
gamma = 1e-4
a,b,c = 5,1,1e-2

def derivativesphere(gamma):
    i=0

    phi = np.loadtxt('phi.dat')
    theta = np.loadtxt('theta.dat')
    E = SkyrmionEnergy.sphere(a,b,c, phi,theta,Lx,Ly)
    phi,theta = minimise(gamma,phi,theta,a,b,c,i)
    #SkyrmionEnergy.sphere(a,b,c, phi,theta)

    for i in range(0,numits):
        print(i)
        Eprev = E
        #phi,theta,dphi, dtheta = minimise(gamma,phi,theta)

        E = SkyrmionEnergy.sphere(a,b,c, phi, theta,Lx,Ly)

        # if (np.sum(dtheta)) < tolerance and np.sum(dtheta) < tolerance:
        #     np.savetxt('phifinal',phi)
        #     np.savetxt('thetafinal',theta)
        #     return(phi,theta)
        print(-(Eprev - E))
        if abs(Eprev - E) < tolerance:
            np.savetxt('phifinal',phi)
            np.savetxt('thetafinal',theta)
            return(phi,theta)
        else:
            phi,theta = minimise(gamma,phi,theta,a,b,c,i)

    np.savetxt('phifinal',phi)
    np.savetxt('thetafinal',theta)
    return(phi,theta)

def minimise(gamma,phi,theta,a,b,c,i):
    dphi = np.zeros((Lx,Ly))
    dtheta = np.zeros((Lx,Ly))
    #gamma = -1e-3
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

            dphidx = min(abs(dphidx1),abs(dphidxm),abs(dphidxp))

            dphidy1 = (phi[x,ya]-phi[x,yb])/(2.*dx)
            dphidym = (phi[x,ya]-phi[x,yb] - 2*math.pi)/(2.*dx)
            dphidyp = (phi[x,ya]-phi[x,yb] + 2*math.pi)/(2.*dx)

            dphidy = min(abs(dphidy1),abs(dphidym),abs(dphidyp))

            dphidxdy1 = (phi[xr,ya] - phi[xl,yb])/(4.*dx*dy)
            dphidxdyp = (phi[xr,ya] - phi[xl,yb] + 2*math.pi)/(4.*dx*dy)
            dphidxdym = (phi[xr,ya] - phi[xl,yb] - 2*math.pi)/(4.*dx*dy)

            dphidxdy = min(abs(dphidxdy1),abs(dphidxdym),abs(dphidxdyp))

            dphidxsq1 = (phi[xr,y]+phi[xl,y] - 2*phi[x,y])/(2.*dx)
            dphidxsqm = (phi[xr,y]+phi[xl,y] - 2*phi[x,y] - 2*math.pi)/(2.*dx)
            dphidxsqp = (phi[xr,y]+phi[xl,y] - 2*phi[x,y] + 2*math.pi)/(2.*dx)

            dphidxsq = min(abs(dphidxsq1),abs(dphidxsqm),abs(dphidxsqp))

            dphidysq1 = (phi[x,ya]+phi[x,yb] - 2*phi[x,y])/(2.*dx)
            dphidysqm = (phi[x,ya]+phi[x,yb] - 2*phi[x,y] - 2*math.pi)/(2.*dx)
            dphidysqp = (phi[x,ya]+phi[x,yb] - 2*phi[x,y] + 2*math.pi)/(2.*dx)

            dphidysq = min(abs(dphidysq1),abs(dphidysqm),abs(dphidysqp))

            dphi[x,y] = gamma*(math.sin(theta[x,y])*(-(a*(dphidxsq + dphidysq)*math.sin(theta[x,y])) + (b*math.sin(theta[x,y])*math.sin(phi[x,y])*(theta[x,ya] - theta[x,yb]))/dy + \
                        (b*math.cos(phi[x,y])*math.sin(theta[x,y])*(-theta[xl,y] + theta[xr,y]))/dx - 2*a*math.cos(theta[x,y])*((dphidy*(theta[x,ya] - theta[x,yb]))/(2.*dy) + (dphidx*(-theta[xl,y] + theta[xr,y]))/(2.*dx)) + \
                        a*math.sin(2*phi[x,y])*(math.sin(theta[x,y])*(-dphidx**2 + 2*dphidxdy + dphidy**2 + (theta[x,ya] - theta[x,yb])**2/(4.*dy**2) - (-theta[xl,y] + theta[xr,y])**2/(4.*dx**2)) + \
                        math.cos(theta[x,y])*((dphidx*(theta[x,ya] - theta[x,yb]))/dy - (-2*theta[x,y] + theta[x,ya] + theta[x,yb])/dx + (dphidy*(-theta[xl,y] + theta[xr,y]))/dx + (-2*theta[x,y] + theta[xl,y] + theta[xr,y])/dx)) + \
                        a*math.cos(2*phi[x,y])*(math.sin(theta[x,y])*(dphidxsq + 2*dphidx*dphidy - dphidysq + ((theta[x,ya] - theta[x,yb])*(-theta[xl,y] + theta[xr,y]))/(2.*dx*dy)) - \
                        2*math.cos(theta[x,y])*((dphidy*(theta[x,ya] - theta[x,yb]))/(2.*dy) - (dphidx*(-theta[xl,y] + theta[xr,y]))/(2.*dx) + (-theta[xl,yb] + theta[xr,ya])/(4.*dx*dy)))))

            dtheta[x,y] = gamma*(-2*b*dphidx*math.cos(phi[x,y])*math.sin(theta[x,y])**2 - 2*b*dphidy*math.sin(theta[x,y])**2*math.sin(phi[x,y]) -
                        a*math.cos(theta[x,y])*math.cos(2*phi[x,y])*(math.sin(theta[x,y])*(-dphidx**2 + 2*dphidxdy + dphidy**2 + (theta[x,ya] - theta[x,yb])**2/(4.*dy**2) - (-theta[xl,y] + theta[xr,y])**2/(4.*dx**2)) + \
                        math.cos(theta[x,y])*((dphidx*(theta[x,ya] - theta[x,yb]))/dy - (-2*theta[x,y] + theta[x,ya] + theta[x,yb])/dx + (dphidy*(-theta[xl,y] + theta[xr,y]))/dx + (-2*theta[x,y] + theta[xl,y] + theta[xr,y])/dx)) +\
                        math.cos(theta[x,y])*(-(a*math.cos(theta[x,y])*((-2*theta[x,y] + theta[x,ya] + theta[x,yb])/dx + (-2*theta[x,y] + theta[xl,y] + theta[xr,y])/dx)) +\
                        math.sin(theta[x,y])*(-2*c + a*(dphidx**2 + dphidy**2 + (theta[x,ya] - theta[x,yb])**2/(4.*dy**2) + (-theta[xl,y] + theta[xr,y])**2/(4.*dx**2)))) +\
                        a*math.cos(theta[x,y])*math.sin(2*phi[x,y])*(math.sin(theta[x,y])*(dphidxsq + 2*dphidx*dphidy - dphidysq + ((theta[x,ya] - theta[x,yb])*(-theta[xl,y] + theta[xr,y]))/(2.*dx*dy)) -\
                        2*math.cos(theta[x,y])*((dphidy*(theta[x,ya] - theta[x,yb]))/(2.*dy) - (dphidx*(-theta[xl,y] + theta[xr,y]))/(2.*dx) + (-theta[xl,yb] + theta[xr,ya])/(4.*dx*dy))))

    phi -= dphi
    theta -= dtheta
    file = open('finaldirectorfield%06d.dat' % (i), 'w')
    for x in range(0,Lx):
        for y in range(0,Ly):
            i = math.sin(theta[x,y])*math.cos(phi[x,y])
            j = math.sin(theta[x,y])*math.sin(phi[x,y])
            k = math.cos(theta[x,y])
            file.write("%f  %f  %f\n" % (i,j,k))
    file.close()

    np.savetxt('dtheta.dat', dtheta)
    np.savetxt('dphi.dat', dphi)
    np.savetxt('theta.dat', theta)
    np.savetxt('phi.dat', phi)
    ##print(dtheta)
    #print(dphi)
    #phi = phi - dphi
    #theta = theta - dtheta
    return(phi,theta)
InitialiseSkyrmion.initialise(Lx,Ly,R,twist)

phi,theta = derivativesphere(gamma)
file = open('finaldirectorfield.dat', 'w')

for x in range(0,Lx):
    for y in range(0,Ly):
        i = math.sin(theta[x,y])*math.cos(phi[x,y])
        j = math.sin(theta[x,y])*math.sin(phi[x,y])
        k = math.cos(theta[x,y])
        file.write("%f  %f  %f\n" % (i,j,k))
file.close()
