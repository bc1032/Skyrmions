import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

import InitialiseSkyrmion
import SkyrmionEnergy
Lx, Ly = 151, 151
R = 30
numtwists = 1.0
twist = 1.0/numtwists
dx,dy = 1,1
numits = 10000

tolerance = -1e-3
gamma = 1e-2
a,b,c = 10,0,1#c is negative

def derivativesphere(gamma):
    i=0
    Energy = []
    phi = np.loadtxt('phi.dat')
    theta = np.loadtxt('theta.dat')
    E = SkyrmionEnergy.sphere(a,b,c, phi,theta,Lx,Ly)
    phi,theta = minimise(gamma,phi,theta,a,b,c,i)
    #SkyrmionEnergy.sphere(a,b,c, phi,theta)
    n=0
    Energy.append(E)
    for i in range(0,numits):
        print(i)
        Eprev = E
        E = SkyrmionEnergy.sphere(a,b,c, phi, theta,Lx,Ly)
        Energy.append(E)
        print(-(Eprev - E))
        # if E-Eprev > -50 and n == 0:
        #     gamma /= 10
        #     n+=1
        # elif E-Eprev > -4:
        #     gamma /= 10
        #     n+=1
        # elif E-Eprev > -1:
        #     gamma *= 10
        if (E - Eprev) > tolerance:
            np.savetxt('phifinal',phi)
            np.savetxt('thetafinal',theta)
            return(phi,theta)
        else:
            phi,theta = minimise(gamma,phi,theta,a,b,c,i)
            phiprev, thetaprev = phi,theta

    Energy = np.asarray(Energy)
    np.savetxt('Energyevolution.dat',Energy)
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

            dthetadx = (theta[xr,y]-theta[xl,y])/(2.*dx)
            dthetady = (theta[x,ya]-theta[x,yb])/(2.*dx)
            dthetadxsq = (theta[xr,y]+theta[xl,y] - 2.*theta[x,y])/(dx)
            dthetadysq = (theta[x,ya]+theta[x,yb] - 2.*theta[x,y])/(dx)

            dphidx1 = (phi[xr,y]-phi[xl,y])/(2.*dx)
            dphidxm = (phi[xr,y]-phi[xl,y] - 2*math.pi)/(2.*dx)
            dphidxp = (phi[xr,y]-phi[xl,y] + 2*math.pi)/(2.*dx)

            if abs(dphidx1) < abs(dphidxm) and abs(dphidx1) < abs(dphidxp):
                dphidx = dphidx1
            elif abs(dphidxm) < abs(dphidx1) and abs(dphidxm) < abs(dphidxp):
                dphidx = dphidxm
            else:
                dphidx = dphidxp
            #dphidx = min(abs(dphidx1),abs(dphidxm),abs(dphidxp))

            dphidy1 = (phi[x,ya]-phi[x,yb])/(2.*dx)
            dphidym = (phi[x,ya]-phi[x,yb] - 2*math.pi)/(2.*dx)
            dphidyp = (phi[x,ya]-phi[x,yb] + 2*math.pi)/(2.*dx)

            if abs(dphidy1) < abs(dphidym) and abs(dphidy1) < abs(dphidyp):
                dphidy = dphidy1
            elif abs(dphidym) < abs(dphidy1) and abs(dphidym) < abs(dphidyp):
                dphidy = dphidym
            else:
                dphidy = dphidyp
            #dphidy = min(abs(dphidy1),abs(dphidym),abs(dphidyp))

            dphidxdy1 = (phi[xr,ya] + phi[xl,yb] - phi[xl,ya] - phi[xr,yb])/(4*dx*dy)
            dphidxdyp = (phi[xr,ya] + phi[xl,yb] - phi[xl,ya] - phi[xr,yb] + 2*math.pi)/(4*dx*dy)
            dphidxdym = (phi[xr,ya] + phi[xl,yb] - phi[xl,ya] - phi[xr,yb] - 2*math.pi)/(4*dx*dy)

            if abs(dphidxdy1) < abs(dphidxdym) and abs(dphidxdy1) < abs(dphidxdyp):
                dphidxdy = dphidxdy1
            elif abs(dphidxdym) < abs(dphidxdy1) and abs(dphidxdym) < abs(dphidxdyp):
                dphidxdy = dphidxdym
            else:
                dphidxdy = dphidxdyp

            #dphidxdy = min(abs(dphidxdy1),abs(dphidxdym),abs(dphidxdyp))

            dphidxsq1 = (phi[xr,y]+phi[xl,y] - 2*phi[x,y])/(dx*dx)
            dphidxsqm = (phi[xr,y]+phi[xl,y] - 2*phi[x,y] - 2*math.pi)/(dx*dx)
            dphidxsqp = (phi[xr,y]+phi[xl,y] - 2*phi[x,y] + 2*math.pi)/(dx*dx)

            if abs(dphidxsq1) < abs(dphidxsqm) and abs(dphidxsq1) < abs(dphidxsqp):
                dphidxsq = dphidxsq1
            elif abs(dphidxsqm) < abs(dphidxsq1) and abs(dphidxsqm) < abs(dphidxsqp):
                dphidxsq = dphidxsqm
            else:
                dphidxsq = dphidxsqp

            #dphidxsq = min(abs(dphidxsq1),abs(dphidxsqm),abs(dphidxsqp))

            dphidysq1 = (phi[x,ya]+phi[x,yb] - 2*phi[x,y])/(dy*dy)
            dphidysqm = (phi[x,ya]+phi[x,yb] - 2*phi[x,y] - 2*math.pi)/(dy*dy)
            dphidysqp = (phi[x,ya]+phi[x,yb] - 2*phi[x,y] + 2*math.pi)/(dy*dy)

            if abs(dphidysq1) < abs(dphidysqm) and abs(dphidysq1) < abs(dphidysqp):
                dphidysq = dphidysq1
            elif abs(dphidysqm) < abs(dphidysq1) and abs(dphidysqm) < abs(dphidysqp):
                dphidysq = dphidysqm
            else:
                dphidysq = dphidysqp

            #dphidysq = min(abs(dphidysq1),abs(dphidysqm),abs(dphidysqp))

            dphi[x,y] = gamma*(2*math.sin(theta[x,y])*(-2*a*(dphidx*dthetadx + dphidy*dthetady)*math.cos(theta[x,y]) - \
                        a*(dphidxsq + dphidysq)*math.sin(theta[x,y]) + b*dthetadx*math.cos(phi[x,y])*math.sin(theta[x,y]) + \
                        b*dthetady*math.sin(theta[x,y])*math.sin(phi[x,y])))

            dtheta[x,y] = gamma*(-2*a*(dthetadxsq + dthetadysq) - 2*b*dphidx*math.cos(phi[x,y])*math.sin(theta[x,y])**2 + 4*c*math.cos(theta[x,y])*math.sin(theta[x,y])**3 +\
                        a*(dphidx**2 + dphidy**2)*math.sin(2*theta[x,y]) - 2*b*dphidy*math.sin(theta[x,y])**2*math.sin(phi[x,y]))
            # dphi[x,y] = 0
            # dtheta[x,y]= gamma*(-c*math.sin(2.0*theta[x,y]))

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
