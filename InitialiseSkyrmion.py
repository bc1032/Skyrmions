import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

if not os.path.exists('results'):
    os.makedirs('results')

os.chdir('results')


Lx, Ly = 51, 51
R = 20

def initialise(Lx, Ly, R):
    theta = np.zeros([Lx, Ly])
    phi = np.zeros([Lx,Ly])
    centrex = int(Lx/2)
    centrey = int(Ly/2)


    for x in range(0,Lx):
        for y in range(0,Ly):
            dx = abs(centrex - x)
            dy = abs(centrey - y)

            #If inside circle of Radius, R, change elements to desired quantities.
            if dx > R or dy > R:
                phi[x,y] = math.pi/2
            elif dx + dy <= R:
                phi[x,y] = math.atan2(y,x)
                theta[x,y] =  math.pi*math.sqrt(x**2 + y**2)/(2*R)
            elif math.sqrt(dx**2 + dy**2) < R:
                phi[x,y] = math.atan2(y,x)
                theta[x,y] =  math.pi*math.sqrt(x**2 + y**2)/(2*R)
            else:
                phi[x,y] = math.pi/2

    print(phi,theta)
    np.savetxt('phi.dat', phi)
    np.savetxt('theta.dat', theta)
    return(0)

initialise(Lx,Ly,R)
