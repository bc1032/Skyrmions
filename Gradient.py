import math
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize

import InitialiseSkyrmion
import SkyrmionEnergy

def derivativesphere():
    dx, dy = 1, 1
    dphi -= Sin(theta[x,y])*((b*Sin(theta[x,y])*Sin(phi[x,y])*(-theta(x,-1 + y) + theta(x,1 + y)))/dy + (b*Cos(phi[x,y])*Sin(theta[x,y])*(-theta[x-1,y] + theta[x+1,y]))/dx - \
                2*a*Cos(theta[x,y])*(((-theta(x,-1 + y) + theta(x,1 + y))*(-phi(x,-1 + y) + phi(x,1 + y)))/(4.*dy**2) + ((-theta[x-1,y] + theta[x+1,y])*(-phi[x-1,y] + phi[x+1,y]))/(4.*dx**2)) -\
                a*Sin(theta[x,y])*((phi(x,-1 + y) - 2*phi[x,y] + phi(x,1 + y))/dx + (phi[x-1,y] - 2*phi[x,y] + phi[x+1,y])/dx) + \
                a*Sin(2*phi[x,y])*(Cos(theta[x,y])*(-((theta(x,-1 + y) - 2*theta[x,y] + theta(x,1 + y))/dx) + (theta[x-1,y] - 2*theta[x,y] + theta[x+1,y])/dx + \
                ((-theta[x-1,y] + theta[x+1,y])*(-phi(x,-1 + y) + phi(x,1 + y)))/(2.*dx*dy) + ((-theta(x,-1 + y) + theta(x,1 + y))*(-phi[x-1,y] + phi[x+1,y]))/(2.*dx*dy)) + \
                Sin(theta[x,y])*((-theta(x,-1 + y) + theta(x,1 + y))**2/(4.*dy**2) - (-theta[x-1,y] + theta[x+1,y])**2/(4.*dx**2) + (-phi(x,-1 + y) + phi(x,1 + y))**2/(4.*dy**2) - \
                (-phi[x-1,y] + phi[x+1,y])**2/(4.*dx**2) + (phi(1,1 + y) - phi[x-1,y] - phi(x,-1 + y) + phi[x+1,y])/(2.*dx*dy))) + \
                a*Cos(2*phi[x,y])*(-2*Cos(theta[x,y])*((theta(1,1 + y) - theta[x-1,y] - theta(x,-1 + y) + theta[x+1,y])/(4.*dx*dy) + \
                ((-theta(x,-1 + y) + theta(x,1 + y))*(-phi(x,-1 + y) + phi(x,1 + y)))/(4.*dy**2) - ((-theta[x-1,y] + theta[x+1,y])*(-phi[x-1,y] + phi[x+1,y]))/(4.*dx**2)) + \
                Sin(theta[x,y])*(((-theta(x,-1 + y) + theta(x,1 + y))*(-theta[x-1,y] + theta[x+1,y]))/(2.*dx*dy) - (phi(x,-1 + y) - 2*phi[x,y] + phi(x,1 + y))/dx + \
                ((-phi(x,-1 + y) + phi(x,1 + y))*(-phi[x-1,y] + phi[x+1,y]))/(2.*dx*dy) + (phi[x-1,y] - 2*phi[x,y] + phi[x+1,y])/dx))

    dtheta -= -((b*Sin(theta[x,y])**2*Sin(phi[x,y])*(-phi(x,-1 + y) + phi(x,1 + y)))/dy) - (b*Cos(phi[x,y])*Sin(theta[x,y])**2*(-phi[x-1,y] + phi[x+1,y]))/dx - \
                a*Cos(theta[x,y])*Cos(2*phi[x,y])*(Cos(theta[x,y])*(-((theta(x,-1 + y) - 2*theta[x,y] + theta(x,1 + y))/dx) + (theta[x-1,y] - 2*theta[x,y] + theta[x+1,y])/dx + \
                ((-theta[x-1,y] + theta[x+1,y])*(-phi(x,-1 + y) + phi(x,1 + y)))/(2.*dx*dy) + ((-theta(x,-1 + y) + theta(x,1 + y))*(-phi[x-1,y] + phi[x+1,y]))/(2.*dx*dy)) + \
                Sin(theta[x,y])*((-theta(x,-1 + y) + theta(x,1 + y))**2/(4.*dy**2) - (-theta[x-1,y] + theta[x+1,y])**2/(4.*dx**2) + (-phi(x,-1 + y) + phi(x,1 + y))**2/(4.*dy**2) - \
                (-phi[x-1,y] + phi[x+1,y])**2/(4.*dx**2) + (phi(1,1 + y) - phi[x-1,y] - phi(x,-1 + y) + phi[x+1,y])/(2.*dx*dy))) + \
                a*Cos(theta[x,y])*Sin(2*phi[x,y])*(-2*Cos(theta[x,y])*((theta(1,1 + y) - theta[x-1,y] - theta(x,-1 + y) + theta[x+1,y])/(4.*dx*dy) + \
                ((-theta(x,-1 + y) + theta(x,1 + y))*(-phi(x,-1 + y) + phi(x,1 + y)))/(4.*dy**2) - ((-theta[x-1,y] + theta[x+1,y])*(-phi[x-1,y] + phi[x+1,y]))/(4.*dx**2)) + \
                Sin(theta[x,y])*(((-theta(x,-1 + y) + theta(x,1 + y))*(-theta[x-1,y] + theta[x+1,y]))/(2.*dx*dy) - (phi(x,-1 + y) - 2*phi[x,y] + phi(x,1 + y))/dx + \
                ((-phi(x,-1 + y) + phi(x,1 + y))*(-phi[x-1,y] + phi[x+1,y]))/(2.*dx*dy) + (phi[x-1,y] - 2*phi[x,y] + phi[x+1,y])/dx)) + \
                Cos(theta[x,y])*(-(a*Cos(theta[x,y])*((theta(x,-1 + y) - 2*theta[x,y] + theta(x,1 + y))/dx + (theta[x-1,y] - 2*theta[x,y] + theta[x+1,y])/dx)) + \
                Sin(theta[x,y])*(-2*c + a*((-theta(x,-1 + y) + theta(x,1 + y))**2/(4.*dy**2) + (-theta[x-1,y] + theta[x+1,y])**2/(4.*dx**2) + (-phi(x,-1 + y) + phi(x,1 + y))**2/(4.*dy**2) +\
                (-phi[x-1,y] + phi[x+1,y])**2/(4.*dx**2)))
