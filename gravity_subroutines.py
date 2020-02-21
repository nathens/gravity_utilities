# gravity_subroutines.py
# Created: February 21st, 2020

"""
These functions are python implementations of fortran subroutines found in:  Blakely, R.J., 1996. Potential Theory in Gravity and Magnetic Applications. Potential Theory in Gravity and Magnetic Applications, by Richard J. Blakely, pp. 461. ISBN 0521575478. Cambridge, UK: Cambridge University Press, September 1996.

Refer to the above citation for details. 

"""

gitimport sys 
import numpy as np 
from scipy.integrate import tplquad

def gbox(x0, y0, z0, x1, y1, z1, x2, y2, z2, rho):
    gamma = 6.670e-11 
    twopi = 2 * np.pi 
    km2m = 1e3
    si2mg = 1e5 
    x = [x0 - x1, x0 - x2]
    y = [y0 - y1, y0 - y2]
    z = [z0 - z1, z0 - z2]
    res = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                rijk = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)
                ijk = (-1)**(i+1) * (-1)**(j+1) * (-1)**(k+1)
                arg1 = np.arctan2(x[i] * y[j], z[k] * rijk)
                if arg1 < 0: 
                    arg1 += twopi
                arg2 = rijk + y[j]
                arg3 = rijk + x[i]
                arg2 = np.log(arg2)
                arg3 = np.log(arg3)
                res += ijk * (z[k] * arg1 - x[i] * arg2 - y[j] * arg3)
    g = rho * gamma * res * si2mg * km2m
    return g

def gbox_by_integration(x0, y0, z0, x1, y1, z1, x2, y2, z2, rho):
    gamma = 6.670e-11 
    km2m = 1e3
    si2mg = 1e5 
    f = lambda z, y, x: (z0-z) / ((x0 - x)**2 + (y0-y)**2 + (z0 - z)**2)**(3/2)
    arg1 = tplquad(f, x1, x2, lambda x: y1, lambda x: y2, lambda x, y: z1, lambda x, y: z2)[0]
    g = -1 * rho * gamma * arg1 * si2mg * km2m
    return g

def gpoly(x0, z0, xcorn, zcorn, rho):
    gamma = 6.670e-11 
    km2m = 1e3
    si2mg = 1e5
    xcorn.append(xcorn[0])
    zcorn.append(zcorn[0])
    res = 0
    for i in range(len(xcorn) - 1):
        x1 = xcorn[i] - x0 
        z1 = zcorn[i] - z0
        x2 = xcorn[i+1] - x0 
        z2 = zcorn[i+1] - z0
        r1sq = x1**2 + z1**2 
        r2sq = x2**2 + z2**2
        denom = z2 - z1 
        if denom == 0:
            denom = 1e-6 
        alpha = (x2 - x1) / denom 
        beta = (x1 * z2 - x2 * z1) / denom 
        factor = beta / (1 + alpha**2)
        term1 = 0.5 * (np.log(r2sq) - np.log(r1sq))
        term2 = np.arctan2(z2, x2) - np.arctan2(z1, x1)
        res += factor * (term1 - alpha * term2)
    g = 2 * rho * gamma * res * si2mg * km2m
    return g

