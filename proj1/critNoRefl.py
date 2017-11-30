#!/usr/bin/env python3
#a basic use of cyl2d, crit search without reflector.
# Finds the critical reactor dimension with no reflect, and r=h/2.
# uses secant method.
from cyl2d import *

nr = 100
nz = 100

# material, takes diff coeff, sig r, nu sig f, group transfer matrix
# note that sig1->2 is its removal minus absorbtion XS
LWRMaterial = Material([900.0], [1.2627, .3543], [.02629, .1210],
        [.008476, .18514], [ [0.0,0.0], [.02619-.01207, 0.0] ])

# geommap must be defined for all cyl2d runs
def geommap(x, y):
    """ For multiregion cores, this maps r,z to
    a Material object that contains cross sections. """
    return LWRMaterial

def reactorRho(height, nr, nh, numi):
    """Returns reactor reactivity as a function of reactor
    height under the constraint R=H/2. It would be nice to
    keep dr or dh under some value to ensure solutions are
    O(something), but the r term kind of complicates this issue."""
    Alist, Blist = makeAMatrix(height/2.0, height, nr, nz, 2, geommap)
    # dbg dbg
    import numpy as np
    np.savetxt('diffusive2DRadial.txt', Alist[0])
    quit()
    # end dbg
    dr = height/2.0 / nr
    dh = height / nz
    eig, vec = inversePower(Alist, Blist, dr, dh, nr, geommap)
    k = 1.0/eig
    rho = (k-1.0)/k
    return rho

# if running from the command line, do this:
if __name__=='__main__':
    x0 = 50.0
    x1 = 100.0
    y1 = 1e6

    # secant method loop:
    while np.abs(y1) > 1e-6:
        y0 = reactorRho(x0, nr, nz, nr)
        y1 = reactorRho(x1, nr, nz, nr)
        newx = secant(x0, x1, y0, y1)
        print('reactivity=',y1)
        print('New reactor height is h={} cm'.format(newx))
        x0 = x1
        x1 = newx

    print("FINAL HEIGHT:", newx, "cm")

