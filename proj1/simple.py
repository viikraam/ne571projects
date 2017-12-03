#!/usr/bin/env python3
#a basic use of cyl2d, crit search without reflector.
# Finds the critical reactor dimension with no reflect, and r=h/2.
# uses secant method.
from cyl2d import *
from xs import *

nr = 100
nz = 100

# material, takes diff coeff, sig r, nu sig f, group transfer matrix
# note that sig1->2 is its removal minus absorbtion XS

void, diffC, sigA, nSigF, gtrans = mat_val(0)

LWRMaterial = Material(void, diffC, sigA, nSigF, gtrans)

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
    dr = height/2.0 / nr
    dh = height / nz
    eig, vec = inversePower(Alist, Blist, dr, dh, nr, geommap)
    k = 1.0/eig
    rho = (k-1.0)/k
    return eig

# if running from the command line, do this:
y1 = reactorRho(100.0, nr, nz, nr)

print(y1)