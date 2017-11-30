#!/usr/bin/env python3
# a simple example of using Cyl2D in a 6-group setting
# these group constants are definitely messed up though.
from cyl2d import *

nr = 100
nz = 100

# material, takes diff coeff, sig r, nu sig f, group transfer matrix
# note that sig1->2 is its removal minus absorbtion XS
gtr = 'gtrs.txt'
gtr = np.loadtxt(gtr)
gtr = np.reshape(gtr, [6,6])
LWRMaterial = Material([0.0280064,0.0184021,0.011311,0.0144786,0.013975,0.0128252],
                       [16.5512,21.7253,31.8009,24.2093,25.0351,27.2159],
                       [0.31780612,0.27663047,0.37039145,0.6277296,1.30504206,3.52626536],
                       gtr)


# geommap must be defined for all cyl2d runs
def geommap(x, y):
    return LWRMaterial

Alist, Blist = makeAMatrix(2.0, 2.0, nr, nz, 2, geommap)
dr = 2.0/ nr
dh = 2.0 / nz
eig, vec = inversePower(Alist, Blist, dr, dh, nr, geommap)
k = 1.0/eig
print("k=",k)
