#!/usr/bin/env python3
from cyl2d import *
from critNoRefl import *

# TIME THE SOLVE:
import time

begintime = time.time()
height = 107.41995
nr = 100
nh = 100
alist, blist = makeAMatrix(height/2.0, height, nr, nh, 2, geommap)
dr = height/2.0 / nr
dh = height/ nh


# inverse power method function now requires dr, dh in order to locate group
# group transfer matrix as function of material at a location
eig, flux = inversePower(alist, blist, dr, dh, nr, geommap)
print(" k is {}".format(1./eig))
print("it worked")
endtime = time.time()
print("solve took {} seconds".format(endtime-begintime))

# plotting stuff 
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
#grid = fluxVector2Grid(flux[0], height/2.0, height, 100, 100)
#grid2 = fluxVector2Grid(flux[1], height/2.0, height, 100, 100)
x = np.linspace(0.0, height/2.0, nr)
y = np.linspace(0.0, height, nh)
xv, yv = np.meshgrid(x, y);
fig  = plt.figure()
ax = fig.add_subplot(1,2,1,projection='3d')
ax.scatter(xv, yv, flux[0],c=flux[0])
ax.set_title("Fast flux in unreflected reactor")
ax.set_xlabel("Radius (cm)")
ax.set_ylabel("Height (cm)")
ax.set_zlabel("Flux (arb)")
ax2 = fig.add_subplot(1,2,2,projection='3d')
ax2.scatter(xv, yv, flux[1],'b', c=flux[1])
ax2.set_title("Thermal flux in unreflected reactor")
ax2.set_xlabel("Radius (cm)")
ax2.set_ylabel("Height (cm)")
ax2.set_zlabel("Flux (arb)")
plt.show()
