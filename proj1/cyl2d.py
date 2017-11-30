#!/usr/bin/env python3
# NE 571 Project 2
# 2 speed diffusion in multiplying cylinder:
# perform eigenvalue search on n-group, arbitrary material
# geometry.
# Units are CGS.
import numpy as np
import scipy.linalg as la
import scipy.interpolate as sci # XS interp
CSpline = sci.CubicSpline

debug = False

class Material(object):
    """ base class for materials. Contains methods for computing material properties
    in anticipation of implementing temperature implementation."""

    def __init__(self, tempvec, diffC, sigA, nSigF, gtrans):
        """ initializing a material. This means going
        ahead and forming splines
        
        Should pass in a vector of temperatures, followed by
        corresponding vector of diff. coefficients, sigA, nSigF,
        GTRANSXS each containing vectors of length <numgroups>"""

        # input checking
        if type(diffC)!=list or type(sigA)!=list or type(nSigF)!=list:
            raise Exception("please provide group constants in the form of a list,"
                            "descending from fastest group to most thermal")
        if not len(diffC)==len(sigA)==len(nSigF):
            raise Exception("Must have same group number for all GCs")
        
        if (len(diffC),len(diffC)) != np.shape(gtrans):
            raise Exception("Group transfer matrix has the wrong shape.")

        if isinstance(diffC[0],list):
            self.diffC = CSpline(tempvec, diffC)
            self.sigA = CSpline(tempvec, sigA)
            self.nSigF = CSpline(tempvec, nSigF)
            self.gtrans = CSpline(tempvec, gtrans)
            self.Tnom = 900.0 # default to 900K
        else:
            # only one value given, use scalar value
            self.diffC  = diffC
            self.sigA   = sigA
            self.nSigF  = nSigF
            self.gtrans = gtrans
            self.Tnom = 900.0 # default to 900K

    def set_Tnom(self, t):
        self.Tnom = t

    # getting cross sections @ nominal temperature
    def getSigA(self, group): 
        return self.sigA[group]

    def getDiffC(self, group):
        return self.diffC[group]

    def getNuSigF(self, group):
        return self.nSigF[group]

    def getGroupTransXS(self, gfrom, gto):
        return self.gtrans[gto][gfrom]

    # interpolation methods
    def interpSigA(self, group, T): 
        return self.sigA(T)[group]

    def interpgetDiffC(self, group, T):
        return self.diffC(T)[group]

    def interpNuSigF(self, group, T):
        return self.nSigF(T)[group]

    def interpGroupTransXS(self, gfrom, gto, T):
        return self.gtrans(T)[gto][gfrom]




def coord2enum(i, j, numi):
    """ Returns the enumeration from the mapping
    f: N x N -> N
    In other words, node coordinates on a uniform mesh
    are enumerated. """
    return (numi)*j + i

def enum2coord(num, numi):
    """ magically inverts the function coord2enum
    seems like magic because it inverts the function
    f: N x N -> N which isn't necessarily always possible."""
    j = num//numi
    i = num - j*numi
    return i,j
    

# next, the system A matrix governing diffusion and absorbtion gets formed.
# This could be done in a smart way, but instead I'll just loop through every
# point on the thing.


def makeAMatrix(r, h, nr, nh, numgroups, geommap):

    Alist = [np.zeros([nr * nh, nr * nh]) for i in range(numgroups)]
    Blist = [np.zeros(nr * nh) for i in range(numgroups)]
    dr = r / nr
    dz = h / nh

    # make the matrix with natural boundary conditions
    for A, B, gnum in zip(Alist, Blist, range(numgroups)):
        for i in range(nr):
            for j in range(nh):
                r = i * dr
                h = j * dz

                if r != 0:

                    # get on-diagonal term:
                    A[coord2enum(i, j, nr), coord2enum(i, j, nr)] = (2*geommap(r, h).getDiffC(gnum) / dz**2 +
                                                                    (r+dr)/(r) * geommap(r, h).getDiffC(gnum) / dr**2 +
                                                                    geommap(r, h).getDiffC(gnum) / dr**2 +
                                                                    geommap(r,h).getSigA(gnum) )

                    # above term:
                    if j!=nh-1:
                        A[coord2enum(i, j, nr), coord2enum(i, j+1, nr)] = -geommap(r, h).getDiffC(gnum) / dz**2
                                #(r+dr/2.0) * dz**2)
                    # below term:
                    if j!=0:
                        A[coord2enum(i, j, nr), coord2enum(i, j-1, nr)] = -geommap(r, h).getDiffC(gnum) / dz**2

                    # left term:
                    if i!=0:
                        A[coord2enum(i, j, nr), coord2enum(i-1, j, nr)] = -geommap(r, h).getDiffC(gnum) / dr**2

                    # right term:
                    if i!=nr-1:
                        A[coord2enum(i, j, nr), coord2enum(i+1, j, nr)] = -(r+dr)/r *geommap(r, h).getDiffC(gnum) / dr**2

                    B[coord2enum(i, j, nr)] = geommap(r, h).getNuSigF(gnum)

                else:
                    A[coord2enum(i, j, nr), coord2enum(i, j, nr)] = 1.0
                    A[coord2enum(i, j, nr), coord2enum(i+1, j, nr)] = -1.0
                    # # get on-diagonal term:
                    # A[coord2enum(i, j, nr), coord2enum(i, j, nr)] =(sigA *
                    #    dr**3 /8. * dz**2 + 2.0 * diffcoeff * dr**3/8.0 +diffcoeff*dr/2.0)

                    # # above term:
                    # if j!=nh-1:
                    #     A[coord2enum(i, j, nr), coord2enum(i, j+1, nr)] = 
                    # # # below term:
                    # if j!=0:
                    #     A[coord2enum(i, j, nr), coord2enum(i, j-1, nr)] =

                    # B[coord2enum(i, j, nr)] = nSigF * dz**2 / 8.0


    # enforce some boundary conditions
    for A, B, gnum in zip(Alist, Blist, range(numgroups)):
        for i in range(nr):
            for j in range(nh):
                if i==0 or j==0 or i==nr-1 or j==nh-1:

                    # or dirichlet zero at reactor edges
                    B[coord2enum(i, j, nr)] = 0.0

    return Alist, Blist

def getInscatterFromOtherGroups(fluxvec, gnum, dr, dh, numi, geommap):
    """ OK, so, loop through other flux groups, then calculate your
    inscatter term from other groups.
    1) if the thing diverges, be sure to not add to boundary stuff
    """
    inscat = np.zeros(len(fluxvec[0]))
    for gnump, flxp in enumerate(fluxvec):
        if gnum == gnump:
            continue
        for ii in range(len(flxp)):
            i,j = enum2coord(ii, numi) # mesh indices
            r = i * dr
            h = j * dh
            inscat[ii] += geommap(r, h).getGroupTransXS(gnump, gnum) * flxp[ii]
    return inscat

def fluxVector2Grid(fluxvec, r, h, nr, nh):
    """ takes a vector of flux (1D), and returns
    rows corresponding to r, h, flux """
    grid = np.zeros([nr*nh, 3])
    for i in range(nr):
        for j in range(nh):
            grid[coord2enum(i, j, nr), 2] = fluxvec[coord2enum(i, j, nr)]
            grid[coord2enum(i, j, nr), 0] = (r / nr) * i
            grid[coord2enum(i, j, nr), 1] = (h / nh) * j

    return grid

def inversePower(Alist, Blist, dr, dh, numi,geommap,adjoint=False):
    """ solves the generalized eigenvalue problem:
    A x = lam B x 
    for its smallest eigenvalue. We assume strictly
    lower triangular group-to-group scatter matrix."""
    eigtol = 1e-9
    print('getting LU')
    lulist = [la.lu_factor(A) for A in Alist]
    print('got LU')
    x0 = [np.ones(np.shape(Alist[0])[1]) for i in range(len(Alist))]
    x1 = [np.ones(np.shape(Alist[0])[1]) for i in range(len(Alist))]
    eig0 = 1.0
    eig1 = 1.0
    maxiter = 300
    thisiter = 0

    while 1:
        # convenient thing about power iteration. You just have to take a ratio
        # of linear mappings in the dual of the solution to get the right eigenvalue,
        # so you can just define k from a ratio of fast flux integrals. (or thermal, or..)
        eig0 = eig1 # save old

        # normalize by L-inf norm of fast flux
        x0 = [np.copy(flux)/np.max(x1[0]) for flux in x1]

        # now go through all the downscatter
        for A, B, gnum in zip(Alist, Blist, range(len(Alist))):

            # assume chi = 1 for fast, 0 all else:
            fissSource = sum([Bfss*flss for Bfss,flss in zip(Blist,x0)])
            fissSource = fissSource if gnum==0 else 0.0

            inscat = getInscatterFromOtherGroups(x1,gnum, dr , dh, numi, geommap)

            if debug:
                print("group {}".format(gnum))
                print("fission source: {}".format(np.sum(fissSource)))
                print("inscatter term: {}".format(inscat))
                print("total inscatter: {}".format(sum(inscat)))
            x1[gnum] = la.lu_solve( lulist[gnum], fissSource + inscat)

        if debug:
            if np.max(x1[0])==0.0 or np.max(x0[0])==0.0:
                print("Broke because debugging, to inspect flux at"
                      "spurious eigenvalue.")
                break

        eig1 = np.max(x1[0]) / np.max(x0[0]) 
        if debug:
            print('k = {}\n'.format(1./eig1))
        thisiter+=1
        if np.abs(eig1-eig0) < 1e-9 or thisiter > maxiter:
            break

    print("Took {} power iterations.".format(thisiter))

    return eig1, x1

def secant(x0, x1, y0, y1):
    return (x0*y1-x1*y0)/(y1-y0)
