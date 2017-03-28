#!/usr/bin/env python
import sys
import numpy as np
#import dd_20140409 as dd, ww_20140403 as ww
#import matplotlib.pyplot as plt; plt.ion()
#import matplotlib as mpl
#import argparse
#from pylab import specgram
#from scipy.interpolate import interp1d
#import warnings
#import random
#import getpass
from IPython import embed
#from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import cdist

sin, cos, norm, dot, cross = np.sin, np.cos, np.linalg.norm, np.dot, np.cross

class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def getLOSfromSPRD(year):
    lines = [l for l in open('/afs/ipp/u/sprd/loscoord/LOS_COORD_%4i'%year).readlines() if 'MSS' in l]
    s0 = np.zeros((60,3))
    s1 = np.zeros((60,3))
    for l in lines:
        R0, phi0, z0, R1, phi1, z1 = np.array(l.split()[1:],float)
        row, column = l.split()[0][5:].split('L')
        column = int(column.replace("'", ""))
        row = int(row)
        i = (row-1)*6+column-1
        s0[i] = R0*cos(phi0/180.*np.pi), R0*sin(phi0/180.*np.pi), z0
        s1[i] = R1*cos(phi1/180.*np.pi), R1*sin(phi1/180.*np.pi), z1
    
    ds = s1-s0; ds = (ds.T/norm(ds, axis=1)).T
    return s0, ds

def getNBIgeo():
    ''' NBI geometry from bgeiger '''
    src = np.array([6.1870004,   -6.8164635,  -0.69000000])
    v = np.array([-1401832.2, 1904877., 229791.33]) 
    v /= norm(v)
    return src+7*v, v # don't need full beam length as we only care for MSE LOS, so move closer

def getApproxIntersections(s0, s1, v0, v1, beamSubdivisionLength):
    ''' determine MSE LsOS interections with NBI '''
    sn = norm(s1 - s0, axis=1)/beamSubdivisionLength
    vn = norm(v1 - v0)/beamSubdivisionLength
    vmat = np.zeros((int(vn), 3))
    for j in xrange(3):
        vmat[:,j] = np.linspace(v0[j], v1[j], vn)
    intersections = np.zeros_like(s0)
    for i in xrange(s0.shape[0]):
        if np.isnan(sn[i]): continue
        smat = np.zeros((int(sn[i]), 3))
        for j in xrange(3):
            smat[:,j] = np.linspace(s0[i,j], s1[i,j], sn[i])
        dmat = cdist(smat, vmat, 'euclidean') # calculates n*bn matrix of all distances
        s_ind, v_ind = np.unravel_index(dmat.argmin(), dmat.shape)
        #cv = vmat[v_ind]
        cs = smat[s_ind]
        intersections[i] = cs #(cv+cs)/2.
    return intersections

def getMSExperp(s0, s):
    ''' determine MSE switchboard x-vector in torus space from MSE LsOS '''
    # nah, we'll just use the 2015 value from now on; code below for reference
    return np.array([-0.89054714,  0.41783235,  0.17983856])

    s1 = s0+2*s
    p0 = np.average(s0, axis=0)
    p1 = np.average(s1, axis=0)
    dp = np.average(s, axis=0)
    dp = dp/norm(dp)

    mps = []
    for i in xrange(s0.shape[0]):
        d = dot((p1-s0[i]), dp) / dot(s[i], dp)
        mp = d*s[i] + s0[i]
        mps.append(mp)
    mps = np.array(mps)

    vecs = []
    for offset in range(6):
        for i, m in enumerate(mps[offset::6]):
            if i+1 > 9: break
            n = mps[offset::6][i+1]
            vecs.append((n-m)/norm(n-m))
    vecs = np.array(vecs)
    vec = np.average(vecs, axis=0)
    vec /= norm(vec)
    return vec

def makeRzAs(year=2014, plot=False, beamSubdivisionLength=1e-3, channels=range(60), withMSEvec=False,
    outputplot=False):
    '''generates MSE R, z, A coords for pi and sigma from a FARO measurement'''

    #year = 2015 if year == 2017 else year

    s0, s = getLOSfromSPRD(int(year))

    s1 = s0 + 2.2*s
    s0 = s0 + 1.0*s

    v0, v = getNBIgeo()
    v1 = v0 + v

    isecs = getApproxIntersections(s0, s1, v0, v1, beamSubdivisionLength)


    toplot = {}
    toplot['LOS'] = []

    if plot or outputplot:
        import kk_abock as kk
        try:
            if plot: from mayavi import mlab # needs intel/12.1!!!!
        except Exception, e:
            print 'mayavi failed to load, are you using intel/12.1?'
            raise e
        for i in xrange(60):
            ts0 = s0 - 1.5*s
            if plot: mlab.plot3d([ts0[i,0], s1[i,0]], [ts0[i,1], s1[i,1]], [ts0[i,2], s1[i,2]], tube_radius=None, color=(1,0,0))
            toplot['LOS'].append([[ts0[i,0], s1[i,0]], [ts0[i,1], s1[i,1]], [ts0[i,2], s1[i,2]]])

        if plot: mlab.points3d(isecs[:,0], isecs[:,1], isecs[:,2])
        if plot: mlab.plot3d([v0[0], v1[0]], [v0[1], v1[1]], [v0[2], v1[2]], color=(0,1,0))
        if plot: mlab.view(0, 0, distance=10, focalpoint=[0,0,0])

        toplot['isecs'] = [isecs[:,0], isecs[:,1], isecs[:,2]]
        toplot['beam'] = [[v0[0], v1[0]], [v0[1], v1[1]], [v0[2], v1[2]]]

        for R in [1.65]: #, 1.8, 1.95]: # axis
            x = R * cos(np.linspace(0,2*np.pi, 128))
            y = R * sin(np.linspace(0,2*np.pi, 128))
            if plot: mlab.plot3d(x,y,[0]*len(x), tube_radius=None)
            toplot['axis'] = [x,y,[0]*len(x)]

        eq = kk.kk()
        eq.Open(31113, 'AUGD', 'EQI')
        sepR, sepz = [], []
        for a in xrange(360):
            res = eq.rhopol_to_Rz(3.0, 0.99, a, True)
            sepR.append(res['R'])
            sepz.append(res['z'])
        sepR = np.array(sepR)
        sepz = np.array(sepz)
        # separatrix:
        phi = -0.37
        if plot: mlab.plot3d(np.cos(phi)*sepR, np.sin(phi)*sepR, sepz, tube_radius=None)
        toplot['sep'] = [np.cos(phi)*sepR, np.sin(phi)*sepR, sepz]

    xperp = getMSExperp(s0, s)

    isecR = (isecs[:,0]**2 + isecs[:,1]**2)**0.5
    isecZ = isecs[:,2]
    isecPhi = np.arctan2(isecs[:,1], isecs[:,0]) 

    eR = np.zeros_like(isecs)
    eR[:,0], eR[:,1] = cos(isecPhi), sin(isecPhi)
    eT = np.zeros_like(isecs)
    eT[:,0], eT[:,1] =  -sin(isecPhi), cos(isecPhi)
    eZ = np.zeros_like(isecs)
    eZ[:,2] = 1.

    sR = np.sum(s * eR, axis=1)
    st = np.sum(s * eT, axis=1)
    sz = np.sum(s * eZ, axis=1)

    xpR = np.sum(xperp * eR, axis=1)
    xpt = np.sum(xperp * eT, axis=1)
    xpz = np.sum(xperp * eZ, axis=1)

    vR = np.sum(v * eR, axis=1)
    vt = np.sum(v * eT, axis=1)
    vz = np.sum(v * eZ, axis=1)

    A1  = -st*vt*xpR - sz*vz*xpR + sR*vt*xpt + sR*vz*xpz
    A2  = -(sR*vR + sz*vz)*xpt + st*(vR*xpR+vz*xpz)
    A3  = sz*(vR*xpR + vt*xpt) - (sR*vR + st*vt)*xpz
    A4  = sz*xpt - st*xpz
    A5  = st*xpR - sR*xpt
    A6  = sR*(sz*vt - st*vz)*xpR - (st*vt + sz*vz)*(-sz*xpt + st*xpz) + (sR**2)*(vz*xpt - vt*xpz)
    A7  = -sz**2*vz*xpR - st*sz*vR*xpt + sR**2*vR*xpz + st**2*(-vz*xpR + vR*xpz) + sR*(-sz*vR*xpR + st*vz*xpt + sz*vz*xpz)
    A8  = st**2*vt*xpR - sR**2*vR*xpt + sz**2*(vt*xpR - vR*xpt) + st*sz*vR*xpz + sR*(st*vR*xpR - st*vt*xpt - sz*vt*xpz)
    A9  = (st**2*xpR + sz**2*xpR - sR*st*xpt - sR*sz*xpz)
    A10 = (-sR*sz*xpR + sR**2*xpz + st*(-sz*xpt + st*xpz))

    Api    = -np.array((A1, A2, A3, A4, A5, A6, A7, A8, A9, A10)).T
    Asigma = -np.array((A6, A7, A8, A9, A10, -A1, -A2, -A3, -A4, -A5)).T

    output = Bunch(
                R = np.zeros(60),
                z = np.zeros(60),
                Asigma = np.zeros((60,10)),
                Api = np.zeros((60,10)),
              )

    # map to our nomenclature, i.e. top right on switchbard = first channel, bottom left last channel
    # remember the directions are flipped in the torus, i.e. top line becomes bottom, rightmost column
    # becomes outermost LOS, etc.
    newOrder = []
    order = isecZ.argsort()
    for i in xrange(5,-1,-1):
        o = order[10*i:10*(i+1)]
        order2 = isecR[o].argsort()
        #print o[order2]
        newOrder += list(o[order2])
    newOrder.reverse()

    output.R = isecR[newOrder]
    output.z = isecZ[newOrder]
    output.Asigma = Asigma[newOrder]
    output.Api = Api[newOrder]

    if withMSEvec:
        Zmse = s
        Ymse = cross(Zmse, xperp)
        Xmse = cross(Ymse, Zmse)
        return output, (Xmse, Ymse, Zmse, xperp)
    elif outputplot:
        return output, toplot
    else:
        return output


if __name__ == '__main__':
    chans1 = [10, 11, 12, 13, 14, 15, 16, 17, 28, 39]
    chans2 = [20, 21, 22, 23, 24, 25, 26, 27, 38, 49]
    #chans1 = [1, 7, 13, 19, 25, 31, 37, 43, 50, 57]
    #chans2 = [2, 8, 14, 20, 26, 32, 38, 44, 51, 58]
    print 'new...'
    asd = makeRzAs(2015, plot=True)

    A = (asd.Asigma[chans1] + asd.Asigma[chans2]) / 2.
    def gammas(Bt, Bz, Ez, A):
        g_m = np.arctan2(A[:,1]*Bt + A[:,2]*Bz + A[:,4]*Ez, A[:,6]*Bt + A[:,7]*Bz + A[:,9]*Ez)
        return g_m/np.pi*180

    #embed()
    #sys.exit()

    print 'old...'
    from RzAmakerOld import makeRzAs as oldmrsa
    asd2 = oldmrsa(mseFaroFile='mse2015.txt', plot=False)
    #sys.exit()
    Asigma = (asd.Asigma[chans1] + asd.Asigma[chans2]) / 2.
    Asigma2 = (asd2.Asigma[chans1] + asd2.Asigma[chans2]) / 2.
    Api = (asd.Api[chans1] + asd.Api[chans2]) / 2.
    Api2 = (asd2.Api[chans1] + asd2.Api[chans2]) / 2.
    #Alist = [('newsig', Asigma), ('oldsig', Asigma2)] #, ('newpi ', Api), ('oldpi ', Api2)]
    #for i in xrange(10):
    #    for name, As in Alist:
    #        print '%s %2i '%(name, i+1),
    #        for j in xrange(10):
    #            print '% 5.4f '%As[i, j],
    #        print 
    #        #break
    newsig = gammas(1,0,0, Asigma)
    oldsig = gammas(1,0,0, Asigma2)
    newpi = gammas(1,0,0, Api)
    oldpi = gammas(1,0,0, Api2)
    # note 4degree variation across minor radius between old and new angles

    embed()



#             A1*Br + A2*Bt + A3*Bz + A4*Er/vb + A5*Ez/vb
# tan g_m =  ---------------------------------------------
#             A6*Br + A7*Bt + A8*Bz + A9*Er/vb + A10*Ez/vb

#             A[0]*Br + A[1]*Bt + A[2]*Bz + A[3]*Er/vb + A[4]*Ez/vb
# tan g_m =  ------------------------------------------------------
#             A[5]*Br + A[6]*Bt + A[7]*Bz + A[8]*Er/vb + A[9]*Ez/vb



