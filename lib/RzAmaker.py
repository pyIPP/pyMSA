#!/usr/bin/env python
#import os, sys
import numpy as np
#import dd_20140409 as dd, ww_20140403 as ww
import kk_abock as kk
#import matplotlib.pyplot as plt
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
from copy import copy

sin, cos = np.sin, np.cos

class Bunch(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def makeRzAs(mseFaroFile='mse2014.txt', plot=False, beamSubdivisionLength=1e-3):
    '''generates MSE R, z, A coords for pi and sigma from a FARO measurement'''
    # shorter = faster, too short = wrong result:
    # check with plot option above
    nbiBeamLength = 1.2
    mseBeamLength = 2. 

    btres = beamSubdivisionLength # beam point interval in m

    # from David Rittich/Felician Mink/Christian Hopf
    Q3 = np.array([[9.248975, 20.3955, -0.600], # R, phi (0 south), z
                   [2.37409,  37.7203,  0]])

    Q3[:,1] = Q3[:,1]/180.*np.pi # rotate into common coordinate system

    Q3c = np.zeros_like(Q3) # create cartesian coordinate array
    # fill it...
    Q3c[:,0] = Q3[:, 0]*cos(Q3[:,1])
    Q3c[:,1] = Q3[:, 0]*sin(Q3[:,1])
    Q3c[:,2] = Q3[:,2]
    vdx = Q3c[1,0] - Q3c[0,0]
    vdy = Q3c[1,1] - Q3c[0,1]
    vdz = Q3c[1,2] - Q3c[0,2]
    vd = (vdx**2+vdy**2+vdz**2)**0.5 / nbiBeamLength
    # use original inner point and another point inside plasma
    Q3c[:] = Q3c[1], Q3c[1] + np.array([vdx/vd, vdy/vd, vdz/vd])

    if plot:
        try:
            from mayavi import mlab # needs intel/12.1!!!!
        except Exception, e:
            print 'mayavi failed to load, are you using intel/12.1?'
            raise e
        mlab.plot3d(Q3c[:,0], Q3c[:,1], Q3c[:,2], color=(0,1,0))
        eq = kk.kk()
        eq.Open(31113, 'AUGD', 'EQI')
        sepR, sepz = [], []
        for a in xrange(360):
            res = eq.rhopol_to_Rz(3.0, 0.99, a, True)
            #print res
            sepR.append(res['R'])# [0])
            sepz.append(res['z'])# [0])
        sepR = np.array(sepR)
        sepz = np.array(sepz)
        mlab.plot3d(np.cos(Q3[1,1])*sepR, np.sin(Q3[1,1])*sepR, sepz, tube_radius=None)


    # load MSE coords from Faro
    msecoords = np.loadtxt(mseFaroFile, delimiter=';', usecols=(1,2,3))/1e3

    # buffer for better suitable coords
    nmc = np.zeros_like(msecoords) # new mse coords
    # fill it...
    for i in xrange(60):
        x = msecoords[[i,i+60], 0]
        y = msecoords[[i,i+60], 1]
        z = msecoords[[i,i+60], 2]
        # rotate into common coordinate system...
        R = (x**2 + y**2)**0.5
        phi = np.arctan2(y, x) + 6*np.pi/16.
        x = R*cos(phi)
        y = R*sin(phi)
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        dz = z[1] - z[0]
        d = (dx**2+dy**2+dz**2)**0.5 / mseBeamLength

        sp = 1
        nmc[i,:] = x[sp]-dx/d*2., y[sp]-dy/d*2., z[sp]-dz/d*2.
        nmc[i+60,:] = x[sp]+dx/d, y[sp]+dy/d, z[sp]+dz/d
        if plot:
            mlab.plot3d(nmc[[i,i+60],0],nmc[[i,i+60],1],nmc[[i,i+60],2], line_width=1.,
                tube_radius=None, color=(1,0,0))
        

    # plot magnetic axis for reference
    if plot:
        x = 1.65 * cos(np.linspace(0,2*np.pi, 128))
        y = 1.65 * sin(np.linspace(0,2*np.pi, 128))
        mlab.plot3d(x,y,[0]*len(x))
        mlab.view(0, 0, distance=10, focalpoint=[0,0,0])

    # calculate intersections and local A's:

    # delta: between v and x,y plane, i.e. independent of MSE
    vxy = (vdx**2 + vdy**2)**0.5
    vd = (vdx**2+vdy**2+vdz**2)**0.5
    delta = np.arccos(vxy/vd)

    # 1st: generate beam points:
    bn = nbiBeamLength / btres
    bx = np.linspace(Q3c[0,0], Q3c[1,0], bn)
    by = np.linspace(Q3c[0,1], Q3c[1,1], bn)
    bz = np.linspace(Q3c[0,2], Q3c[1,2], bn)
    # make a matrix from points
    bmat = np.zeros((bn, 3))
    bmat[:,0], bmat[:,1], bmat[:,2] = bx, by, bz

    # preliminary output buffers
    R = np.zeros(60)
    Z = np.zeros(60)
    Asigma = np.zeros((60,10))
    Api = np.zeros((60,10))

    for i in xrange(60):
        n = mseBeamLength / btres
        x = np.linspace(nmc[i,0],nmc[i+60,0], n)
        y = np.linspace(nmc[i,1],nmc[i+60,1], n)
        z = np.linspace(nmc[i,2],nmc[i+60,2], n)
        # make a matrix from points
        pmat = np.zeros((n, 3))
        pmat[:,0], pmat[:,1], pmat[:,2] = x, y, z

        dmat = cdist(bmat, pmat, 'euclidean') # calculates n*bn matrix of all distances
        # get indices of points closest to each other
        beamindex, mseindex = np.unravel_index(dmat.argmin(), dmat.shape)

        r, z = (x[mseindex]**2 + y[mseindex]**2)**0.5, z[mseindex] # R, z

        R[i] = r 
        Z[i] = z

        # now: A1-10 for pi and sigma
        # alpha: between (vr, vt, 0) and (sr, st, 0) (i.e. between v and s in x,y plane)
        s = nmc[i+60,:]-nmc[i,:]
        sxy = (s[0]**2+s[1]**2)**0.5
        vz0 = np.array([vdx, vdy])
        alpha = np.arccos(np.dot(vz0,s[:2])/vxy/sxy)
        # theta: between s and x,y plane
        sz0 = copy(s)
        sz0[2] = 0 
        theta = np.arccos(np.dot(s, sz0)/np.linalg.norm(sz0)/np.linalg.norm(s))
        # Omega: between (sr, st, 0) and e_T
        phi = np.arctan2(y[mseindex], x[mseindex])
        eT = np.array([-sin(phi),cos(phi),0])
        Omega = np.arccos(np.dot(sz0, eT)/np.linalg.norm(sz0))

        # sigma
        A1  =  sin(delta) * sin(Omega)
        A2  = -cos(Omega) * sin(delta)
        A3  =  cos(delta) * cos(alpha-Omega)
        A4  =  cos(Omega)
        A5  =  0
        A6  =  (-cos(alpha) * cos(delta)) * cos(theta) - (cos(Omega) * sin(delta)) * sin(theta)
        A7  =  (-cos(delta) * cos(theta)) * sin(alpha) - (sin(delta) * sin(theta)) * sin(Omega)
        A8  =  (-cos(delta) * sin(theta)) * sin(alpha-Omega)
        A9  =   (sin(theta) * sin(Omega))
        A10 = cos(theta)
        Asigma[i, :] = A1, A2, A3, A4, A5, A6, A7, A8, A9, A10

        # pi
        A1  =   (cos(alpha) * cos(delta)) * cos(theta) + (cos(Omega) * sin(delta)) * sin(theta)
        A2  =   (cos(delta) * cos(theta)) * sin(alpha) + (sin(delta) * sin(theta)) * sin(Omega)
        A3  =   (cos(delta) * sin(theta)) * sin(alpha-Omega)
        A4  = -((sin(theta) * sin(Omega)))
        A5  = -(cos(theta))
        A6  =   sin(delta) * sin(Omega)
        A7  =  -cos(Omega) * sin(delta)
        A8  =   cos(delta) * cos(alpha-Omega)
        A9  =   cos(Omega)
        A10 = 0
        Api[i, :] = A1, A2, A3, A4, A5, A6, A7, A8, A9, A10

    # prepare final output buffers
    output = Bunch(
                    R = np.zeros(60),
                    z = np.zeros(60),
                    Asigma = np.zeros((60,10)),
                    Api = np.zeros((60,10))
                  )
    # map to our nomenclature, i.e. top right first channel, bottom left last channel
    zOrder = Z.argsort() # get z order, then iterate from top row to bottom row:
    for i in xrange(5,-1,-1):
        rs = R[zOrder[i*10:(i+1)*10]]
        zs = Z[zOrder[i*10:(i+1)*10]]
        asi = Asigma[zOrder[i*10:(i+1)*10]]
        api = Api[zOrder[i*10:(i+1)*10]]
        rOrder = rs.argsort() # get R order in current row, and apply it:
        rs = rs[rOrder][::-1]
        zs = zs[rOrder][::-1]
        asi = asi[rOrder][::-1]
        api = api[rOrder][::-1]
        # put into output buffer
        output.R[(-i+5)*10:(-i+6)*10] = rs
        output.z[(-i+5)*10:(-i+6)*10] = zs
        output.Asigma[(-i+5)*10:(-i+6)*10] = asi
        output.Api[(-i+5)*10:(-i+6)*10] = api

    return output


if __name__ == '__main__':
    asd = makeRzAs(mseFaroFile='mse2015.txt', plot=True)
    embed()








