#!/usr/bin/python
import os, os.path, sys
import numpy as np
from copy import deepcopy

np.seterr(all='raise')
def lagrange(points):
    """
    return lagrange interpolating polynomial
    usually we use len(points) = 3, i.e. quadratic interpolation
    """
    def P(x):
        total = 0
        n = len(points)
        for i in xrange(n):
            xi, yi = points[i]

            def g(i, n):
                tot_mul = 1
                for j in xrange(n):
                    if i == j:
                        continue
                    xj, yj = points[j]
                    tot_mul *= (x - xj) / float(xi - xj)

                return tot_mul

            total += yi * g(i, n)
        return total
    return P

def distance(x12):
    """
    Calculate atom pair distance.
    x12: [[x1,y1,z1],[x2,y2,z2]]
    """

    dx = x12[0][0]-x12[1][0]
    dy = x12[0][1]-x12[1][1]
    dz = x12[0][2]-x12[1][2]
    result = sqrt(dx*dx+dy*dy+dz*dz)

    return result

def translate(vec, dvec):
    return [vec[i]+dvec[i] for i in range(3)]

def rotate(vec, axis, theta):
    """
    According to "Euler-Rodrigues formula"
    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    There was a negative sign '-' in front of the right term in "a = -np.cos(theta*0.5)"
    but not in Wikipedia formula.
    It seems, without the '-' sign results in "right hand role".
    """
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta*0.5)
    b,c,d = axis*np.sin(theta*0.5)
    rotmat = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                       [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                       [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
    return np.dot( rotmat, vec )

def get_normal(vec1, vec2):
    return np.cross(vec1, vec2)

def get_unit(vec):
    vec = np.array(vec)
    return vec/np.sqrt((vec*vec).sum())

def get_bisect(vec1, vec2):
    v1 = get_unit(vec1)
    v2 = get_unit(vec2)
    return [(v1[i]+v2[i])*0.5 for i in range(3)]

def get_normal_unit(vec1, vec2):
    vec = np.cross(vec1, vec2)
    return get_unit(vec)

def get_bisect_unit(vec1, vec2):
    vec = np.array( [(vec1[i]+vec2[i])*0.5 for i in range(3)] )
    vec /= np.sqrt((vec*vec).sum())
    return vec

def length(vec):
    return np.sqrt(np.dot(vec, vec))

def angle(vec1, vec2):
    """
    unit: radius
    """
    dotprod = np.dot(vec1, vec2)
    return np.arccos(dotprod / (length(vec1) * length(vec2)))

## for water probe
Mol2BisectDataNdx = { '111': {0:[ 4,20,19, 3], 1:[25,18,19,20], 2:[ 2, 3,19,18]},
                     '100': {0:[ 4,12,13, 5], 1:[24,14,13,12], 2:[ 6, 5,13,14]},
                     '010': {0:[ 0,16,23, 7], 1:[25,22,23,16], 2:[ 6, 7,23,22]},
                     '001': {0:[ 0, 8, 9, 1], 1:[24,10, 9, 8], 2:[ 2, 1, 9,10]},
                     '110': {0:[ 4, 5,21,20], 1:[25,20,21,22], 2:[ 6,22,21, 5]},
                     '101': {0:[ 4, 3,11,12], 1:[24,12,11,10], 2:[ 2,10,11, 3]},
                     '011': {0:[ 0, 1,17,16], 1:[25,16,17,18], 2:[ 2,18,17, 1]},
                     '000': {0:[ 0, 7,15, 8], 1:[24, 8,15,14], 2:[ 6,14,15, 7]}}

AxisInQuadrant = {  '111': [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]],
                  '100': [[ 1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]],
                  '010': [[-1, 0, 0],[ 0, 1, 0],[ 0, 0,-1]],
                  '001': [[-1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]],
                  '110': [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0,-1]],
                  '101': [[ 1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]],
                  '011': [[-1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]],
                  '000': [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]]}


def lagrange_interp(points,x):
    # points = [(ang1,ener1),(ang2,ener2),(ang3,ener3)]
    # x: the query angle
    # return the energy at point x
    P = lagrange(points)
    return P(x)


def linear1(dxx0, dx1x0, y0, y1):
    return y0 + (y1-y0)*dxx0/dx1x0

def bilinear_gen(y0,y1,y2,y3,angx1,angx2,angy,label=1): #for general assignment of subsection
    """
    Bilinear (Lifeng), label = 0:

      y3<----y2
      ^
      ya  | ---
      |  |angy
      y0---->y1---
      |<->|
      angx

    """

    """
    Bilinear (general, Lan), label = 1:

          angx2
          |<->|
          y2<---y3
          ^
          ya  | ---
          |  |angy
          y0------>y1---
          |<--->|
          angx1

    """

    y01 = linear1( angx1, 1, y0, y1 )
    if label == 1:
        y23 = linear1( angx2, 1, y2, y3 )
    else:
        y23 = linear1( angx2, 1, y3, y2 )

    ya  = linear1( angy, 1, y01, y23 )

    return ya


def quadrant_ndx(vec):
    ndx = ['1','1','1']
    for i in range(3):
        if vec[i]<0: ndx[i] = '0'
        return ndx[0]+ndx[1]+ndx[2]


def weights_for_normal_general(norm_vec, config_vecs, cut=0.0000001):
    """
    Return the weights of the proper two configurations at one bisector direction.
Linear interpolation is used.
"""

    nvec = len(config_vecs)
    vec1 = config_vecs[0]
    vec2 = config_vecs[1]
    ## treat vec1 as the x-axis, the normal of (vec1,vec2) as the z-axis.
    vx = vec1
    vz = get_normal_unit(vec1,vec2)
    vy = np.cross(vz,vx)

    ## project the "norm_vec" to the plane of vec1 and vec2, then use symmetry operation.
    px = np.dot(norm_vec,vx)
    py = np.dot(norm_vec,vy)
    if abs(px)<cut:
        if py>0: w2=1.0
        else: w2=-1.0
        w1 = 0.0
        ndx2 = (np.pi/2)/(np.pi/nvec)
        ndx2 = int(ndx2)
        ndx1 = ndx2 - 1
    else:
        ang = np.arctan(py/px)
        if px<0: ang += np.pi
        elif py<0: ang += 2*np.pi
        period = np.pi
        nrot = ang//period
        rrot = ang%period

    ndx1 = rrot//(np.pi/nvec)
    ang2 = rrot%(np.pi/nvec)
    w2 = ang2/(np.pi/nvec)
    w1 = 1-w2
    if ndx1 == (nvec-1):
        ndx2 = 0
    else:
        ndx2 = ndx1 + 1

    return w1, w2, int(ndx1), int(ndx2)


def get_neighors_for_normal(norm_vec, config_vecs, cut=0.0000001):
    """
    return the angle for query normal vector and the angles of its neighbor (n=3)
"""

    nvec = len(config_vecs)
    vec1 = config_vecs[0]
    vec2 = config_vecs[1]
    ## treat vec1 as the x-axis, the normal of (vec1,vec2) as the z-axis.
    vx = vec1
    vz = get_normal_unit(vec1,vec2)
    vy = np.cross(vz,vx)

    period = np.pi
    da = period/nvec

    ## project the "norm_vec" to the plane of vec1 and vec2, then use symmetry operation.
    px = np.dot(norm_vec,vx)
    py = np.dot(norm_vec,vy)
    if abs(px)<cut:
        ang0 = np.pi/2
        ndx2 = ang0/da
        ndx2 = int(ndx2)
        ndx1 = ndx2 - 1
        if ndx2 == (nvec-1):
            ndx3 = 0 + nvec
            if ndx1 == 0: sys.exit() # lagrange interpolation needs three different points
        else:
            ndx3 = ndx2 + 1
    else:
        ang = np.arctan(py/px)
        if px<0: ang += np.pi
        elif py<0: ang += 2*np.pi
        nrot = ang//period
        rrot = ang%period
        ang0 = rrot

        ndx1 = rrot//da
        if ndx1 == (nvec-1):
            ndx2 = 0 + nvec
            ndx3 = ndx1 - 1
            if ndx3 == 0: sys.exit()
        else:
            ndx2 = ndx1 + 1
            tmp1 = (ndx2+1)*da - ang0
            tmp2 = ang0 - (ndx1-1)*da
            if abs(tmp1) < abs(tmp2):
                if ndx2 == (nvec-1):
                    ndx3 = 0 + nvec
                else:
                    ndx3 = ndx2 + 1
            else:
                ndx3 = ndx1 - 1

    return ang0, int(ndx1), int(ndx2), int(ndx3)


def weights_in_subsection(bisvec, cutoff=0.9999):
    """
    As the bisector directions are fixed, it is not necessary to
    specify indices of layer, grid, or others.

    "bisvec" is the bisector of mole2, and should be a unit vector.
    1. Find the sub-section in quadrant 'quad_ndx';
    2. Calculate the angle distances to the four conner vectors in the sub-section;
    """
    quad_ndx = quadrant_ndx(bisvec)
    angcos = [np.dot(bisvec, AxisInQuadrant[quad_ndx][0]),
              np.dot(bisvec, AxisInQuadrant[quad_ndx][1]),
              np.dot(bisvec, AxisInQuadrant[quad_ndx][2])]
    maxcos = angcos[0]
    subndx = 0
    ## The angle with max cos is the smallest:
    for i in [1,2]:
        if angcos[i]>maxcos:
            maxcos = angcos[i]
            subndx = i

    ## Bilinear interpolation is used
    ##     for the interpolation in a sub-section of a quadrant:
    ## Three SUBSECTIONS: below is a quadrant with three sub-sections:
    #       z
    #      / \
        #     /   \
        #    c     b
    #   /  \o/  \
        #  /    |    \
        # x-----a-----y
    #
    # subsection 0: x-a-o-c (0-1-4-3)
    #            1: y-b-o-a
    #            2: z-c-o-b
    ## indices of axis: 0, 1, 2, 3, 4 (3,4 for cyclic call)
    ax = {  '111':{0:[0,1,2], 1:[1,2,0], 2:[2,0,1]},
          '100':{0:[0,1,2], 1:[1,2,0], 2:[2,0,1]},
          '010':{0:[0,1,2], 1:[1,2,0], 2:[2,0,1]},
          '001':{0:[0,1,2], 1:[1,2,0], 2:[2,0,1]},
          '110':{0:[0,2,1], 1:[1,0,2], 2:[2,1,0]},
          '101':{0:[0,2,1], 1:[1,0,2], 2:[2,1,0]},
          '011':{0:[0,2,1], 1:[1,0,2], 2:[2,1,0]},
          '000':{0:[0,2,1], 1:[1,0,2], 2:[2,1,0]}}

    ## angx is calculated by arctan(dy/dx):
    tempx = bisvec[ax[quad_ndx][subndx][0]]
    tempy = bisvec[ax[quad_ndx][subndx][1]]
    angx = tempy / tempx
    angx = np.arctan( abs(angx) )

    ## angx is calculated by arctan(dy/dx):
    tempz = bisvec[ ax[quad_ndx][subndx][2] ] # i+1 is the index of "z" axis
    angy = np.arcsin( abs(tempz) )

    n0 = Mol2BisectDataNdx[quad_ndx][subndx][0]
    n1 = Mol2BisectDataNdx[quad_ndx][subndx][1]
    n2 = Mol2BisectDataNdx[quad_ndx][subndx][2]
    n3 = Mol2BisectDataNdx[quad_ndx][subndx][3]

    return quad_ndx, [n0,n1,n2,n3], angx, angy

