#!/usr/bin/python

import numpy as np
from copy import deepcopy

np.seterr(all='raise')

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


def convert_coords( coord ):
    """
    Reorient the coords as to set mole1 to origon and in xy face.
    """
    ## translate:
    dx1 = np.array(coord[0]) * -1
    tcoord = []
    for ii in range(6): tcoord.append(translate(coord[ii], dx1))

    ## rot to z-axis
    norm_m1 = get_normal(tcoord[1], tcoord[2])
    ang1 = angle(norm_w1, [0,0,1])
    norm_r1 = get_normal(norm_m1, [0,0,1])
    for ii in range(6): tcoord[ii] = rotate(tcoord[ii], norm_r1, ang1)

    ## rot to x-axis
    bis_m1 = get_bisect(tcoord[1], tcoord[2])
    ang2 = angle(bis_m1, [1,0,0])
    for ii in range(6): tcoord[ii] = rotate(tcoord[ii], [0,0,1], ang2)

    ## calt norm and bisector of mole2
    dx2 = np.array(tcoord[3])
    mol2h1 = translate(tcoord[4], -dx2)
    mol2h2 = translate(tcoord[5], -dx2)
    norm_m2 = get_normal(mol2h1, mol2h2)
    bis_m2 = get_bisect( mol2h1, mol2h2 )

    ## calt r,theta,phi of the second molecule
    r = length( dx2 )
    theta = angle( dx2, [dx2[0],dx2[1],0] )
    phi = angle( [1,0,0], [dx2[0],dx2[1],0] )


## '1' means >0, '0' means <0. So '111' means the 1st quadrant in all the 8 quadrants
##    in a space, while '000' means the 8th quadrant in the space.
## order of vectors in each quadrant (right hand, xy plane always first): 
#Mol2BisectDataNdx = {'111': [0,8,24,1,9,10,2],   # x,xy,y,z 
#              '100': [0,16,25,7,23,22,6], # x,x(-y),-y,-z 
#              '010': [4,12,24,5,13,14,6], # -x,(-x)y,y,-z  
#              '001': [4,20,25,3,19,18,2], # -x,(-x)(-y),-y,z 
#              '110': [24,8,0,14,15,7,6],  # y,yx,x,-z  
#              '101': [25,16,0,18,17,1,2], # -y,(-y)x,x,z  
#              '011': [24,12,4,10,11,3,2], # y,y(-x),-x,z 
#              '000': [25,20,4,22,21,5,6]} # -y,(-y)(-x),-x,-z 
## for water probe
Mol2BisectDataNdx = { '111': {0:[ 4,20,19, 3], 1:[25,18,19,20], 2:[ 2, 3,19,18]},
              '100': {0:[ 4,12,13, 5], 1:[24,14,13,12], 2:[ 6, 5,13,14]},
              '010': {0:[ 0,16,23, 7], 1:[25,22,23,16], 2:[ 6, 7,23,22]},
              '001': {0:[ 0, 8, 9, 1], 1:[24,10, 9, 8], 2:[ 2, 1, 9,10]},
              '110': {0:[ 4, 5,21,20], 1:[25,20,21,22], 2:[ 6,22,21, 5]},
              '101': {0:[ 4, 3,11,12], 1:[24,12,11,10], 2:[ 2,10,11, 3]},
              '011': {0:[ 0, 1,17,16], 1:[25,16,17,18], 2:[ 2,18,17, 1]},
              '000': {0:[ 0, 7,15, 8], 1:[24, 8,15,14], 2:[ 6,14,15, 7]}}
## for nh4,ch4 probe
#Mol2BisectDataNdx = { '111': {0:[ 0, 8,15, 7], 1:[24,14,15, 8], 2:[ 6, 7,15,14]},
#              '100': {0:[ 0,16,17, 1], 1:[16,25,18,17], 2:[ 2, 1,17,18]},
#              '010': {0:[ 4,12,11, 3], 1:[24,10,11,12], 2:[ 2, 3,11,10]},
#              '001': {0:[ 4,20,21, 5], 1:[25,22,21,20], 2:[ 6, 5,21,22]},
#              '110': {0:[ 0, 1, 9, 8], 1:[24, 8, 9,10], 2:[ 2,10, 9, 1]},
#              '101': {0:[ 0, 7,23,16], 1:[25,16,23,22], 2:[ 6,22,23, 7]},
#              '011': {0:[ 4, 5,13,12], 1:[24,12,13,14], 2:[ 6,14,13, 5]},
#              '000': {0:[ 4, 3,19,20], 1:[25,20,19,18], 2:[ 2,18,19, 3]}}

AxisInQuadrant = {  '111': [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]],
            '100': [[ 1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]],
            '010': [[-1, 0, 0],[ 0, 1, 0],[ 0, 0,-1]],
            '001': [[-1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]],
            '110': [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0,-1]],
            '101': [[ 1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]],
            '011': [[-1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]],
            '000': [[-1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]]}


## Linear interpolation:
#def linear(x0, y0, x1, y1, x):
#    return y0 + (y1-y0)*(x-x0)/(x1-x0)

def linear1(dxx0, dx1x0, y0, y1):
    return y0 + (y1-y0)*dxx0/dx1x0

def bilinear(y0,y1,y2,y3,angx,angy):
    """
    Bilinear:
    
      y3<----y2
             ^
         ya  | ---
             |  |angy
      y0---->y1---
      |<->|
       angx
    
    """

    ## np.pi/4.0:
    pi4 = 0.78539816339744817
    
    y01 = linear1( angx, pi4, y0, y1 )
    y32 = linear1( angx, pi4, y3, y2 )
    ya  = linear1( angy, pi4, y01, y32 )

    return ya


def quadrant_ndx(vec):
    ndx = ['1','1','1']
    for i in range(3):
        if vec[i]<0: ndx[i] = '0'
    return ndx[0]+ndx[1]+ndx[2]

def weights_for_2_configs(norm_vec, config1, config2, cut=0.0000001):
    """
    Return the weights of the two configurations at one bisector direction.
    Linear interpolation is used.
    NOTE: The two vectors of the two configurations are calculated here
        but in future, they should be retrieve from a pre-calculated 
        constant data list or dictionary of fixed directions.
    """
    vec1 = get_normal_unit(np.array(config1[1])-np.array(config1[0]),
              np.array(config1[2])-np.array(config1[0]))
    vec2 = get_normal_unit(np.array(config2[1])-np.array(config2[0]),
              np.array(config2[2])-np.array(config2[0]))


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
    else:
        ang = np.arctan(py/px)
        if px<0: ang += np.pi
        elif py<0: ang += 2*np.pi
        period = np.pi
        nrot = ang//period
        rrot = ang%period
        if (rrot-np.pi/2.0)>0.0:
            ang2 = np.pi-rrot                
        else: ang2 = rrot
        w2=ang2/(np.pi/2.0)
        w1=1-w2
    return w1, w2        


    #pi2 = np.pi*0.5
    #if abs(w2)<cut: 
    #    if w1<0.0: w1 = -pi2
    #    else: w1 = pi2
    #    w2 = 0.0
    #else: 
    #    t1 = np.arctan(abs(w1/w2))
    #    t2 = (pi2-t1)/pi2
    #    t1 = t1/pi2
    #    if w1<0: w1 = -t1
    #    else: w1 = t1
    #    if w2<0: w2 = -t2
    #    else: w2 = t2
    ## this sckeme is using "norm_vec" directly.
    # temp = abs(w1) + abs(w2)
    # w1 = w1 / temp
    # w2 = w2 / temp
    #return w1, w2

def weights_for_2_norms(nrmvec, vec1, vec2):
    """
    Return the weights of the two configurations at one bisector direction.
    Linear interpolation is used.
    NOTE: The two vectors of the two configurations are calculated here
        but in future, they should be retrieve from a pre-calculated 
        constant data list or dictionary of fixed directions.
    """
    w1 = abs(np.dot(nrmvec, vec1))
    w2 = abs(np.dot(nrmvec, vec2))
    temp = w1 + w2
    ## according to linear interpolation:
    # y = y2*w1+y1*w2 
    #    (w1 = (x-x1)/(x2-x1); w2 = (x2-x)/(x2-x1))
    # Thus, for the following data, y = y2*w1+y1*w2
    w1 = w1 / temp
    w2 = w2 / temp
    return w1, w2

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

    ## Bilinear:
    #
    #  y3<----y2
    #         ^
    #     ya  | ---
    #         |  |angy
    #  y0---->y1---
    #  |<->|
    #   angx
    #
    ### np.pi/4.0:
    #pi4 = 0.78539816339744817
    #
    #y01 = linear1( angx, pi4, y0, y1 )
    #y32 = linear1( angx, pi4, y3, y2 )
    #ya  = linear1( angy, pi4, y01, y32 )

    n0 = Mol2BisectDataNdx[quad_ndx][subndx][0]
    n1 = Mol2BisectDataNdx[quad_ndx][subndx][1]
    n2 = Mol2BisectDataNdx[quad_ndx][subndx][2]
    n3 = Mol2BisectDataNdx[quad_ndx][subndx][3]

    return quad_ndx, [n0,n1,n2,n3], angx, angy

if __name__=='__main__':
    vec = np.array([0.9, -0.2, 0.1])
    print weights_in_subsection(vec)
