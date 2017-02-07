# Atomic coordinate translation and rotation functions;
# Functions of calculations of distance of atom pair, angle of three atoms and torsion of four atoms were also included.
#     -- by Lifeng Zhao, 02,4th,2009

import sys
from math import *
import numpy as np
from math import sqrt
from copy import deepcopy
class bcolors:
    """
    Usage: "  print bcolors.BLUE + 'text you want to print' + bcolors.ENDC"
    """
    PINK = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.PINK = ''
        self.BLUE = ''
        self.GREEN = ''
        self.YELLOW = ''
        self.RED = ''
        self.ENDC = ''

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


def middle(x12):
    """
    Calculate the center of atom pair.
    x12: [[x1,y1,z1],[x2,y2,z2]]
    """

    x = [0.0, 0.0, 0.0]
    for i in range(len(x12)):
        x = [ x[k]+ x12[i][k] for k in range(3)]
    x = [x[k]/len(x12) for k in range(3)]
    #dx = (x12[0][0]+ x12[1][0])/2.0
    #dy = (x12[0][1]+ x12[1][1])/2.0
    #dz = (x12[0][2]+ x12[1][2])/2.0
    #result = [dx,dy,dz]

    return x



def distance2(x12):
    """
    Calculate atom pair distance.
    x12: [[x1,y1,z1],[x2,y2,z2]]
    """

    dx = x12[0][0]-x12[1][0]
    dy = x12[0][1]-x12[1][1]
    dz = x12[0][2]-x12[1][2]
    result = dx*dx+dy*dy+dz*dz

    return result


def arraypointmatrix(array, matrix):
    """
    Point multiplication function for an array to a matrix.
    array:[x1,x2,...,xm]
    matrix:[[y11,y12,...,y1n],...,[ym1,ym2,...,ymn]]
    """

    array = np.array(array)
    matrix = np.array(matrix)
    TempRes = np.dot(array, matrix)

    return list(TempRes)


def SpecialAngle(x0,x1):
    """
    Calculate these two angles:
        alpha: vector (x0--x1) with y axle
        beta : dihedral between face 0 (vector (x0--x1) with y axle) and face 1 (yz face)
    return: radians
    """
    dx = [x1[0]-x0[0],x1[1]-x0[1],x1[2]-x0[2]]
    r  = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])

    # Angle between y axle and r:
    if x0[0]==x1[0] and x0[2]==x1[2]: alpha = 0.0
    else: alpha = acos(dx[1]/r)
    # Angle between face (with vector and y axle) and face YZ
    if x0[0]==x1[0]: beta = 0.0
    else:
        tempr = dx[0]/(r*sin(alpha))
        if tempr>=1.0 and tempr<1.00001: tempr = 1.0
        if tempr<=-1.0 and tempr>-1.00001: tempr = -1.0
        beta  = asin(tempr)
        #beta  = asin(dx[0]/(r*sin(alpha)))

        if   dx[2]<0 and dx[0]>0: beta =  pi-beta
        elif dx[2]<0 and dx[0]<0: beta = -pi-beta

    return [alpha,beta]


def translate(AtomCoord, dxyz):
    """
    AtomCoord: [x,y,z]
    dxyz: [dx,dy,dz]
    """

    tempx = AtomCoord[0] + dxyz[0]
    tempy = AtomCoord[1] + dxyz[1]
    tempz = AtomCoord[2] + dxyz[2]

    return [tempx,tempy,tempz]


def transmole(MoleCoord,dxyz):
    """
    MoleCoord: [[x0,y0,z0],[x1,y1,z1],...]
         dxyz: [dx,dy,dz]
    """

    tempmole = deepcopy(MoleCoord)
    for ii in range(len(tempmole)):
        tempmole[ii] = translate(tempmole[ii],dxyz)

    return tempmole

###############################
# Rotate a single atom around an "x", "y" or "z" by "angle" using *RIGHT* hand rule:
def rotate(AtomCoord, axle, angle):
    """
    Rotate a single atom around an "x", "y" or "z" by "angle" using *RIGHT* hand rule:
    AtomCoord: [x,y,z]
         axle: 'x','y' or 'z'
        angle: Unit: degree
    """

    arrayorign = [AtomCoord[0], AtomCoord[1], AtomCoord[2], 1.0]
    angle = radians(angle)

    if axle == 'x':
        matrix = [[1.0, 0.0, 0.0, 0.0],\
                  [0.0, cos(angle), sin(angle), 0.0],\
                  [0.0,-sin(angle), cos(angle), 0.0],\
                  [0.0, 0.0, 0.0, 1.0]]
    elif axle == 'y':
        matrix = [[cos(angle), 0.0,-sin(angle), 0.0],\
                  [0.0, 1.0, 0.0, 0.0],\
                  [sin(angle), 0.0, cos(angle), 0.0],\
                  [0.0, 0.0, 0.0, 1.0]]
    elif axle == 'z':
        matrix = [[cos(angle), sin(angle), 0.0, 0.0],\
                  [-sin(angle), cos(angle),0.0, 0.0],\
                  [0.0, 0.0, 1.0, 0.0],\
                  [0.0, 0.0, 0.0, 1.0]]
    else:
        err = 2
        print 'Error(%d)' % err
        sys.exit(err)

    arraynew = arraypointmatrix(arrayorign, matrix)
    #for ii in range(3): AtomCoord[ii] = arraynew[ii]
    temp = arraynew[:3]

    return temp


##############################
def rotmole(MoleCoord,point,axle,angle):
    """
    Rotate molecule according to "point" and "axle"(x,y or z) by "angle":
        MoleCoord: [[x0,y0,z0],[x1,y1,z1],...]
            point: [x,y,z]
    Unit: angle: degree
    """
    tempcoord = deepcopy(MoleCoord)
    dxyz = [-point[0],-point[1],-point[2]]

    rescoord = []
    for ii in tempcoord[:]:
        temp = deepcopy(translate(ii,dxyz))
        temp = deepcopy(rotate(temp,axle,angle))
        rescoord.append(deepcopy(translate(temp,point)))

    return rescoord


##############################
def rotmole2impose(MoleCoord,Vector):
    """
    Rotate molecule according to a vector which was constructed by two point r1 and r2 (r1->r2)
        and translate the molecule to r1.
        MoleCoord: [[x0,y0,z0],[x1,y1,z1],...]
               x0: [x,y,z]
               x1: [x,y,z]
    """
    angles = SpecialAngle(x0,x1)
    temp_mole = rotmole(MoleCoord,[0,0,0],'y',angles[0])
    temp_mole = rotmole(temp_mole,[0,0,0],'x',angles[1])

    return temp_mole


##############################
def rotmole1(MoleCoord, x0, x1, angle):
    """
    Rotate molecule according to a vector (x0->x1) by "angle"
    InputUnit: angle: degree; x: angstrom
    """
    ## the backward direction of the vector:
    dxyz = [-x0[0], -x0[1], -x0[2]]
    tempx1 = [x1[0]-x0[0], x1[1]-x0[1], x1[2]-x0[2]]

    temp_mole = deepcopy(MoleCoord)
    ## translate all atoms backward with -x0:
    temp_mole = transmole(temp_mole, dxyz)

    angles = SpecialAngle(x0,x1)
    tempangle = torsion([[0.0, 0.0, 1.0],
                         [0.0,-1.0, 0.0],
                         [0.0, 1.0, 0.0],
                         tempx1])
    angles = [degrees(angles[0]), tempangle]

    ## aline the vector (x0->x1) to "y":
    temp_mole = rotmole(temp_mole, [0.0, 0.0, 0.0], 'y', -angles[1])
    temp_mole = rotmole(temp_mole, [0.0, 0.0, 0.0], 'x', -angles[0])

    ## rotate "angle" according to "y":
    temp_mole = rotmole(temp_mole, [0.0, 0.0, 0.0], 'y', angle)

    ## move all atoms back to vector (x0->x1) from "y":
    temp_mole = rotmole(temp_mole, [0.0, 0.0, 0.0], 'x', angles[0])
    temp_mole = rotmole(temp_mole, [0.0, 0.0, 0.0], 'y', angles[1])

    ## translate back to x0 from the origin:
    temp_mole = transmole(temp_mole, x0)

    return temp_mole


def angle(x123):
    """
    Calculate angle of three atoms: 2 is the center atom.
    An angle value in degree will return.
    x123: [[x1,y1,z1],...,[x3,y3,z3]]
    """

    xa = x123[0]
    xb = x123[1]
    xc = x123[2]

    dxA = np.array(xb) - np.array(xc) #[xb[0]-xc[0],xb[1]-xc[1],xb[2]-xc[2]]
    dxB = np.array(xc) - np.array(xa) #[xc[0]-xa[0],xc[1]-xa[1],xc[2]-xa[2]]
    dxC = np.array(xa) - np.array(xb) #[xa[0]-xb[0],xa[1]-xb[1],xa[2]-xb[2]]

    rA = sqrt(dxA[0]*dxA[0]+dxA[1]*dxA[1]+dxA[2]*dxA[2])
    rB2 = dxB[0]*dxB[0]+dxB[1]*dxB[1]+dxB[2]*dxB[2]
    rC = sqrt(dxC[0]*dxC[0]+dxC[1]*dxC[1]+dxC[2]*dxC[2])

    result = acos((rA*rA + rC*rC - rB2)/(2*rA*rC))*180.0/pi

    return result



#####################################
def torsion(x1234):
    """
    Calculate torsion of four atoms:
    x1234: [[x1,y1,z1],...,[x4,y4,z4]]
    """
    temp = deepcopy(x1234)

    [alpha01, beta01] = SpecialAngle(temp[1],temp[2])
    for ii in temp[:]:
        ii[:] = rotate(ii,'y',-degrees(beta01))
        ii[:] = rotate(ii,'x',-degrees(alpha01))
        
    if (temp[3][1]==temp[2][1] and temp[3][2]==temp[2][2]) or \
       (temp[0][1]==temp[1][1] and temp[0][2]==temp[1][2]):
        print 'Error: linear triad.'
        sys.exit(1)

    [alpha01, beta01] = SpecialAngle(temp[2],temp[3])
    for ii in temp[:]:
        ii[:] = rotate(ii,'y',-beta01*180.0/pi)

    [alpha01, beta01] = SpecialAngle(temp[1],temp[0])
    result = -degrees(beta01)

#    temp = []
#    for ii in range(4):
#        temp0 = copy.deepcopy(x1234[ii][:])
#        temp1 = [str(ii),str(ii)] + temp0
#        temp.append(temp1)
#
#    TorsMole = molecule('torstemp',temp)
#    [alpha01, beta01] = SpecialAngle(x1234[1],x1234[2])
#
#    TorsMole = TorsMole.rotate('y',-beta01)
#    TorsMole = TorsMole.rotate('x',-alpha01)
#
#    # The 1234.5 is used for capping of NAC -C(=O)-N- bond
#    if x1234[3][0]!=1234.5:
#        x3 = TorsMole.getxyz(2)
#        x4 = TorsMole.getxyz(3)
#
#        [alpha02, beta02] = SpecialAngle(x3,x4)
#        TorsMole = TorsMole.rotate('y',-beta02)
#
#    x2 = TorsMole.getxyz(1)
#    x1 = TorsMole.getxyz(0)
#    
#    [alpha03, beta03] = SpecialAngle(x2,x1)
#    result = -beta03*180.0/pi

    return result

