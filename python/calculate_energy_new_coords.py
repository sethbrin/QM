import os,sys
from lib_atomcoords import distance
from trans_rot_coords import *
import numpy as np
from likely import Grid3DInterpolator
from read_energy_force_new import *
from grids_structures import DS,Grid_Quarts

__version__ = 0.03

AU2KCAL = 23.0605*27.2116


class new_atom():
    def __init__(self, line, ftype='gjf'):
        #self.add = {'gjf': self.addgjf,
        #        'gms': self.addinp,
        #        'pdb': self.addpdb}
        #self.add[ftype](line)
        if ftype=='gjf': self.addgjf(line)
        elif ftype=='gms': self.addinp(line)
        elif ftype=='pdb': self.addpdb(line)

        self.f = [0.0,0.0,0.0]

    def addgjf(self, line):
                line = line.split()
                self.a_nam = line[0]
                self.x = [float(line[1]), float(line[2]), float(line[3])]

    def addpdb(self, line):
        self.line = line
        self.i_atm = int(line[6:11])
        self.a_nam = line[11:16].strip()
        self.a_res = line[16:20].strip()
        self.a_chn = line[20:22].strip()
        self.i_res = int(line[22:26])
        self.x = []
        self.x.append(float(line[30:38]))
        self.x.append(float(line[38:46]))
        self.x.append(float(line[46:54]))

    def addinp(self, line):
        line = line.split()
        self.a_nam = line[0]
        self.x = [float(line[2]), float(line[3]), float(line[4])]

class coordinates():
    def __init__(self, n1, n2, FragType, name=''):
        ## n1,n2 is the number of atoms in mole1 and mole2:
        self.n1 = n1
        self.n2 = n2

        ## records of operations of translation and rotation:
        self.OperateNdx = []
        self.Operation = []

        ## fragment type:
        self.FT = FragType

        ## symmetry faces:
        self.symface = DS[self.FT].symface

        self.IsOriented = False
        self.facendx = {'yx':2, 'xy':2,
                        'yz':0, 'zy':0,
                        'zx':1, 'xz':1,
                'zarg':5,
                'zben':6}
        self.symm = [1,1,1]
        self.center = 0
        self.natoms = 0
        self.original_atoms = []

        self.name = name

    def addatom(self, line, ftype='pdb'):
        temp = new_atom(line, ftype)
        self.original_atoms.append(temp)
        self.natoms += 1

    def addpdbatom(self, line):
        self.original_atoms.append(new_atom(line, 'pdb'))
        self.natoms += 1

    def set_atom(self, i, atom):
        if i>=len(self.original_atoms):
            self.original_atoms.append( deepcopy(atom) )
            self.natoms += 1
        else: self.original_atoms[i] = deepcopy(atom)

    #def calt_com(self, lst):
    #    res = [0.0, 0.0, 0.0]
    #    lenlst = len(lst)
    #    for i in range(3):
    #        for j in range(lenlst):
    #            res[i] += self.original_atoms[lst[j]][i]
    #        res[i] /= lenlst
    #    self.com = res

    def MirrorAll(self):
        """
        According to the coords of the 1st atom in mole2.
        """
        for face in self.symface:
            fndx = self.facendx[face]
            ## ben:
            if fndx==6:
                ## Rotate to xz (+-)30.0:
                if self.ang2<0: tempang2 = self.ang2+2.0*np.pi
                else: tempang2 = self.ang2
                period = np.pi/3.0
                nrot = tempang2//period
                rrot = tempang2%period
                if (rrot-np.pi/6.0)>0.0:
                    nrot += 1
                    IsMirrorXZ = True
                else: IsMirrorXZ = False

                for i in range(self.n1, self.natoms):
                    self.atoms[i].x = rotate(self.atoms[i].x, [0,0,-1], nrot*period)

                ## mirror xz:
                if IsMirrorXZ:
                    for i in range(self.n1, self.natoms): self.atoms[i].x[1] *= -1
                continue
            elif fndx==5:
                ## Rotate to xz (+-)60.0:
                if self.ang2<0: tempang2 = self.ang2+2.0*np.pi
                else: tempang2 = self.ang2
                period = np.pi*0.66666666666667
                nrot = tempang2//period
                rrot = tempang2%period
                if (rrot-np.pi/3.0)>0.0:
                    nrot += 1
                    IsMirrorXZ = True
                else: IsMirrorXZ = False

                for i in range(self.n1, self.natoms):
                    self.atoms[i].x = rotate(self.atoms[i].x, [0,0,-1], nrot*period)

                ## mirror xz:
                if IsMirrorXZ:
                    for i in range(self.n1, self.natoms): self.atoms[i].x[1] *= -1
                continue

            #if self.atoms[self.n1].x[fndx]<0.0:
            if self.center2[fndx] < 0.0:
                self.symm[ fndx ] = -1
                for i in range(self.n1, self.natoms): 
                    self.atoms[i].x[fndx] *= -1
                    self.atoms[i].f[fndx] *= -1

        self._spherical_x()
        #print self.r, self.ang1*180/np.pi, self.ang2*180/np.pi
        #self.writegjf()

    def ImMirrorAll(self):
        for i in range(3):
            if self.symm[i]<0:
                self.symm[i] *= -1
                for j in range(self.natoms): 
                    self.atoms[j].x[fndx] *= -1
                    self.atoms[j].f[fndx] *= -1

    def ImMirrorForce(self):
        if self.symm[1]<0:
            for i in range(3):
                name1 = 'f%03d%03d'%(1,i)
                name2 = 'f%03d%03d'%(2,i)
                temp = self.properties[name1]
                self.properties[name1] = self.properties[name2]
                self.properties[name2] = temp
            
        for i in range(3):
            if self.symm[i]<0:
                ## print 'sym',i
                for j in range(self.natoms):
                    name = 'f%03d%03d'%(j,i)
                    self.properties[name] *= -1

    def ReorientToOrigin(self, cut=0.0000001):
        #t = self.original_atoms
        #from copy import deepcopy
        #a = deepcopy(t)
        #for i in range(len(a)):
        #    print t[i].x
        self.atoms = deepcopy(self.original_atoms)
        ## translate:
        ## For new coords that the bck is tuned to the situation that 
        ##   the oxygen is on x-axis.
        dvec = DS[self.FT].calt_dvec( self.atoms[0].x, self.atoms[1].x, self.atoms[2].x )

        for i in range(self.natoms):
            self.atoms[i].x = translate(self.atoms[i].x, dvec)
        self.OperateNdx.append(0)
        self.Operation.append(np.array(dvec))

        ## rotate oxygen to x:
        vec, ax0 = DS[self.FT].calt_vec1( self.atoms[0].x, self.atoms[1].x, self.atoms[2].x )
        ang = angle(vec, ax0)
        ax = get_normal(vec, ax0)
        if ax[0]==0.0 and ax[1]==0.0 and ax[2]==0.0: pass
        else:
            for i in range(self.natoms):
                self.atoms[i].x = rotate(self.atoms[i].x, ax, ang)
            self.OperateNdx.append(1)
            self.Operation.append([ax, ang])
        ## rotate C5 to y:
        vec, ax0 = DS[self.FT].calt_vec2( self.atoms[0].x, self.atoms[1].x, self.atoms[2].x )
        ang = angle(vec, ax0)

        ## Such a condition is used for 0 degree or 180 degree test:
        #if ax[0]==0.0 and ax[1]==0.0 and ax[2]==0.0: pass

        if abs(ang)<cut: pass
        else:
            if abs(ang-np.pi)<cut: ax = [1,0,0]
            else: ax = get_normal(vec, ax0)
            for i in range(self.natoms):
                self.atoms[i].x = rotate(self.atoms[i].x, ax, ang)
            self.OperateNdx.append(1)
            self.Operation.append([ax, ang])

        self.IsOriented = True
        self._spherical_x()

        #self.writegjf()

    def _spherical_x(self):
        """
        Calculate the coords in spherical coordination system for molecule 2.
        """

        ## the center of molecule 2 :i
        x = self.atoms[self.n1].x
        r = np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])

        ## phi   of principal coords:
        ang1 = np.pi*0.5 - np.arccos(x[2]/r)

        ## theta of principal coords (from -pi to pi):
        if abs(x[0])<0.000001: 
            if x[1]>0: ang2 = np.pi*0.5
            else:      ang2 = np.pi*1.5
        else:
            ang2 = np.arctan(x[1]/x[0])
            if x[0]<0: ang2 += np.pi
            elif x[1]<0: ang2 += np.pi*2

        self.r = r
        self.ang1 = ang1
        self.ang2 = ang2
        self.center2 = x

    def indexing(self):
        if not self.IsOriented: raise Exception, "Error: indexing beforce reorientation."

        r = self.r
        ang1 = self.ang1
        ang2 = self.ang2

        ## ndx of r:
        ir = 10001
        if r<DS[self.FT].R_NDX[0]: ir = -1
        else:
            for i in range(1,DS[self.FT].nDist):
                if r<=DS[self.FT].R_NDX[i]:
                    ir = i-1
                    break
        #print 'ir',ir

        if ir>10000:
            self.ir = ir
            self.ig = 0
            self.vbis = [0,0,0]
            self.vnrm = [0,0,0]
            return 10000,0,0
        elif ir<0:
            self.ir = ir
            self.ig = 0
            self.vbis = [0,0,0]
            self.vnrm = [0,0,0]
            return -1, 0,0
        #print "r=%.1f"%r, ir
        ## ndx of ang1 (Phi):
        if ang1<DS[self.FT].PHI_angles[0]: ih = -1
        for i in range(1, DS[self.FT].nPhi):
            if ang1<=DS[self.FT].PHI_angles[i]:
                ih = i-1
                #ang1 = ang1 - DS[self.FT].PHI_angles[ih]
                break
        #print 'ih', ih

        ## ndx_of_ang2 (Theta):
        ip = -1
        #print DS[self.FT].NTheta[ih]
        for i in range(1, DS[self.FT].NTheta[ih]):
            #print i,ang2, DS[self.FT].THETA_angles[ih][i]
            if ang2<=DS[self.FT].THETA_angles[ih][i]:
                ip = i-1
                #ang2 = ang2 - DS[self.FT].THETA_angles[ih][ip]
                break
        if ip==-1: ip = DS[self.FT].NTheta[ih]-1
        #print 'ip',ip

        ## ig is the index of all the grids:
        ig = 0
        for i in range(ih): ig += DS[self.FT].NTheta[i]
        ig += ip

        ## calculate the vectors of bisector and normal of mole2:
        a20 = self.atoms[self.n1].x
        a21 = self.atoms[self.n1+1].x
        a22 = self.atoms[self.n1+2].x
        a20 = np.array(a20)
        a21 = np.array(a21)
        a22 = np.array(a22)

        v0 = a21 - a20
        v1 = a22 - a20
        ## These two vectors must be unit vector:
        bisect = get_bisect_unit(v0,v1)
        normal = get_normal_unit(v0,v1)

        ## relative coords of mole2 in the cubic of the 8 corners:
        r -= DS[self.FT].R_NDX[ir]
        ang1 -= DS[self.FT].PHI_angles[ih]
        ang2 -= DS[self.FT].THETA_angles[ih][ip]
        xr = r/DS[self.FT].DR[ir]
        #print 'r,r0,a1,a2', r,DS[self.FT].DR[ir],ang1,ang2
        tr = ang1/(DS[self.FT].PHI_angles[ih+1] - DS[self.FT].PHI_angles[ih])
        if ip==DS[self.FT].NTheta[ih]-1:
            pr = ang2/(2*np.pi+DS[self.FT].THETA_angles[ih][0]-DS[self.FT].THETA_angles[ih][ip])
        else:
            pr = ang2/(DS[self.FT].THETA_angles[ih][ip+1]-DS[self.FT].THETA_angles[ih][ip])

        self.ir = ir
        self.ig = ig
        ## xr: distance
        ## tr: Phi angle (between xr to z axis)
        ## pr: Theta angle (between xr to xz face)
        self.rel_x = [xr,tr,pr]
        self.vbis = bisect
        self.vnrm = normal

        #print 'r, ir, g, ig, xr, tr, pr'
        #print self.ir, self.ig, 'relative x', xr, tr, pr
        #print self.vbis, self.vnrm

    def SetPdbMoles(self, lines1, lines2):
        if self.natoms!=0: raise Exception, "Error, add atom to existed mole."
        for line in lines1: self.addpdbatom(line)
        for line in lines2: self.addpdbatom(line)

    def readgjf(self,fname):
        Ifile = open(fname)
        while 1:
            line = Ifile.readline()
            if line.find('0 1')==0: break
        while 1:
            line = Ifile.readline()
            if line.strip()== '': break
            self.addatom(line)
        Ifile.close()

    def ReadGmsInp(self, Iname):
        Ifile = open(Iname)
        while 1:
            line = Ifile.readline()
            if line=='': break
            if line.find(' $DATA')==0:
                Ifile.readline()
                Ifile.readline()
                while 1:
                    line = Ifile.readline()
                    if line.find(' $END')==0: break
                    self.addatom(line, 'gms')
                break
        Ifile.close()

    def writegjf(self):
        Ofile = open( self.name+'.gjf', 'w' )
        Ofile.write('#am1\n\ntemp\n\n0 1\n')
        for ii in range(self.natoms):
            Ofile.write('%5s'% self.atoms[ii].a_nam)
            for j in range(3): Ofile.write('%10.4f'%self.atoms[ii].x[j])
            Ofile.write('\n')
        Ofile.write('\n\n')
        Ofile.close()

    def writeinp(self):
        tempgrids = """ H  1.0    4.0000000   0.0000000   0.0000000
 H  1.0    3.8637033   1.0352762   0.0000000
 H  1.0    3.4641016   2.0000000   0.0000000
 H  1.0    2.8284271   2.8284271   0.0000000
 H  1.0    2.0000000   3.4641016   0.0000000
 H  1.0    1.0352762   3.8637033   0.0000000
 H  1.0    0.0000000   4.0000000   0.0000000
 H  1.0   -1.0352762   3.8637033   0.0000000
 H  1.0   -2.0000000   3.4641016   0.0000000
 H  1.0   -2.8284271   2.8284271   0.0000000
 H  1.0   -3.4641016   2.0000000   0.0000000
 H  1.0   -3.8637033   1.0352762   0.0000000
 H  1.0   -4.0000000   0.0000000   0.0000000
 H  1.0    3.8637033   0.0000000   1.0352762
 H  1.0    3.7320508   1.0000000   1.0352762
 H  1.0    3.3460652   1.9318517   1.0352762
 H  1.0    2.7320508   2.7320508   1.0352762
 H  1.0    1.9318517   3.3460652   1.0352762
 H  1.0    1.0000000   3.7320508   1.0352762
 H  1.0    0.0000000   3.8637033   1.0352762
 H  1.0   -1.0000000   3.7320508   1.0352762
 H  1.0   -1.9318517   3.3460652   1.0352762
 H  1.0   -2.7320508   2.7320508   1.0352762
 H  1.0   -3.3460652   1.9318517   1.0352762
 H  1.0   -3.7320508   1.0000000   1.0352762
 H  1.0   -3.8637033   0.0000000   1.0352762
 H  1.0    3.4641016   0.0000000   2.0000000
 H  1.0    3.3460652   0.8965755   2.0000000
 H  1.0    3.0000000   1.7320508   2.0000000
 H  1.0    2.4494897   2.4494897   2.0000000
 H  1.0    1.7320508   3.0000000   2.0000000
 H  1.0    0.8965755   3.3460652   2.0000000
 H  1.0    0.0000000   3.4641016   2.0000000
 H  1.0   -0.8965755   3.3460652   2.0000000
 H  1.0   -1.7320508   3.0000000   2.0000000
 H  1.0   -2.4494897   2.4494897   2.0000000
 H  1.0   -3.0000000   1.7320508   2.0000000
 H  1.0   -3.3460652   0.8965755   2.0000000
 H  1.0   -3.4641016   0.0000000   2.0000000
 H  1.0    2.8284271   0.0000000   2.8284271
 H  1.0    2.6578521   0.9673791   2.8284271
 H  1.0    2.1667009   1.8180779   2.8284271
 H  1.0    1.4142136   2.4494897   2.8284271
 H  1.0    0.4911512   2.7854570   2.8284271
 H  1.0   -0.4911512   2.7854570   2.8284271
 H  1.0   -1.4142136   2.4494897   2.8284271
 H  1.0   -2.1667009   1.8180779   2.8284271
 H  1.0   -2.6578521   0.9673791   2.8284271
 H  1.0   -2.8284271   0.0000000   2.8284271
 H  1.0    2.0000000   0.0000000   3.4641016
 H  1.0    1.7320508   1.0000000   3.4641016
 H  1.0    1.0000000   1.7320508   3.4641016
 H  1.0    0.0000000   2.0000000   3.4641016
 H  1.0   -1.0000000   1.7320508   3.4641016
 H  1.0   -1.7320508   1.0000000   3.4641016
 H  1.0   -2.0000000   0.0000000   3.4641016
 H  1.0    1.0352762   0.0000000   3.8637033
 H  1.0    0.7320508   0.7320508   3.8637033
 H  1.0    0.0000000   1.0352762   3.8637033
 H  1.0   -0.7320508   0.7320508   3.8637033
 H  1.0   -1.0352762   0.0000000   3.8637033
 H  1.0    0.0000000   0.0000000   4.0000000
"""
        Ofile = open( 'temp2.inp', 'w' )
        Ofile.write(' $DATA\ntemp\nC1\n')
        Ofile.write(tempgrids)
        for ii in range(self.natoms):
            if ii == 0 or ii==3:
                Ofile.write(' O  8.0 ')
            else: Ofile.write(' H  1.0 ')
            for j in range(3): Ofile.write(' %10.5f'%self.atoms[ii].x[j])
            Ofile.write('\n')
        Ofile.write(' $END\n')
        Ofile.close()

    def writemole2origin(self):
        temprots = """ H  1.0    3.0000000   0.0000000   0.0000000
 H  1.0    2.1213203   0.0000000   2.1213203
 H  1.0   -0.0000000   0.0000000   3.0000000
 H  1.0   -2.1213203   0.0000000   2.1213203
 H  1.0   -3.0000000   0.0000000   0.0000000
 H  1.0   -2.1213203   0.0000000  -2.1213203
 H  1.0   -0.0000000   0.0000000  -3.0000000
 H  1.0    2.1213203   0.0000000  -2.1213203
 H  1.0    2.1213203   2.1213203   0.0000000
 H  1.0    1.5000000   2.1213204   1.5000000
 H  1.0   -0.0000000   2.1213203   2.1213203
 H  1.0   -1.5000000   2.1213204   1.5000000
 H  1.0   -2.1213203   2.1213203   0.0000000
 H  1.0   -1.5000000   2.1213204  -1.5000000
 H  1.0   -0.0000000   2.1213203  -2.1213203
 H  1.0    1.5000000   2.1213204  -1.5000000
 H  1.0    2.1213203  -2.1213204   0.0000000
 H  1.0    1.5000000  -2.1213204   1.5000000
 H  1.0   -0.0000000  -2.1213204   2.1213203
 H  1.0   -1.5000000  -2.1213204   1.5000000
 H  1.0   -2.1213203  -2.1213204   0.0000000
 H  1.0   -1.5000000  -2.1213204  -1.5000000
 H  1.0   -0.0000000  -2.1213204  -2.1213203
 H  1.0    1.5000000  -2.1213204  -1.5000000
 H  1.0    0.0000000   3.0000000   0.0000000
 H  1.0    0.0000000  -3.0000000   0.0000000
"""
        Ofile = open('mole2origin.inp', 'w')
        Ofile.write(' $data\n temp\nC1\n')
        Ofile.write(temprots)
        dx = -np.array(self.atoms[3].x)

        for i in range(3,6):
            x = self.atoms[i].x
            x = translate(x, dx)
            if i==3: Ofile.write(' O  8.0 ')
            else: Ofile.write(' H  1.0 ')
            for j in range(3):
                Ofile.write(' %12.8f'%x[j])
            Ofile.write('\n')
        tempb = self.vbis*3.0
        tempn = self.vnrm*3.0
        Ofile.write(' H  1.0 ')
        for i in range(3): Ofile.write(' %12.8f'%tempb[i])
        Ofile.write('\n')
        Ofile.write(' H  1.0 ')
        for i in range(3): Ofile.write(' %12.8f'%tempn[i])
        Ofile.write('\n')

        Ofile.write(' $end\n')
        Ofile.close()

    def calt_conf_energy(self, allconfigs, IsForce=False, ehigh=100.0):
        ri = self.ir
        if ri>100:
            self.properties = {'E':0.0}
            return
        elif ri<0:
            self.properties = {'E':ehigh}
            return

        gi = self.ig
        relative_x = self.rel_x
        bisv = self.vbis
        nrmv = self.vnrm
        
        #print 'ri', ri, gi, relative_x, bisv, nrmv
        cart_ndx, grids_sub_ndx, wghx, wghy = weights_in_subsection( bisv )
        #print bisv
        #print 'cart_ndx', cart_ndx, grids_sub_ndx, wghx, wghy

        properties = {'E':[]}
        propnames = ['E']

        if IsForce:
            for i in range(self.natoms):
                for j in range(3):
                    name = 'f%03d%03d'%(i,j)
                    properties[name] = []
                    propnames.append(name)
        properties = {'E':[]}
        propnames = ['E']
        tempprop = deepcopy(properties)

        ndx = 0

        ## tempf ('tempbox.gjf') is used for testing the grids and orientations:
        #tempf = open('tempbox.gjf', 'w')
        #IsWrite = True
        ## end of tempf 1/5

        for i in range(ri, ri+2):
            #print 'i',i
            for j in Grid_Quarts[self.FT][gi]:
                #print 'j',j
                ## Ei[i]: average energy of the two ortho-normal pair of each of
                #    the four grids in a sub-section.

                prop = deepcopy(tempprop)
                for ni in grids_sub_ndx:
                    #xconf0 = allconfigs.wtrcfg[ni][0].xmole2
                    #xconf1 = allconfigs.wtrcfg[ni][1].xmole2
                    inpfile0 = 'r%3.1f/tempconf_d%3.1f_g%02d_c%02d.inp'%(DS[self.FT].R_NDX[i],DS[self.FT].R_NDX[i],j,ni)
                    inpfile1 = 'r%3.1f/tempconf_d%3.1f_g%02d_c%02d.inp'%(DS[self.FT].R_NDX[i],DS[self.FT].R_NDX[i],j,ni+26)
                        
                    #print i, j, ni, inpfile0, inpfile1
                    #if os.path.exists(inpfile0) and os.path.exists(inpfile1):
                    xconf0 = allconfigs.allcfg[i][j][ni][0].xmole2
                    xconf1 = allconfigs.allcfg[i][j][ni][1].xmole2
                    w0, w1 = weights_for_2_configs( nrmv, xconf0, xconf1 )
                    #print allconfigs.get_prop(i,j,ni,0,'E',w0, ehigh=ehigh), w0, allconfigs.get_prop(i,j,ni,1,'E',w1, ehigh=ehigh), w1
                    #elif not os.path.exists(inpfile0) and os.path.exists(inpfile1):
                    #    print "%s does not exist"%inpfile0
                    #    w0 = 0.0
                    #    w1 = 1.0
                    #    print allconfigs.get_prop(i,j,ni,0,'E',w0, ehigh=ehigh), allconfigs.get_prop(i,j,ni,1,'E',w1, ehigh=ehigh)
                    #elif os.path.exists(inpfile0) and not os.path.exists(inpfile1):
                    #    print "%s does not exist"%inpfile1
                    #    w0 = 1.0
                    #    w1 = 0.0
                    #    print allconfigs.get_prop(i,j,ni,0,'E',w0, ehigh=ehigh), allconfigs.get_prop(i,j,ni,1,'E',w1, ehigh=ehigh)
                    #else:
                    #    print "%s does not exist"%inpfile0
                    #    print "%s does not exist"%inpfile1
                    #    w0 = 0.0
                    #    w1 = 0.0
                    #    print allconfigs.get_prop(i,j,ni,0,'E',w0, ehigh=ehigh), allconfigs.get_prop(i,j,ni,1,'E',w1, ehigh=ehigh)
                        
                    for pp in propnames:
                        p0 = allconfigs.get_prop(i,j,ni,0,pp,w0, ehigh=ehigh) #allconf[i][j][ni][0].get_prop(pp, w0)
                        #except: p0 = ehigh
                        p1 = allconfigs.get_prop(i,j,ni,1,pp,w1, ehigh=ehigh) #allconf[i][j][ni][1].get_prop(pp, w1)
                        #except: p1 = ehigh
                        #if p0==ehigh: p = p1
                        #elif p1 == ehigh: p = p0
                        #else: p = p1*abs(w1) + p0*abs(w0)
                        p = p1*abs(w1) + p0*abs(w0)
                        #print i,j,ni,w0,p0,p1,p
                        prop[pp].append(p)
                    #if i==2 and j==94: print ni, allconf[i][j][ni][0].a_nam, 'p0,p1', p0, p1

                ## tempf ('tempbox.gjf') is used for testing the grids and orientations:
                #IsWrite = False
                ## end of tempf 4/5

                ## Esub: bilinear interpolated energy of the sub-section.
                #print prop
                for pp in propnames:
                    psub = bilinear(prop[pp][0], prop[pp][1], prop[pp][2], prop[pp][3], wghx, wghy)
                    #print prop,wghx,wghy,psub
                    #if i==1 and j==77: print 'bilinear',j, prop[pp][0], prop[pp][1], prop[pp][2], prop[pp][3], '@', psub, wghx, wghy
                    properties[pp].append(psub)
    
        ## tempf ('tempbox.gjf') is used for testing the grids and orientations:
        #tempf.close()
        ## end of tempf 5/5
        #print '@', properties[pp]

        ## construct cubic matrix for tricubic interpolation:
        self.properties = {}
        for pp in propnames:
            pCubic = [[[properties[pp][0], properties[pp][1]],
                   [properties[pp][2], properties[pp][3]]],
                  [[properties[pp][4], properties[pp][5]],
                   [properties[pp][6], properties[pp][7]]]]
    
            #print np.array(pCubic)
            #print relative_x
            pCubic = np.array(pCubic, dtype=float, order='F')
            grids8 = Grid3DInterpolator( pCubic, 1.,1.,1. )
            self.properties[pp] = grids8.test( relative_x[0], relative_x[1], relative_x[2] )
            #self.properties[pp] = grids8.test( 1.0, 0.0, 1.0 )
            #print pp,self.properties[pp]
        #if self.FT=='gln': print pCubic
        #print self.properties[pp]

        if IsForce: self.ImMirrorForce()

    def write_wtr_gjf(self, x, Oname, ty='n'):
        if ty=='n': 
            Ofile = open(Oname, 'w')
            Ofile.write('#am1\n\ntemp\n\n0 1\n')
        elif ty=='a': Ofile = open(Oname, 'a')
        for i in range(3):
            if i == 0: Ofile.write(' O  ')
            else: Ofile.write(' H  ')
            for j in range(3):
                Ofile.write(' %12.8f'%x[i][j])
            Ofile.write('\n')
        if ty=='a':
            Ofile.write('\n')
        Ofile.close()

    def get_interp_force(self, i, j):
        return self.properties[ 'f%03d%03d'%(i,j) ]
        if i<3: return self.properties[ 'f%03d%03d'%(i,j) ]
        else: return self.properties[ 'f%03d%03d'%(i-3,j) ]

    def get_interp_energy(self):
        return self.properties['E']

