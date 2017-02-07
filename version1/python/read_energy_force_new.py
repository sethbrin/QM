import os,sys
import numpy as np
from lib_atomcoords import distance
from grids_structures import DS
import gzip
import pdb

"""
    0.02:    The order of mole1 and mole2 is corrected:
        mole1: fixed molecule to origin
        mole2: probe molecule
        No.-1: mole1
        No. 0: mole2
    0.03:    No. -1: mole1
        No.  0: mole2
        All mole1 coords in dimers were not stored
        Energies of monomers are in a.u.
        Energies of dimers are deltaE and in kcal.
"""

__version__ = "0.03"

AU2KCAL = 23.0605*27.2116

## data density:
#R_NDX = [2.0, 2.2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
#PHI_angles = [i*15.0 for i in range(0, 7)]
#THETA_angles = {0: [float(i) for i in range(-165, 181, 15)],
#              1: [float(i) for i in range(-165, 181, 15)],
#              2: [float(i) for i in range(-165, 181, 15)],
#              3: [float(i) for i in range(-160, 181, 20)],
#              4: [float(i) for i in range(-150, 181, 30)],
#              5: [float(i) for i in range(-135, 181, 45)],
#              6: [0.0]}
#NTheta = {}
#nGrid = 0 # No. of grids for one dist
#for i in range(len(THETA_angles)):
#    NTheta[i] = len(THETA_angles[i])
#    nGrid += NTheta[i]
#nDist = len(R_NDX) # No. of dists
#nConf = 26 # No. of orientations for one grid
#nNorm = 2  # No. of confs for one orientation
#
### Convert degree to radians:
#for ii in range(len(PHI_angles)): PHI_angles[ii] *= D2R
#for k in THETA_angles.keys():
#        for ii in range(len(THETA_angles[k])):
#                THETA_angles[k][ii] *= D2R
#
#Phi_Ndx = {}
#temp = 0
#for i in range(len(PHI_angles)):
#    Phi_Ndx[i] = range(temp, temp+NTheta[i])
#    temp += NTheta[i]


class atom():
    def __init__(self, line):
        line = line.split()
        self.name = line[0]
        self.x = [float(line[2]), float(line[3]), float(line[4])]
        #self.f = [float(line[5]), float(line[6]), float(line[7])]

class molecule():
    def __init__(self):
        self.natoms = 0
        self.atoms = []

    def addatom(self, line):
        self.atoms.append(atom(line))
        self.natoms += 1

    def set_energy(self, e): self.E = e

SwapH = {0:0, 1:1, 2:2, 3:3, 4:5, 5:4}

class pair_configuration():
    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2
        self.fmole1 = []
        self.xmole2 = []
        self.fmole2 = []
        self.E = 0.0
        pass

    def set_energy(self, e):
        self.E = e

    def add_atom_force_in_mole1(self, line):
        return
        line = line.split()
        #temp = [float(line[5]),float(line[6]),float(line[7])]
        #self.fmole1.append( temp )

    def add_atom_in_mole2(self, line):
        line = line.split()
        temp = [float(line[2]),float(line[3]),float(line[4])]
        self.xmole2.append( temp )
        #temp = [float(line[5]),float(line[6]),float(line[7])]
        #self.fmole2.append( temp )
        
    def get_prop(self, name, w):
        if name=='E': return self.E
        else:
            import pdb
            pdb.set_trace()
            i = int(name[1:4])
            j = int(name[4:7])
            if w<0: 
                #print 'negative',name
                i = SwapH[i]
            if i<self.n1: return self.fmole1[i][j]
            else: return self.fmole2[i-self.n1][j]

    def writemole2(self):
        Ofile = open('temp3.inp', 'w')
        Ofile.write(' $DATA\ntemp\nC1\n')
        for ii in range(3):
            if ii==0: Ofile.write(' O  8.0 ')
            else: Ofile.write(' H  1.0 ')
            for j in range(3):
                Ofile.write(' %10.5f'%self.xmole2[ii][j])
            Ofile.write('\n')
        Ofile.write(' $END\n')
        Ofile.close()

    def set_name(self, name):
        self.name = name

class energy_force_data():
    def __init__(self, Iname, FragType):
        self.filename = Iname
        self.FT = FragType

        self.allcfg = []
        for i in range(DS[self.FT].nDist):
            self.allcfg.append([])
            for j in range(DS[self.FT].nGrid):
                self.allcfg[i].append([])
                for m in range(DS[self.FT].nConf):
                    self.allcfg[i][j].append([])
                    for k in range(DS[self.FT].nNorm):
                        self.allcfg[i][j][m].append('None')
        self.wtrcfg = []
        for m in range(DS[self.FT].nConf):
            self.wtrcfg.append([])
            for k in range(DS[self.FT].nNorm):
                self.wtrcfg[m].append('None')
        self.natoms1 = DS[self.FT].n1
        self.natoms2 = DS[self.FT].n2

        self.ReadFile()

    def ReadFile(self):
        datafile = gzip.open( self.filename, 'r' )
        ## Read Mole 1 (config No.-1):
        while 1:
            line = datafile.readline()
            if line=='':
                raise Exception, "No data found."
            if line.find('# No.')==0:
                line = line.split()
                if line[-1]!='-1':
                    raise Exception, "Mole 1 not found."

                for ii in range(3): line = datafile.readline()
                ## energy of mole 1:
                line = line.split()
                self.mole1 = molecule()
                self.mole1.set_energy( float(line[1]) )
                line = datafile.readline()
                ## coords of mole 1:
                while 1:
                    line = datafile.readline()
                    if line.find('# End of No.')==0: break
                    self.mole1.addatom(line)
                line = datafile.readline()
                break

        ## Read Mole 2 (config No. 0):
        while 1:
            line = datafile.readline()
            if line=='':
                raise Exception, "No data found."
            if line.find('# No.')==0:
                line = line.split()
                if line[-1]!='0':
                    raise Exception, "Mole 2 not found."

                for ii in range(3): line = datafile.readline()
                ## energy of mole 2:
                line = line.split()
                self.mole2 = molecule()
                self.mole2.set_energy( float(line[1]) )
                line = datafile.readline()
                ## coords of mole 2:
                while 1:
                    line = datafile.readline()
                    if line.find('# End of No.')==0: break
                    self.mole2.addatom(line)
                line = datafile.readline()
                break

        self.distndx = []
        self.distname = []
        self.distnum = []
        ## Read configs:
        while 1:
            line = datafile.readline() # No.
            line = datafile.readline() # file name
            if line=='': break
            
            name = line.split()[1]
            d,g,c = self._parse_name(line)
            if c>25:
                ic = 1
                c = c-26
            else:
                ic = 0
            datafile.readline()

            line = datafile.readline() # energy
            if line.find('Failed')>-1:
                for nn in range(3): datafile.readline()
                continue

            config = pair_configuration(self.natoms1, self.natoms2)
            config.set_name(name)
            line = line.split()

            E = float(line[1])
            config.set_energy( E )
            datafile.readline()
            ## coords & forces of mole1
            #for nn in range(self.natoms1):
            #    line = datafile.readline()
            #    #config.add_atom_force_in_mole1(line)

            ## coords & forces of mole2
            for nn in range(self.natoms2):
                line = datafile.readline()
                config.add_atom_in_mole2(line)
            H_E = self.H_energy(config)
            #H_E = 0.0
            #config.set_energy( (E-self.mole1.E-self.mole2.E)*AU2KCAL - H_E )
            config.set_energy( E - H_E )
            #print "%d %d %d %d %.2lf"%(d,g,c,ic,(E-H_E))

            datafile.readline() # End of This Config.
            datafile.readline()

            #if d ==2 and g==80 and c==25: print '#',E,H_E

            self.allcfg[d][g][c][ic] = config
            if self.wtrcfg[c][ic]=='None': self.wtrcfg[c][ic] = config

        datafile.close()

    def charge(self, q0,q1,r):
            return q0*q1/r*332.5

    def H_energy(self, config):
        res = 0.0
        if not DS[self.FT].IsHmm: return res

        for i in DS[self.FT].Hmm:
            try:eps = DS[self.FT].param['HC-OW'][0]
            except:
                print "Error:",DS[self.FT].param
                raise
            sig = DS[self.FT].param['HC-OW'][1]
            r = distance([self.mole1.atoms[i].x, config.xmole2[0]])
            res += DS[self.FT].lj(r, eps, sig)
            q0 = DS[self.FT].param['qHC']
            q1 = DS[self.FT].param['qOW']
            res += DS[self.FT].coulomb(q0, q1, r)
            #res += self.charge(DS[self.FT].param['qHC'], DS[self.FT].param['qOW'], r)

            eps = DS[self.FT].param['HC-HW'][0]
            sig = DS[self.FT].param['HC-HW'][1]
            for j in range(1,3):
                r = distance([self.mole1.atoms[i].x, config.xmole2[j]])
                res += DS[self.FT].lj(r, eps, sig)
                res += self.charge(DS[self.FT].param['qHC'], DS[self.FT].param['qHW'], r)
        return res

    def get_energy(self, i,j,k):
        return self.allcfg[i][j][k]

    def get_prop(self, i,j,si,ni,name, w, ehigh=100.0):
        if i<0:
            return ehigh
        #elif i>=DS[self.FT].nDist-1:
        elif i>DS[self.FT].nDist-1:
            return 0.0
        cfg = self.allcfg[i][j][si][ni]
        if cfg=='None': return ehigh
        else:
            return self.allcfg[i][j][si][ni].get_prop(name, w)


    def _parse_name(self, line):
        line = line.split()
        line = line[1].rsplit('.',2)
        line = line[0].split('_')
        d = DS[self.FT].R_NDX.index(float(line[1][1:]))
        g = int(line[2][1:])
        c = int(line[3][1:])
        return d,g,c

if __name__=='__main__':

    a = energy_force_data( sys.argv[1],'gln')
    print 'a.allcfg[ 6][61][25][1].xmole2',a.allcfg[ 6][61][24][1].xmole2,a.allcfg[6][61][24][1].name
    print a.get_prop(5,109,24,0,'E',0.8)
    print a.allcfg[5][109][24][0].get_prop('E', 1)
    print a.allcfg[2][80][25][0].get_prop('E', -0.8)

    ## test grids: correct.
    Ofile = open( 'grids_coords1.inp', 'w' )
    Ofile.write(' $DATA\n temp\nC1\n')
    for i in range(nGrid):
        Ofile.write(' H  1.0 ')
        for j in range(3):
            Ofile.write('%12.7f'%a.allcfg[6][i][12][1].xmole2[0][j])
        Ofile.write('\n')
    Ofile.write(' $end\n')
    Ofile.close()

    ## test rotations: correct
    from trans_rot_coords import translate, get_bisect_unit
    Ofile = open( 'rots_on_grid.inp', 'w')
    Ofile.write(' $data\n temp\nC1\n')
    for i in range(DS[a.FT].nConf):
        Ofile.write(' H  1.0 ')
        x0 = a.allcfg[4][5][i][1].xmole2[0]
        x1 = a.allcfg[4][5][i][1].xmole2[1]
        x2 = a.allcfg[4][5][i][1].xmole2[2]

        x1 = translate(x1, -np.array(x0))
        x2 = translate(x2, -np.array(x0))
        bis = get_bisect_unit(x1, x2) * 3.0
        for j in range(3):
            Ofile.write('%12.7f'%bis[j])
        Ofile.write('\n')
    Ofile.write(' $end\n')
    Ofile.close()

