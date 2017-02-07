import os,sys
import numpy as np
from trans_rot_coords import distance
from grids_structures_general import DS
import gzip


AU2KCAL = 23.0605*27.2116

class atom():
    def __init__(self, line):
        line = line.split()
        self.name = line[0]
        self.x = [float(line[2]), float(line[3]), float(line[4])]

class molecule():
    def __init__(self):
        self.natoms = 0
        self.atoms = []

    def addatom(self, line):
        self.atoms.append(atom(line))
        self.natoms += 1

    def set_energy(self, e): self.E = e


class pair_configuration():
    def __init__(self, n1, n2):
        self.n1 = n1
        self.n2 = n2
        self.fmole1 = []
        self.xmole2 = []
        self.E = 0.0
        pass

    def set_energy(self, e):
        self.E = e

    def set_force(self, ff):
        self.Fx = ff[0]
        self.Fy = ff[1]
        self.Fz = ff[2]

    def set_torque(self,tq):
        self.Tx = tq[0]
        self.Ty = tq[1]
        self.Tz = tq[2]

    def add_atom_in_mole2(self, line):
        line = line.split()
        temp = [float(line[2]),float(line[3]),float(line[4])]
        self.xmole2.append( temp )

    def get_prop(self, name, w):
        if name=='E': return self.E
        elif name == 'Fx': return self.Fx
        elif name == 'Fy': return self.Fy
        elif name == 'Fz': return self.Fz
        elif name == 'Tx': return self.Tx
        elif name == 'Ty': return self.Ty
        elif name == 'Tz': return self.Tz
        else: return 'NONE'

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
                for m in range(DS[self.FT].nConf[i]):
                    self.allcfg[i][j].append([])
                    for k in range(DS[self.FT].nNorm[i]):
                        self.allcfg[i][j][m].append('None')

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
                for ii in range(3): line = datafile.readline()
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
                for ii in range(3): line = datafile.readline()
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
            ic = c//DS[self.FT].nConf[d]
            c = c%DS[self.FT].nConf[d]
            datafile.readline()

            line = datafile.readline() # energy
            if line.find('Failed')>-1:
                for nn in range(3): datafile.readline()
                continue

            config = pair_configuration(self.natoms1, self.natoms2)
            config.set_name(name)

            # energy
            line = line.split()
            E = float(line[1])
            config.set_energy( E )
            datafile.readline()

            # force
            line = datafile.readline()
            tmp = line.strip().split()
            ff = [float(x) for x in tmp[1:4]]
            config.set_force( ff )

            # torque
            line = datafile.readline()
            tmp = line.strip().split()
            tq = [float(x) for x in tmp[1:4]]
            config.set_torque( tq )

            ## coords of mole2
            for nn in range(self.natoms2):
                line = datafile.readline()
                config.add_atom_in_mole2(line)

            datafile.readline() # End of This Config.
            datafile.readline()

            self.allcfg[d][g][c][ic] = config

        datafile.close()


    def get_energy(self, i,j,k):
        return self.allcfg[i][j][k]

    def get_prop(self, i,j,si,ni,name, w, ehigh=100.0):
        if i<0:
            return ehigh
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

