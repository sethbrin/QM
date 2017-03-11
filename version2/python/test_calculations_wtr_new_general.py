import os,sys
import numpy as np
from read_energy_force_new import *
from grids_structures_general import DS,Grid_Quarts
from calculate_energy_new_coords_general import coordinates, new_atom
from trans_rot_coords import get_normal

AU2KCAL = 23.0605*27.2116
HperB2toque = 1185.82 # 1Hartree/Bohr = 1185.82 kcal/mol/Angstrom
tMass = [15.999, 1.008, 1.008]


class QM_test_interpolation():
    def __init__(self, resname):
        self.FT = resname

    def get_com(self,coords):
        x = [0,0,0]
        totalM = 0
        for i in range(len(coords)):
            x = [ x[k]+ coords[i][k]*tMass[i] for k in range(3)]
            totalM += tMass[i]
        x = [x[k]/totalM for k in range(3)]
        return x


    def extract_coordtxt(self, inpfile):
        coords = []
        Inf = open(inpfile,'r')
        while 1:
            line = Inf.readline()
            if line.strip() == '': break
            if line.startswith(' $data') or line.startswith(' $DATA'):
                for i in range(2): Inf.readline()
                break

        while 1:
            line = Inf.readline()
            if line.strip() == '' or line.startswith(' $end') or line.startswith(' $END'): break
            line = line.strip()
            coords.append(' '+line)

        Inf.close()
        return coords


    def run_compare_polycen(self, allconfig):
        ## atom indices in residues for fragment locating:
        ## Start from 1 !!!
        aa_ndx = [1,2,3]
        prob_ndx = [4,5,6]

        ## Monomer energies:
        E1 = allconfig.mole1.E
        E2 = allconfig.mole2.E

        ## filenames to be tested:
        temp = os.popen('ls random3/???.pdb.inp')
        fnames = []
        while 1:
            line = temp.readline()
            if line=='': break
            fnames.append(line.strip())
        temp.close()

        ## Read the QM energies:
        energies = {}

        ## Read coords of the .inp files:
        Ofile = open('wtr01.txt', 'w')
        fout1 = open('wtr_wtr_force.txt','w')
        fout1.write('# name  Fx  Fy  Fz\n')
        fout2 = open('wtr_wtr_torque.txt','w')
        fout2.write('# name  Tx  Ty  Tz\n')
        for name in fnames:
            #print name
            # extract coordinates
            lines = self.extract_coordtxt(name)

            interp = coordinates(3,3, self.FT, name+('%02d'%aa_ndx[0]))
            ## Read center molecule atoms:
            for iline in aa_ndx:
                interp.addatom(lines[iline-1], 'gms')

            ## Read probe:
            for iline in prob_ndx:
                interp.addatom(lines[iline-1], 'gms')

            ## Reorient and interpolation:
            interp.ReorientToOrigin()
            interp.MirrorAll()
            try:interp.indexing_auto3()
            except:
                #print name
                raise
            interp.calt_conf_energy( allconfig )
            e_interp = interp.get_interp_energy()
            dist = interp.r
            interp.reverse_force_toque ()
            force = interp.get_interp_force()
            torque = interp.get_interp_torque()

            Ofile.write('%s %12.7f %.3f\n'%(name, e_interp, dist))
            Ofile.flush()
            fout1.write('%s %12.7f %12.7f %12.7f\n'%(name, force[0],force[1],force[2]))
            fout1.flush()
            fout2.write('%s %12.7f %12.7f %12.7f\n'%(name, torque[0],torque[1],torque[2]))
            fout2.flush()

        Ofile.close()
        fout1.close()
        fout2.close()


if __name__=='__main__':
    resname = 'wtr'
    allconf = energy_force_data('Dimer_deltaEForceTorque_data_mp2_wtr_wtr.txt.gz', resname)

    test = QM_test_interpolation(resname)
    test.run_compare_polycen(allconf)

