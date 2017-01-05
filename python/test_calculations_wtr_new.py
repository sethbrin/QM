import os,sys
from lib_atomcoords import distance
import numpy as np
from read_energy_force_new import *
from grids_structures import DS,Grid_Quarts
from calculate_energy_new_coords import coordinates, new_atom

AU2KCAL = 23.0605*27.2116


class QM_test_interpolation():
    def __init__(self, resname):
        self.FT = resname

        def charge(self, q0,q1,r):
            return q0*q1/r*332.5

        def H_energy(self, config1, config2):
            res = 0.0
            if not DS[self.FT].IsHmm: return res
            # for i in range(len(config1)):
            #         eps = DS[self.FT].param['HC-OW'][0]
            #         sig = DS[self.FT].param['HC-OW'][1]
            #         r = distance([config1[i].x, config2[0].x])
            #         res += DS[self.FT].lj(r, eps, sig)
            #         res += self.charge(DS[self.FT].param['qHC'], DS[self.FT].param['qOW'], r)

            #         eps = DS[self.FT].param['HC-HW'][0]
            #         sig = DS[self.FT].param['HC-HW'][1]
            #         for j in range(1,3):
            #                 r = distance([config1[i].x, config2[j].x])
            #                 res += DS[self.FT].lj(r, eps, sig)
            #                 res += self.charge(DS[self.FT].param['qHC'], DS[self.FT].param['qHW'], r)

            # return res
    def run_compare_polycen(self, allconfig):
        ## atom indices in residues for fragment locating:
        ## Start from 1 !!!
                #aa_ndx = [[1,3,5]]
        aa_ndx = [[1,2,3]]
        C_ndx = [1] #old, not correct
        prob_ndx = [4,5,6]
        H_ndx = [2,3,4] # old, not correct

        ## Monomer energies:
        E1 = allconfig.mole1.E
        E2 = allconfig.mole2.E

        ## filenames to be tested:
        temp = os.popen('ls random/test????.inp')
        fnames = []
        while 1:
            line = temp.readline()
            if line=='': break
            fnames.append(line.strip())
        temp.close()

        ## Read the QM energies:
        # energies = {}
        # for name in fnames:
        #     #temp = os.popen('grep "TOTAL ENERGY =" '+name+'.log')
        #     temp = os.popen('grep "TOTAL ENERGY =" '+name+'.log')
        #     e = temp.read()
        #                       print e
        #     print name

        #     temp.close()
        #     e = e.split()
        #                       print e
        #     try:energies[name] = (float(e[-1]) - E1 - E2) * AU2KCAL
        #     except:
        #         print energies, e, name
        #         raise

        ## Read coords of the .inp files:
        IsCorrect = 0
        alltests = {}
        all_dist = {}
        for name in fnames:
            alltests[name] = 0.0
            all_dist[name] = []  
            #ecurr = 0.0
            #icoord = -1
            lines = []
            ## "icoord": beginning location of the coords in the file:
            Ifile = open(name, 'r')
            while 1:
                line = Ifile.readline()
                if line=='': raise Exception, "no coords found"
                if line.find(' $DATA')==0:
                    for i in range(2): Ifile.readline()
                    icoord = Ifile.tell()
                    break
            while 1:
                line = Ifile.readline()
                if line.strip() == '' or line.startswith(' $end') or line.startswith(' $END'): break
                line = line.strip()
                lines.append(line)
                #ce shi
                # fil_object=open('line.txt','w')
                # for line in lines:
                #     fil_object.write('%s\n\n' % line)
                # fil_object.close()    

            Ifile.close()

            ## For H force field energy test:
            # H_atoms = []
            # for i in H_ndx:
            #     H_atoms.append(new_atom(lines[i-1], 'gms'))

            # prob_atoms = []
            # for i in prob_ndx:
            #     prob_atoms.append(new_atom(lines[i-1],'gms'))


            #H_E = self.H_energy(H_atoms, prob_atoms)
            # print name,energies[name]
            original_atoms = []
            for mndx in aa_ndx:
                interp = coordinates(3,3, self.FT, name+('%02d'%mndx[0]))
                ## Read center molecule atoms:
                for iline in mndx:
                    #print iline,lines[iline-1]
                    interp.addatom(lines[iline-1], 'gms')
                #for i in range(len(mndx)):
                    #print interp.original_atoms[i].a_nam, interp.original_atoms[i].x
                
                ## Read probe:
                for iline in prob_ndx:
                    #print iline, lines[iline-1]
                    interp.addatom(lines[iline-1], 'gms')
                #for i in range(len(prob_ndx)):
                                #        print interp.original_atoms[i+len(mndx)].a_nam, interp.original_atoms[i+len(mndx)].x

                #from copy import deepcopy
                #tt = deepcopy(original_atoms)
                #class test():
                #    def __init__(self,x):
                #        self.x = x*x
                #tt = []
                #for ii in range(5):
                #    temp = test(ii)
                #    tt.append(temp)
                #ttt = deepcopy(tt)    
                #sys.exit()
                ## Subtract overlap CH3:


                # if IsCorrect:
                #     for i in C_ndx:
                #         line = lines[i-1]
                #         line = line.split()
                #         x1 = []
                #         for j in range(2,5):
                #             x1.append(float(line[j]))
                #         alltests[name] -= lennard(x1,interp.atoms[3].x)

                ## Reorient and interpolation:
                interp.ReorientToOrigin()
                interp.MirrorAll()
                try:interp.indexing()
                except:
                    print name
                    raise
                interp.calt_conf_energy( allconfig )
                tempe = interp.get_interp_energy()
                #print mndx, tempe
                alltests[name] += tempe
                all_dist[name].append(interp.r)

            #alltests[name] += H_E

            Ifile.close()

        Ofile = open('wtr01.txt', 'w')
        for name in fnames:
            #if all_dist[name][0] <= 6.0 or all_dist[name][0] >= 8.0: continue
            Ofile.write('%s %12.7f %.3f\n'%(name,alltests[name],all_dist[name][0]))
            #if name.find('0391')>-1: print energies[name], alltests[name]
        Ofile.close()

# def lennard(x1,x2):
#     Csig = 2.7
#     Ceps = 0.11128/4.184
#     Osig = 3.15061
#     Oeps = 0.636386/4.184
#     sig = np.sqrt(Csig*Osig)
#     eps = np.sqrt(Ceps*Oeps)
#     r = distance([x1,x2])
#     tr = (sig/r)**6.0
#     res = 4.0*eps*(tr*tr-tr)
#     return res

class QMInterpolation():
  def __init__(self, ftype):
    self.ftype = ftype
    cur_dir = os.path.split(os.path.realpath(__file__))[0]
    database_name = 'Dimer_deltaE_data_mp2_wtr_wtr.txt.gz'
    database_name = cur_dir + os.sep + database_name
    self.allconf = energy_force_data(database_name, ftype)

  #
  # ('O', [8.0,0.10293690, 0.11093940, 0.07435482],
  #  'H', [8.0,0.10293690, 0.11093940, 0.07435482],
  #  'H', [8.0,0.10293690, 0.11093940, 0.07435482])
  def calculate(self, lhs, rhs):
    interp = coordinates(3,3, self.ftype, 'test')
    for item in lhs:
      interp.addgmsatom(item[0], item[1])
    for item in rhs:
      interp.addgmsatom(item[0], item[1])

    interp.ReorientToOrigin()
    interp.MirrorAll()
    try:interp.indexing()
    except:
      raise
    interp.calt_conf_energy(self.allconf)
    tempe = interp.get_interp_energy()
    return (tempe, interp.r)

def test():
  qm_interpolation = QMInterpolation('wtr')
  lhs = [('O', [8.0, 0.10293690, 0.11093940, 0.07435482]),
         ('H', [1.0, 0.73097151, 0.83572755, 0.02040124]),
         ('H', [1.0, -0.22532655, 0.13689786, 0.97669930])]
  rhs = [('O', [8.0, 2.85534030, 0.60573860, -1.03993027]),
         ('H', [1.0, 1.93176673, 0.51127346, -0.79346626]),
         ('H', [1.0, 3.16776562, -0.29302760, -1.17133059])]
  print qm_interpolation.calculate(lhs, rhs)

if __name__=='__main__':
    test()
    import sys
    sys.exit()
    resname = 'wtr'
    allconf = energy_force_data('Dimer_deltaE_data_mp2_wtr_wtr.txt.gz', resname)
    import pickle

    #f = open('water.db', 'wb')
    #pickle.dump(allconf, f)
    #f = open('water.db', 'rb')
    #allconfig = pickle.load(f)
    
    test = QM_test_interpolation(resname)
    test.run_compare_polycen(allconf)

