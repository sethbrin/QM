import numpy as np
from trans_rot_coords import get_normal, get_unit


D2R = 3.14159265358979/180.0
Mass_table = {'O':15.999, 'H':1.008}


## Indices of one square:
Grid_Quarts_wtr = { 0: [ 0, 1, 13, 14] ,
                   1: [ 1, 2, 14, 15] ,
                   2: [ 2, 3, 15, 16] ,
                   3: [ 3, 4, 16, 17] ,
                   4: [ 4, 5, 17, 18] ,
                   5: [ 5, 6, 18, 19] ,
                   6: [ 6, 7, 19, 20] ,
                   7: [ 7, 8, 20, 21] ,
                   8: [ 8, 9, 21,  22] ,
                   9: [ 9, 10, 22, 23] ,
                   10: [10, 11, 23, 24] ,
                   11: [11, 12, 24, 25] ,
                   13: [13, 14, 26, 27] ,
                   14: [14, 15, 27, 28] ,
                   15: [15, 16, 28, 29] ,
                   16: [16, 17, 29, 30] ,
                   17: [17, 18, 30, 31] ,
                   18: [18, 19, 31, 32] ,
                   19: [19, 20, 32, 33] ,
                   20: [20, 21, 33, 34] ,
                   21: [21, 22, 34, 35] ,
                   22: [22, 23, 35, 36] ,
                   23: [23, 24, 36, 37] ,
                   24: [24, 25, 37, 38] ,
                   26: [26, 27, 39, 40] ,
                   27: [27, 28, 40, 41] ,
                   28: [28, 29, 41, 41] ,
                   29: [29, 30, 41, 42] ,
                   30: [30, 31, 42, 43] ,
                   31: [31, 32, 43, 44] ,
                   32: [32, 33, 44, 44] ,
                   33: [33, 34, 44, 45] ,
                   34: [34, 35, 45, 46] ,
                   35: [35, 36, 46, 47] ,
                   36: [36, 37, 47, 47] ,
                   37: [37, 38, 47, 48] ,
                   39: [39, 40, 49, 50] ,
                   40: [40, 41, 50, 50] ,
                   41: [41, 42, 50, 51] ,
                   42: [42, 43, 51, 52] ,
                   43: [43, 44, 52, 52] ,
                   44: [44, 45, 52, 53] ,
                   45: [45, 46, 53, 54] ,
                   46: [46, 47, 54, 54] ,
                   47: [47, 48, 54, 55] ,
                   49: [49, 50, 56, 57] ,
                   50: [50, 51, 57, 57] ,
                   51: [51, 52, 57, 58] ,
                   52: [52, 53, 58, 59] ,
                   53: [53, 54, 59, 59] ,
                   54: [54, 55, 59, 60] ,
                   56: [56, 57, 61, 61] ,
                   57: [57, 58, 61, 61] ,
                   58: [58, 59, 61, 61] ,
                   59: [59, 60, 61, 61]}

class data_structure():

    def __init__(self):
        self.set_symmetry()
        self.set_R()
        self.set_phi()
        self.set_theta()
        self.set_nConf()
        self.degree2radius()
        self.set_num_of_atoms()

    def set_theta(self):
        self.THETA_angles = {0: [float(i) for i in range(0, 359, 15)],
                             1: [float(i) for i in range(0, 359, 15)],
                             2: [float(i) for i in range(0, 359, 15)],
                             3: [float(i) for i in range(0, 359, 20)],
                             4: [float(i) for i in range(0, 359, 30)],
                             5: [float(i) for i in range(0, 359, 45)],
                             6: [0.0]}

    def set_symmetry(self):
        self.symface = ['xy']

    def set_num_of_atoms(self):
        self.n1 = 12
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: O, CA, next CA
        """
        ratio_CA1 = 0.5618541311379575
        ratio_CA2 = 0.4381458688620425
        return -np.array(a1)*ratio_CA1 - np.array(a2)*ratio_CA2

    def calt_vec1(self, a0, a1, a2): return a0, [1,0,0]

    def calt_vec2(self, a0, a1, a2): return a1, [0,1,0]

    def set_R(self):
        self.R_NDX = [2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,4.0,4.2,4.5,4.7,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0]
        self.nDist = len(self.R_NDX)
        self.DR = []
        for i in range(len(self.R_NDX)-1):
            self.DR.append(self.R_NDX[i+1]-self.R_NDX[i])
        self.DR.append( 0.0 )

    def set_phi(self):
        """
        default: 0~90.0, d_ang = 15.0
        """
        self.PHI_angles = [i*15.0 for i in range(0,7)]
        self.nPhi = len(self.PHI_angles)

    def set_nConf(self): # depending on the distance
        self.nConf = []
        self.nNorm = []
        for i in range(self.nDist):
            if self.R_NDX[i] < 2.51:
                self.nConf.append(54)
                self.nNorm.append(4)
            elif self.R_NDX[i] > 2.51 and self.R_NDX[i] < 3.51:
                self.nConf.append(210)
                self.nNorm.append(4)
            elif self.R_NDX[i] > 3.51 and self.R_NDX[i] < 5.51:
                self.nConf.append(54)
                self.nNorm.append(2)
            else:
                self.nConf.append(26)
                self.nNorm.append(2)


    def degree2radius(self):
        for i in range(self.nPhi): self.PHI_angles[i] *= D2R
        self.NTheta = {}
        self.nGrid = 0
        for i in range(self.nPhi):
            self.NTheta[i] = len(self.THETA_angles[i])
            for j in range(self.NTheta[i]):
                self.THETA_angles[i][j] *= D2R
                self.nGrid += self.NTheta[i]


class wtr_structure(data_structure):

    def set_symmetry(self):
        self.symface = ['xy','xz']

    def set_num_of_atoms(self):
        self.n1 = 3
        self.n2 = 3

    def set_theta(self):
        self.THETA_angles = {0: [float(i) for i in range(0, 181, 15)],
                             1: [float(i) for i in range(0, 181, 15)],
                             2: [float(i) for i in range(0, 181, 15)],
                             3: [float(i) for i in range(0, 181, 20)],
                             4: [float(i) for i in range(0, 181, 30)],
                             5: [float(i) for i in range(0, 181, 45)],
                             6: [0.0]}

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: 'O','H','H' in H2O.
        """
        totalM = Mass_table['O'] + Mass_table['H'] + Mass_table['H']
        com = [a0[k]*Mass_table['O']+a1[k]*Mass_table['H']+a2[k]*Mass_table['H'] for k in range(3)]
        com = [com[k]/totalM for k in range(3)]
        dvec = np.array(com)
        return -dvec

    def calt_vec1(self, a0, a1, a2):
        return 0.5*(np.array(a1)+np.array(a2)), [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        vec = get_normal(np.array(a2), np.array(a1))
        return vec, [0,0,1]


Grid_Quarts = { 'wtr': Grid_Quarts_wtr}
## Data Structures of the fragments:
DS = {  'wtr': wtr_structure()}
