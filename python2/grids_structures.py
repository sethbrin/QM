import numpy as np
from trans_rot_coords import get_normal, get_unit

__version__ = 0.3

D2R = 3.14159265358979/180.0

## Indices of one square:
Grid_Quarts_ben = { 0: [ 0, 1, 3, 4] ,
         1: [ 1, 2, 4, 5] ,
         3: [ 3, 4, 6, 7] ,
         4: [ 4, 5, 7, 8] ,
         6: [ 6, 7, 9, 10] ,
         7: [ 7, 8, 10, 11] ,
         9: [ 9, 10, 12, 13] ,
        10: [10, 11, 13, 14] ,
        12: [12, 13, 15, 15] ,
        13: [13, 14, 16, 16] ,
        15: [15, 16, 17, 17] }
Grid_Quarts_arg = { 0: [0, 1, 5, 6],
         1: [ 1,  2,  6, 7],
         2: [ 2,  3,  7, 8],
         3: [ 3,  4,  8, 9],
         5: [ 5,  6, 10, 11],
         6: [ 6,  7, 11, 12],
         7: [ 7,  8, 12, 13],
         8: [ 8,  9, 13, 14],
        10: [10, 11, 15, 16],
        11: [11, 12, 16, 17],
        12: [12, 13, 17, 18],
        13: [13, 14, 18, 19],
        15: [15, 16, 20, 21],
        16: [16, 17, 21, 21],
        17: [17, 18, 22, 22],
        18: [18, 19, 22, 23],
        20: [20, 21, 24, 25],
        21: [21, 22, 25, 25],
        22: [22, 23, 25, 26],
        24: [24, 25, 27, 27],
        25: [25, 26, 27, 27] }
Grid_Quarts_nh4 = {0:[0,  0,  1,  3],
         1: [ 1,  2,  4,  5],
         2: [ 2,  3,  5,  7],
         4: [ 4,  5,  8,  9],
         5: [ 5,  6,  9, 11],
         6: [ 6,  7, 11, 12],
         8: [ 8,  9, 13, 14],
         9: [ 9, 10, 14, 15],
        10: [10, 11, 15, 16],
        11: [11, 12, 16, 17],
        13: [13, 14, 18, 19],
        14: [14, 15, 19, 20],
        15: [15, 16, 20, 21],
        16: [16, 17, 21, 22],
        18: [18, 19, 23, 24],
        19: [19, 20, 24, 25],
        20: [20, 21, 25, 26],
        21: [21, 22, 26, 27],
        23: [23, 24, 28, 29],
        24: [24, 25, 29, 30],
        25: [25, 26, 30, 31],
        26: [26, 27, 31, 32],
        28: [28, 29, 33, 34],
        29: [29, 30, 34, 35],
        30: [30, 31, 35, 36],
        31: [31, 32, 36, 37],
        33: [33, 34, 38, 39],
        34: [34, 35, 39, 40],
        35: [35, 36, 40, 41],
        36: [36, 37, 41, 42],
        38: [38, 39, 43, 44],
        39: [39, 40, 44, 44],
        40: [40, 41, 44, 45],
        41: [41, 42, 45, 46],
        43: [43, 44, 47, 48],
        44: [44, 45, 48, 48],
        45: [45, 46, 48, 49],
        47: [47, 48, 50, 50],
        48: [48, 49, 50, 50]}
Grid_Quarts_prp = { 0: [ 0, 1, 13, 14] ,
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
        52: [52, 53, 58, 58] ,
        53: [53, 54, 58, 59] ,
        54: [54, 55, 59, 59] ,
        56: [56, 57, 60, 60] ,
        57: [57, 58, 60, 60] ,
        58: [58, 59, 60, 60] }
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

Grid_Quarts_bck = { 0: [ 0, 1, 24, 25] ,
         1: [ 1, 2, 25, 26] ,
         2: [ 2, 3, 26, 27] ,
         3: [ 3, 4, 27, 28] ,
         4: [ 4, 5, 28, 29] ,
         5: [ 5, 6, 29, 30] ,
         6: [ 6, 7, 30, 31] ,
         7: [ 7, 8, 31, 32] ,
         8: [ 8, 9, 32, 33] ,
         9: [ 9, 10, 33, 34] ,
        10: [10, 11, 34, 35] ,
        11: [11, 12, 35, 36] ,
        12: [12, 13, 36, 37] ,
        13: [13, 14, 37, 38] ,
        14: [14, 15, 38, 39] ,
        15: [15, 16, 39, 40] ,
        16: [16, 17, 40, 41] ,
        17: [17, 18, 41, 42] ,
        18: [18, 19, 42, 43] ,
        19: [19, 20, 43, 44] ,
        20: [20, 21, 44, 45] ,
        21: [21, 22, 45, 46] ,
        22: [22, 23, 46, 47] ,
        23: [23,  0, 47, 24] ,
        24: [24, 25, 48, 49] ,
        25: [25, 26, 49, 50] ,
        26: [26, 27, 50, 51] ,
        27: [27, 28, 51, 52] ,
        28: [28, 29, 52, 53] ,
        29: [29, 30, 53, 54] ,
        30: [30, 31, 54, 55] ,
        31: [31, 32, 55, 56] ,
        32: [32, 33, 56, 57] ,
        33: [33, 34, 57, 58] ,
        34: [34, 35, 58, 59] ,
        35: [35, 36, 59, 60] ,
        36: [36, 37, 60, 61] ,
        37: [37, 38, 61, 62] ,
        38: [38, 39, 62, 63] ,
        39: [39, 40, 63, 64] ,
        40: [40, 41, 64, 65] ,
        41: [41, 42, 65, 66] ,
        42: [42, 43, 66, 67] ,
        43: [43, 44, 67, 68] ,
        44: [44, 45, 68, 69] ,
        45: [45, 46, 69, 70] ,
        46: [46, 47, 70, 71] ,
        47: [47, 24, 71, 48] ,
        48: [48, 49, 72, 73] ,
        49: [49, 50, 73, 73] ,
        50: [50, 51, 73, 74] ,
        51: [51, 52, 74, 75] ,
        52: [52, 53, 75, 76] ,
        53: [53, 54, 76, 76] ,
        54: [54, 55, 76, 77] ,
        55: [55, 56, 77, 78] ,
        56: [56, 57, 78, 79] ,
        57: [57, 58, 79, 79] ,
        58: [58, 59, 79, 80] ,
        59: [59, 60, 80, 81] ,
        60: [60, 61, 81, 82] ,
        61: [61, 62, 82, 82] ,
        62: [62, 63, 82, 83] ,
        63: [63, 64, 83, 84] ,
        64: [64, 65, 84, 85] ,
        65: [65, 66, 85, 86] ,
        66: [66, 67, 86, 86] ,
        67: [67, 68, 86, 87] ,
        68: [68, 69, 87, 88] ,
        69: [69, 70, 88, 88] ,
        70: [70, 71, 88, 89] ,
        71: [71, 48, 89, 72] ,
        72: [72, 73, 90, 91] ,
        73: [73, 74, 91, 91] ,
        74: [74, 75, 91, 92] ,
        75: [75, 76, 92, 93] ,
        76: [76, 77, 93, 93] ,
        77: [77, 78, 93, 94] ,
        78: [78, 79, 94, 95] ,
        79: [79, 80, 95, 95] ,
        80: [80, 81, 95, 96] ,
        81: [81, 82, 96, 97] ,
        82: [82, 83, 97, 97] ,
        83: [83, 84, 97, 98] ,
        84: [84, 85, 98, 99] ,
        85: [85, 86, 99, 99] ,
        86: [86, 87, 99,100] ,
        87: [87, 88,100,101] ,
        88: [88, 89,101,101] ,
        89: [89, 72,101, 90] ,
        90: [90, 91,102,103] ,
        91: [91, 92,103,103] ,
        92: [92, 93,103,104] ,
        93: [93, 94,104,105] ,
        94: [94, 95,105,105] ,
        95: [95, 96,105,106] ,
        96: [96, 97,106,107] ,
        97: [97, 98,107,107] ,
        98: [98, 99,107,108] ,
        99: [99, 100,108,109] ,
        100:[100,101,109,109],
        101:[101, 90,109,102],
        102:[102,103,110,110],
        103:[103,104,110,110],
        104:[104,105,110,110],
        105:[105,106,110,110],
        106:[106,107,110,110],
        107:[107,108,110,110],
        108:[108,109,110,110],
        109:[109,102,110,110]}

class data_structure():
    #param = {'CT-OW':[0.1590, 3.8210],
    #        'CT-HW':[0.2499, 2.2993],
    #        'HC-OW':[0.2972, 2.4279],
    #        'HC-HW':[0.0074, 2.9165],
    #        'qCT': -0.0951,
    #        'qHC':  0.0317,
    #        'qOW': -0.6888,
    #        'qHW':  0.3444}
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0000, 3.8279],
    'HC-HW':[0.0000, 1.4165],
    'qCT':  0.0000,
    'qHC':  0.0000,
    'qOW':  0.0000,
    'qHW':  0.0000}

    def __init__(self, fftype='14-7'):
        self.set_symmetry()
        self.set_R()
        self.set_phi()
        self.set_theta()
        self.set_nConf()
        self.degree2radius()
        self.set_H_correction(False)
        self.set_num_of_atoms()

        if fftype=='none':
            self.param['qCT'] = 0.0
            self.param['qHC'] = 0.0
            self.param['qOW'] = 0.0
            self.param['qHW'] = 0.0
        temp = {'none':self.nocorr,
            '9-6': self.lj_9_6,
            '12-6': self.lj_12_6,
            '14-7': self.lj_b_14_7}
        self.lj = temp[fftype]

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
        self.R_NDX = [2.0,2.2,2.5,2.7,3.0,3.2,3.5,3.7,4.0,4.5,4.7,5.0,5.5,6.0,6.5,7.0,8.0]
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

    def set_nConf(self):
        self.nConf = 26
        self.nNorm =  2

    def degree2radius(self):
        for i in range(self.nPhi): self.PHI_angles[i] *= D2R
        self.NTheta = {}
        self.nGrid = 0
        for i in range(self.nPhi):
            self.NTheta[i] = len(self.THETA_angles[i])
            for j in range(self.NTheta[i]):
                self.THETA_angles[i][j] *= D2R
            self.nGrid += self.NTheta[i]

    def set_H_correction(self, IsCorr):
        self.IsHmm = IsCorr
        self.Hmm = []

    def lj_12_6(self, r, eps, sig):
        tempr = (sig/r)**6
        return 4.*eps*(tempr*tempr-tempr)

    def lj_9_6(self, r, eps, sig):
        tempr = (sig/r)**3
        tempr2 = tempr*tempr
        return eps*(2.*tempr*tempr2 - 3.*tempr2)

    def lj_b_14_7(self, r, eps, sig):
        rou = r/sig
        res = eps*(1.07/(rou+0.07))**7 * (1.12/(rou**7+0.12)-2.0)
        return res

    def coulomb(self, q0,q1,r):
            return q0*q1/r*332.5

    def nocorr(self,r,eps,sig): return 0.0

class bck_structure(data_structure):
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0100, 2.9279],
    'HC-HW':[0.0100, 2.2165],
    'qCT':  0.0000,
    'qHC':  0.0000,
    'qOW': -0.6800,
    'qHW':  0.3400}
    #param = {'CT-OW':[0.0000, 3.8210],
    #    'CT-HW':[0.0000, 2.2993],
    #    'HC-OW':[0.0100, 3.8279],
    #    'HC-HW':[0.0000, 2.2165],
    #    'qCT':  0.0000,
    #    'qHC':  0.0000,
    #    'qOW': -0.6800,
    #    'qHW':  0.3400}

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

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [5,6,7,9,10,11]

class prp_structure(data_structure):
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0100, 2.5279],
    'HC-HW':[0.0100, 2.2165],
    'qCT':  0.0000,
    'qHC':  0.0000,
    'qOW': -0.6800,
    'qHW':  0.3400}

    def set_theta(self):
        self.THETA_angles = {0: [float(i) for i in range(0, 181, 15)],
                          1: [float(i) for i in range(0, 181, 15)],
                          2: [float(i) for i in range(0, 181, 15)],
                          3: [float(i) for i in range(0, 181, 20)],
                          4: [float(i) for i in range(0, 181, 30)],
                          5: [float(i) for i in range(0, 181, 60)],
                          6: [0.0]}
    def set_symmetry(self):
        self.symface = ['xy','xz']

    def set_num_of_atoms(self):
        self.n1 = 11
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C1, C2, C3 in C1-C2-C3
        """
        dvec = np.array(a0) + np.array(a2)
        dvec = 0.5*dvec + np.array(a1)
        dvec *= 0.5
        return -dvec

    def calt_vec1(self, a0, a1, a2): return a1, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a0) - np.array(a1)
        v2 = np.array(a2) - np.array(a1)
        vec = get_normal(v1,v2)
        return vec, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [1,2,3,8,9,10]


class ben_structure(data_structure):
    def set_theta(self):
        self.THETA_angles = {0: [float(i) for i in range(0, 31, 15)], 
                          1: [float(i) for i in range(0, 31, 15)],
                          2: [float(i) for i in range(0, 31, 15)],
                          3: [float(i) for i in range(0, 31, 15)],
                          4: [float(i) for i in range(0, 31, 15)],
                          5: [float(i) for i in range(0, 31, 30)],
                          6: [0.0]}
    def set_symmetry(self):
        self.symface = ['xy','zben']

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: Cx(0), C2, and C4 in C0(x)-C1-C2-C3-C4-C5
            or CG, CE1 and CE2 in PHE
        """
        dvec = np.array(a1) + np.array(a2) + np.array(a0)
        dvec /= 3.0
        return -dvec

    def calt_vec1(self, a0, a1, a2): return a0, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a1) - np.array(a0)
        v2 = np.array(a2) - np.array(a0)
        vec = get_normal(v1,v2)
        return vec, [0,0,1]

    def set_num_of_atoms(self):
        self.n1 = 12
        self.n2 = 3

class arg_structure(data_structure):
    def set_theta(self):
        self.THETA_angles = { 0: [float(i) for i in range(0, 61, 15)],
                      1: [float(i) for i in range(0, 61, 15)],
                      2: [float(i) for i in range(0, 61, 15)],
                      3: [float(i) for i in range(0, 61, 15)],
                      4: [float(i) for i in range(0, 61, 20)],
                      5: [float(i) for i in range(0, 61, 30)],
                      6: [0.0]}
    def set_symmetry(self):
        self.symface = ['xy', 'zarg']

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C0, N1(x), N3(xy); [CZ, NH1, NH2] in pdb files.
        """
        return -np.array(a0)

    def calt_vec1(self, a0, a1, a2): return a1, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a1)
        v2 = np.array(a2)
        vec = get_normal(v1,v2)
        return vec, [0,0,1]

    def set_num_of_atoms(self):
        self.n1 = 10
        self.n2 = 3

class nh4_structure(data_structure):
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0100, 3.1279],
    'HC-HW':[0.0100, 2.4165],
    'qCT':  0.0000,
    'qHC':  0.2900,
    'qOW': -0.6800,
    'qHW':  0.3400}

    def set_phi(self):
        """
        -90.0~90.0, d_ang = 15.0
        """
        self.PHI_angles = [i*15.0 for i in range(-6,7)]
        self.nPhi = len(self.PHI_angles)

    def set_theta(self):
        self.THETA_angles = { 0: [0.0],
                      1: [float(i) for i in range(0, 61, 30)],
                      2: [float(i) for i in range(0, 61, 20)],
                      3: [float(i) for i in range(0, 61, 15)],
                      4: [float(i) for i in range(0, 61, 15)],
                      5: [float(i) for i in range(0, 61, 15)],
                      6: [float(i) for i in range(0, 61, 15)],
                      7: [float(i) for i in range(0, 61, 15)],
                      8: [float(i) for i in range(0, 61, 15)],
                      9: [float(i) for i in range(0, 61, 15)],
                     10: [float(i) for i in range(0, 61, 20)],
                     11: [float(i) for i in range(0, 61, 30)],
                     12: [0.0]}

    def set_symmetry(self):
        self.symface = ['zarg']

    def calt_dvec(self, a0,a1,a2):
        """
        a0, a1, a2: CH3, NH3, H in CH3-NH2-H;
        [CE, NZ, HZ1] in 'LYS'.
        """
        return -0.5*(np.array(a0)+np.array(a1))

    def calt_vec1(self, a0, a1, a2):
        return a1, [0,0,1]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a2)-np.array(a1)
        return np.array([v1[0],v1[1],0.0]), [1,0,0]

    def set_num_of_atoms(self):
        self.n1 = 8
        self.n2 = 3

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [5,6,7]

class alc_structure(bck_structure):
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0100, 3.0279],
    'HC-HW':[0.0100, 2.2165],
    'qCT':  0.0000,
    'qHC':  0.0200,
    'qOW': -0.6800,
    'qHW':  0.3400}

    def set_num_of_atoms(self):
        self.n1 = 6
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        NOTE: negative vector of the center!
        a0, a1, a2: CA, O, H
        """
        return -(np.array(a0)+np.array(a1))*0.5

    def calt_vec1(self, a0, a1, a2):
        return a1, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a2)-np.array(a1)
        v2 = np.array(a0)-np.array(a1)
        vn = np.cross(v1, v2)
        return vn, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [1,2,3]

class cys_structure(alc_structure):
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0100, 3.0279],
    'HC-HW':[0.0100, 1.5165],
    'qCT':  0.0000,
    'qHC':  0.0800,
    'qOW': -0.6800,
    'qHW':  0.3400}

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [1,2,3]

class co2_structure(prp_structure):
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0100, 3.0279],
    'HC-HW':[0.0100, 1.0165],
    'qCT':  0.0000,
    'qHC':  0.0900,
    'qOW': -0.6800,
    'qHW':  0.3400}

    def set_num_of_atoms(self):
        self.n1 = 7
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C4, O5, O6; [CD, OE1, OE2] in 'GLU'.
        """
        dvec = np.array(a0)
        return -dvec

    def calt_vec1(self, a0, a1, a2):
        return 0.5*(np.array(a1)+np.array(a2)), [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        vec = get_normal(np.array(a2), np.array(a1))
        return vec, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [1,2,3]

class hip_structure(prp_structure):
    def set_num_of_atoms(self):
        self.n1 = 10
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C0, N1, and N4
        """
        dvec = np.array(a1) + np.array(a2)
        dvec = 0.5*dvec
        return -dvec

    def calt_vec1(self, a0, a1, a2): return a0, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a1) - np.array(a0)
        v2 = np.array(a2) - np.array(a0)
        vec = get_normal(v1,v2)
        return vec, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm = False
        self.Hmm = []

class met_structure(prp_structure):
    def set_num_of_atoms(self):
        self.n1 = 9
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C0, S4, C5 in C0-S4-C5; [CG, SD, CE] in 'MET'.
        """
        dvec = np.array(a0) + np.array(a2)
        dvec = 0.5*dvec + np.array(a1)
        dvec *= 0.5
        return -dvec

    def calt_vec1(self, a0, a1, a2): return a1, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a0) - np.array(a1)
        v2 = np.array(a2) - np.array(a1)
        vec = get_normal(v1,v2)
        return vec, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [1,2,3,6,7,8]

class hid_structure(bck_structure):
    def set_num_of_atoms(self):
        self.n1 = 9
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C0, N1, N4; [CE1, NE2(with H), ND1(without H)] in 'HIS'.
        """
        dvec = np.array(a1) + np.array(a2)
        dvec = 0.5*dvec
        return -dvec

    def calt_vec1(self, a0, a1, a2): return a0, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a1) - np.array(a0)
        v2 = np.array(a2) - np.array(a0)
        vec = get_normal(v1,v2)
        return vec, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm = False
        self.Hmm = []

class trp_structure(bck_structure):
    def set_num_of_atoms(self):
        self.n1 = 16
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C3, C5, C2; [CD2, CE2, CD1] in 'TRP'.
        """
        dvec = 0.5*(np.array(a0)+np.array(a1))
        return -dvec

    def calt_vec1(self, a0, a1, a2):
        return a1, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        return a2, [0,1,0]

    def set_H_correction(self, IsCorr):
        self.IsHmm = False
        self.Hmm = []

class gln_structure(bck_structure):
    param = {'CT-OW':[0.0000, 3.8210],
    'CT-HW':[0.0000, 2.2993],
    'HC-OW':[0.0100, 3.3279],
    'HC-HW':[0.0100, 1.6165],
    'qCT':  0.0000,
    'qHC':  0.1300,
    'qOW': -0.6800,
    'qHW':  0.3400}

    def set_num_of_atoms(self):
        self.n1 = 9
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: CH3,C, O; [CB,CG,OD1] in 'ASN'.
        """
        dvec = np.array(a1)
        return -dvec

    def calt_vec1(self, a0, a1, a2): return a2, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        v1 = np.array(a0)
        v2 = np.array(a2)
        vec = get_normal(v1,v2)
        return vec, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [4,5,6]

class tyr_structure(bck_structure):
    #def set_R(self):
    #    self.R_NDX = [2.0, 2.2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0]
    #    self.nDist = len(self.R_NDX)
    #    self.DR = []
    #    for i in range(len(self.R_NDX)-1):
    #        self.DR.append(self.R_NDX[i+1]-self.R_NDX[i])
    #    self.DR.append( 0.0 )

    def set_num_of_atoms(self):
        self.n1 = 13
        self.n2 = 3

    def calt_dvec(self, a0, a1, a2):
        """
        a0, a1, a2: C1,C5,C6; [CE1(torsion is 0.0 to HH), CE2, CZ] in 'PHE'.
        """
        dvec = (np.array(a0)+np.array(a1))*0.5
        return -dvec

    def calt_vec1(self, a0, a1, a2): return a2, [1,0,0]

    def calt_vec2(self, a0, a1, a2):
        return a0, [0,1,0]

    def set_H_correction(self, IsCorr):
        self.IsHmm = False
        self.Hmm = []


class wtr_structure(prp_structure):
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
        dvec = np.array(a0)
        return -dvec

    def calt_vec1(self, a0, a1, a2):
        #v1 = get_unit(a1)
        #v2 = get_unit(a2) 
        #return 0.5*(v1+v2), [1,0,0]
        return 0.5*(np.array(a1)+np.array(a2)), [1,0,0]

    def calt_vec2(self, a0, a1, a2):

        vec = get_normal(np.array(a2), np.array(a1))
        return vec, [0,0,1]

    def set_H_correction(self, IsCorr):
        self.IsHmm =  False
        self.Hmm = [1,2,3]


Grid_Quarts = { 'bck': Grid_Quarts_bck,
        'prp': Grid_Quarts_prp,
        'ben': Grid_Quarts_ben,
        'alc': Grid_Quarts_bck,
        'co2': Grid_Quarts_prp,
        'cys': Grid_Quarts_bck,
        'arg': Grid_Quarts_arg,
        'hip': Grid_Quarts_prp,
        'hid': Grid_Quarts_bck,
        'met': Grid_Quarts_prp,
        'gln': Grid_Quarts_bck,
        'nh4': Grid_Quarts_nh4,
        'trp': Grid_Quarts_bck,
        'tyr': Grid_Quarts_bck,
        'wtr': Grid_Quarts_wtr}
## Data Structures of the fragments:
#DS = {  'bck': bck_structure(),
#    'ben': ben_structure(),
#    'prp': prp_structure(),
#    'alc': alc_structure(),
#    'co2': co2_structure(),
#    'cys': cys_structure(),
#    'arg': arg_structure(),
#    'hip': hip_structure(),
#    'hid': hid_structure(),
#    'met': met_structure(),
#    'gln': gln_structure(),
#    'nh4': nh4_structure(),
#    'trp': trp_structure(),
#    'tyr': tyr_structure(),
#    'wtr': wtr_structure()}
DS = {'wtr': wtr_structure()}
