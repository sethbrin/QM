import os, os.path, sys
import numpy as np

D2R = 3.14159265358979/180.0

class orient_data_structure():
    def __init__(self):
        self.set_phi()
        self.set_theta()
        self.degree2radius()

    def set_theta(self):
        self.THETA_angles = { 0: [float(i) for i in range(0, 359, 30)],
                             1: [float(i) for i in range(0, 359, 30)],
                             2: [float(i) for i in range(0, 359, 45)],
                             3: [0.0]}

    def set_phi(self):
        self.PHI_angles = [-i*30.0 for i in range(0,4)]
        self.nPhi = len(self.PHI_angles)

    def degree2radius(self):
        for i in range(self.nPhi): self.PHI_angles[i] *= D2R
        self.NTheta = {}
        self.nGrid = 0
        for i in range(self.nPhi):
            self.NTheta[i] = len(self.THETA_angles[i])
            for j in range(self.NTheta[i]):
                self.THETA_angles[i][j] *= D2R
            self.nGrid += self.NTheta[i]



class orient_wtr_structure(orient_data_structure):
    def set_theta(self):
        self.THETA_angles = { 0: [0.0],
                             1: [float(i) for i in range(0, 359, 45)],
                             2: [float(i) for i in range(0, 359, 30)],
                             3: [float(i) for i in range(0, 359, 15)],
                             4: [float(i) for i in range(0, 359, 15)],
                             5: [float(i) for i in range(0, 359, 15)],
                             6: [float(i) for i in range(0, 359, 15)],
                             7: [float(i) for i in range(0, 359, 15)],
                             8: [float(i) for i in range(0, 359, 15)],
                             9: [float(i) for i in range(0, 359, 15)],
                             10: [float(i) for i in range(0, 359, 30)],
                             11: [float(i) for i in range(0, 359, 45)],
                             12: [0.0]}

    def set_phi(self):
        self.PHI_angles = [i*15.0 for i in range(-6,7)]
        self.nPhi = len(self.PHI_angles)


OrientDS = {  'wtr': orient_wtr_structure()}
