#!/usr/bin/env python

import numpy as np
from calculate_energy_new_coords import coordinates
from read_energy_force_new import energy_force_data

## Get the energy for two water molecules with coordinates coors0 and coors1.
## coors0 and coors1 are lists of floats or numpy arrays
def interp_coors(coors0, coors1, allconfig):
    interp = coordinates(3,3, 'wtr', None) #create a coordinates object
    ## Add center molecule atoms:
    interp.addatomfloat('O', coors0[0])
    interp.addatomfloat('H', coors0[1])
    interp.addatomfloat('H', coors0[2])
    ## Add probe molecule atoms:
    interp.addatomfloat('O', coors1[0])
    interp.addatomfloat('H', coors1[1])
    interp.addatomfloat('H', coors1[2])
    ## Reorient and interpolate:
    interp.ReorientToOrigin() # reorient
    interp.MirrorAll() # mirror
    interp.indexing() #???
    interp.calt_conf_energy(allconfig) # do the interpolation E = interp.get_interp_energy()
    E = interp.get_interp_energy()
    return E

## Evaluate the interaction between two water molecules, the coordinates of which
## are stored in allcoors as a list of 18 floats in the order of [XO0, YO0, ZO0, XH0, 
## YH0, ZH0, XH0, YH0, ZH0, XO1, YO1, ZO1, XH1, YH1, ZH1, XH1, YH1, ZH1].
## Return energy, force and torque acting on molecule0 as a 7 components list
def calc_energy_force(allcoors, allconfig):
    mass = np.array([15.99900, 1.00800, 1.00800])
    coors0, coors1 = np.array(allcoors).reshape((2,3,3))
    E0 = interp_coors(coors0, coors1, allconfig)
    dr = 0.00001
    ## Translate along X
    newcoors0 = coors0.copy()
    newcoors0[:,0] += dr
    E1 = interp_coors(newcoors0, coors1, allconfig)
    fx = (E0 - E1) / dr
    ## Translate along Y
    newcoors0 = coors0.copy()
    newcoors0[:,1] += dr
    E1 = interp_coors(newcoors0, coors1, allconfig)
    fy = (E0 - E1) / dr
    ## Translate along Z
    newcoors0 = coors0.copy()
    newcoors0[:,2] += dr
    E1 = interp_coors(newcoors0, coors1, allconfig) 
    fz = (E0 - E1) / dr
    dtheta = 0.00001
    com0 = np.dot(mass, coors0) / mass.sum()
    dcoors0 = coors0 - com0
    c, s = np.cos(dtheta), np.sin(dtheta)
    ## Rotate molecule 0 about X passing COM
    R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    newcoors0 = np.dot(R, dcoors0) + com0
    E1 = interp_coors(newcoors0, coors1, allconfig)
    tx = (E0 - E1) / dtheta
    ## Rotate molecule 0 about Y passing COM
    R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    newcoors0 = np.dot(R, dcoors0) + com0
    E1 = interp_coors(newcoors0, coors1, allconfig)
    ty = (E0 - E1) / dtheta
    ## Rotate molecule 0 about Z passing COM
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]) 
    newcoors0 = np.dot(R, dcoors0) + com0
    E1 = interp_coors(newcoors0, coors1, allconfig)
    tz = (E0 - E1) / dtheta
    return E0, fx, fy, fz, tx, ty, tz

class QMInterpolation():
  def __init__(self, ftype):
    self.ftype = ftype
    import os
    cur_dir = os.path.split(os.path.realpath(__file__))[0]
    database_name = 'Dimer_deltaE_data_mp2_wtr_wtr.txt.gz'
    database_name = cur_dir + os.sep + database_name
    self.allconf = energy_force_data(database_name, ftype)

  #
  # ('O', [8.0,0.10293690, 0.11093940, 0.07435482],
  #  'H', [8.0,0.10293690, 0.11093940, 0.07435482],
  #  'H', [8.0,0.10293690, 0.11093940, 0.07435482])
  def calculate(self, lhs, rhs):
    allcoors = []
    for item in lhs:
      allcoors.extend(item[1][1:])
    for item in rhs:
      allcoors.extend(item[1][1:])
    print "allcoors:", allcoors
    tempe = calc_energy_force(allcoors, self.allconf)
    print "res:", tempe
    return tempe[1:]

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
    allconf = energy_force_data('Dimer_deltaE_data_mp2_wtr_wtr.txt.gz', 'wtr')
    allcoors = [0.10293690, 0.11093940, 0.07435482, 0.73097151, 0.83572755, 0.02040124, -0.22532655, 0.13689786, 0.97669930,
        2.85534030, 0.60573860, -1.03993027, 1.93176673, 0.51127346, -0.79346626, 3.16776562, -0.29302760, -1.17133059]
    allcoors = [-5.86, -8.331, -8.296, -5.323191, -8.138266, -9.117395, -5.247876, -8.606731,
        -7.554867, -8.151347,
        -8.087137, -2.63456, -8.151347, -8.087137, -2.63456, -8.151347, -8.087137, -2.63456]
    print calc_energy_force(allcoors, allconf)

