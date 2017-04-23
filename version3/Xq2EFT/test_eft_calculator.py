#!/usr/bin/env python2

import numpy as np
from time import time

from eft_calculator import EFT_calculator
import tools


def load_coordinates(name):
    lines = open('random/'+name).readlines()[-7:-1]
    coors = [[float(item) for item in line.split()[2:5]] for line in lines]
    return np.array(coors)

def test_random_set():
    ener = []
    force = []
    torque = []
    t0 = time()
    Ofile = open('wtr01.txt', 'w')
    fout1 = open('wtr_wtr_force.txt','w')
    fout2 = open('wtr_wtr_torque.txt','w')
    for i in range(2, 2000):
    #for i in range(2, 3):
        name = 'test%04d.inp' % i
        coors = load_coordinates(name)
        eft = calculator.eval(coors[:3], coors[3:])
        Ofile.write('%s %12.7f\n'%(name, eft[0]))
        Ofile.flush()
        fout1.write('%s %12.7f %12.7f %12.7f\n'%(name, eft[1], eft[2], eft[3]))
        fout1.flush()
        fout2.write('%s %12.7f %12.7f %12.7f\n'%(name, eft[4], eft[5], eft[6]))
        fout2.flush()
        ener.append(eft[0])
        force.append(eft[1:4])
        torque.append(eft[4:7])
    t1 = time()
    print 'took %.1f s to evaluate the random set' % (t1 - t0)
    Ofile.close()
    fout1.close()
    fout2.close()
    return ener, force, torque


if __name__ == '__main__':
    order = 2
    calculator = EFT_calculator(order)
    t0 = time()
    calculator.setup('grid.dat')
    t1 = time()
    print 'took %.1f s to fill the grid' % (t1 - t0)
    ener, force, torque = test_random_set()

