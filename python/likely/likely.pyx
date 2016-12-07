# distutils: language = c++
# distutils: sources = TriCubicInterpolator.cpp

cimport numpy as np
import numpy as np

cdef extern from "TriCubicInterpolator.h" namespace "likely": 
    cdef cppclass TriCubicInterpolator:
        TriCubicInterpolator(double* data, double spacing, int n1, int n2, int n3) except +
        double test(double x, double y, double z)

cdef class Grid3DInterpolator:
    cdef TriCubicInterpolator* thisptr
    cdef double stepx, stepy, stepz
    def __cinit__(self, np.ndarray[np.float64_t, ndim=3] data, double stepx=1.0, double stepy=1.0, double stepz=1.0):
        self.thisptr = new TriCubicInterpolator(<double*>(data.data), 1, data.shape[0], data.shape[1], data.shape[2])
        self.stepx = stepx
        self.stepy = stepy
        self.stepz = stepz
    
    def __dealloc__(self):
        del self.thisptr
    
    def test(self, double x, double y, double z):
        return self.thisptr.test(x/self.stepx, y/self.stepy, z/self.stepz)

