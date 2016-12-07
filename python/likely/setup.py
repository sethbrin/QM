from distutils.core import setup,  Extension
from Cython.Distutils import build_ext
import numpy

   
ext = Extension(
    "likely",                 # name of extension
    ["likely.pyx", "TriCubicInterpolator.cpp"],     # filename of our Cython source
    language="c++",              # this causes Cython to create C++ source
    include_dirs=[numpy.get_include()]
    #include_dirs=[r'C:\Python27\Lib\site-packages\numpy\core\include'],          # usual stuff
    # libraries=["stdc++", ...],             # ditto
    # cmdclass = {'build_ext': build_ext}
)
    
setup(name='likely', cmdclass={'build_ext':build_ext}, ext_modules=[ext]) 
