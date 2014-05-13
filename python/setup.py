#!/usr/bin/env python

# To build the scafacos python wrapper,
# * build ScaFaCoS with the configure option "--enable-single-lib --with-pic"
# * adapt the following paths (builddir might be the same as srcdir)
srcdir = '/home/olenz/projects/scafacos/src'
builddir = '/home/olenz/projects/scafacos/obj.lancre/p3m-debug'

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import os

ext_params = {}
ext_params['include_dirs'] = ['{}/src'.format(srcdir), builddir, numpy.get_include()] 
ext_params['library_dirs'] = ['{}/package/.libs'.format(builddir), '/usr/lib64']
ext_params['libraries'] = ['fcs', 'mpi', 'mpi_cxx', 'fftw3', 'gsl', 'gslcblas']

ext_modules=[
    Extension("scafacos", ["scafacos.pyx"], **ext_params),
]

setup(
    name = 'scafacos',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
