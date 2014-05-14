#!/usr/bin/env python

# To build the scafacos python interface
# * build ScaFaCoS with the configure option "--with-pic"
# * first *install* the ScaFaCoS library
# * make sure pkg-config can find scafacos (execute "pkg-config --libs scafacos")
# * if not, add /path/to/scafacos/lib/pkgconfig to the environment
#   variable PKG_CONFIG_PATH
# * run "python setup.py build_ext -fi" to build the python module

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import os
import commands

def pkgconfig(package, **kw):
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    cmd = "pkg-config --libs --cflags {}".format(package)
    status, output = commands.getstatusoutput(cmd)
    if status != 0:
        raise Exception(output)
    tokens = output.split()
    for token in tokens:
        key = token[:2]
        value = token[2:]
        if flag_map.has_key(key):
            kw.setdefault(flag_map[key], []).append(value)
        else: # throw others to extra_link_args
            kw.setdefault('extra_link_args', []).append(token)
    return kw

scafacos_params = pkgconfig("scafacos")
scafacos_params.setdefault('include_dirs', []).append(numpy.get_include())

print(scafacos_params.setdefault('libraries', []))

ext_modules=[
    Extension("scafacos", ["scafacos.pyx"], **scafacos_params),
]

setup(
    name = 'scafacos',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
