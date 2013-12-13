import mpi4py.MPI
import scafacos
import numpy

s = scafacos.scafacos("p3m")

positions = numpy.array([[0.5, 0.5], [2.0, 2.0]], dtype='double')
charges = numpy.array([1.0, -1.0], dtype='double')
s.tune(positions, charges)
