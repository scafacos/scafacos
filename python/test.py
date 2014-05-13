import mpi4py.MPI
import scafacos
import numpy

p3m = scafacos.scafacos("p3m")
print p3m.method

p3m.box = numpy.array([[2.0, 0.0, 0.0],
                        [0.0, 2.0, 0.0],
                        [0.0, 0.0, 2.0]])
print p3m.box

p3m.periodicity = (True, True, True)
print p3m.periodicity

p3m.near_field_flag = False

positions = numpy.array([[1.0, 1.0], [0.5, 1.5], [1.0, 1.0]], dtype='double')
charges = numpy.array([1.0, -1.0], dtype='double')
p3m.tune(positions, charges)

fields, potentials = p3m(positions, charges)
