from __future__ import print_function
import mpi4py.MPI
import scafacos
import numpy

p3m = scafacos.scafacos("p3m")
p3m.box = numpy.array([[2.0, 0.0, 0.0],
                        [0.0, 2.0, 0.0],
                        [0.0, 0.0, 2.0]])
p3m.periodicity = (True, True, True)
p3m.near_field_flag = True
p3m.total_particles = 2

positions = numpy.array([[0.5, 1.0, 1.0], [1.5, 1.0, 1.0]], dtype='double')
print(positions)
charges = numpy.array([1.0, -1.0], dtype='double')

print(positions.reshape(6))

p3m.tune(positions, charges)
fields, potentials = p3m(positions, charges)

print("fields={}".format(fields))
print("potentials={}".format(potentials))
