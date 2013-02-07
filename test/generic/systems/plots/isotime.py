import benchmarks
import matplotlib.pylab as plt
from matplotlib.image import NonUniformImage
import numpy

tolerance = 1e-3
interpolation='bicubic'
levels = numpy.arange(-4,4)

all_charges, all_tolerances, all_cores, timing, speedup, efficiency = \
  benchmarks.read('3d-periodic/bm_cloud_wall_jugene.xml')

ixtol = all_tolerances.index(tolerance)
charges = numpy.array(all_charges)[::-1]
cores = numpy.array(all_cores)
Z = numpy.log(timing['p2nfft'][:,ixtol,:])

plt.subplot(221)
plt.imshow(Z, origin='lower', interpolation='nearest')
cbar = plt.colorbar(ticks=levels)
cbar.ax.set_yticklabels(10.0**levels)
plt.xticks(numpy.arange(len(cores)), cores)
plt.xlabel('#cores')
plt.yticks(numpy.arange(len(charges)), charges[::-1])
plt.ylabel('#charges')
plt.grid(True)

plt.subplot(223)
plt.contourf(Z, levels=levels, origin='lower', interpolation=interpolation)
cbar = plt.colorbar(ticks=levels)
cbar.ax.set_yticklabels(10.0**levels)
plt.xticks(numpy.arange(len(cores)), cores)
plt.xlabel('#cores')
plt.yticks(numpy.arange(len(charges)), charges[::-1])
plt.ylabel('#charges')
plt.grid(True)

Z = numpy.log(timing['fmm'][:,ixtol,:])
plt.subplot(222)
plt.imshow(Z, origin='lower', interpolation='nearest')
cbar = plt.colorbar(ticks=levels)
cbar.ax.set_yticklabels(10.0**levels)
plt.xticks(numpy.arange(len(cores)), cores)
plt.xlabel('#cores')
plt.yticks(numpy.arange(len(charges)), charges[::-1])
plt.ylabel('#charges')
plt.grid(True)

plt.subplot(224)
plt.contourf(Z, levels=levels, origin='lower', interpolation=interpolation)
cbar = plt.colorbar(ticks=levels)
cbar.ax.set_yticklabels(10.0**levels)
plt.xticks(numpy.arange(len(cores)), cores)
plt.xlabel('#cores')
plt.yticks(numpy.arange(len(charges)), charges[::-1])
plt.ylabel('#charges')
plt.grid(True)

# ax = plt.subplot(313)
# im = NonUniformImage(ax, extent=(1,10000000,1,1000000))
# im.set_data(cores, charges, Z)
# ax.images.append(im)

plt.show()
