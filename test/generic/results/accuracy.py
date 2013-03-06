import benchmarks
from matplotlib.pyplot import *

all_charges, all_tolerances, all_cores, \
        timing, speedup, efficiency = \
        benchmarks.read('accuracy.xml')

methods = timing.keys()
methods.sort()

cores = 1
ixcor = all_cores.index(cores)
for charges in all_charges:
    figure('{} charges'.format(charges))
    ixcha = all_charges.index(charges)
    for method in methods:
        loglog(all_tolerances, timing[method][ixcha, :, ixcor], 
               label=method, **benchmarks.fmt(method))

xlabel('tolerance')
ylabel('Time [s]')
legend()
        
show()
