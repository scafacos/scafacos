import benchmarks
from matplotlib.pyplot import *

all_charges, all_tolerances, all_cores, \
        timing, speedup, efficiency = \
        benchmarks.read('complexity.xml')

methods = timing.keys()
methods.sort()

cores = 1
ixcor = all_cores.index(cores)
tolerance = 1e-3
ixtol = all_tolerances.index(tolerance)

for method in methods:
    loglog(all_charges, timing[method][:, ixtol, ixcor], 'o-', label=method)

xlabel('#charges')
ylabel('Time [s]')
legend()
        
show()
