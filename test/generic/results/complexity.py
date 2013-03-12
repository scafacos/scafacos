import benchmarks
from matplotlib.pyplot import *

all_charges, all_tolerances, all_cores, timing = \
        benchmarks.read('complexity.xml')

methods = timing.keys()
methods.sort()

cores = 1
ixcor = all_cores.index(cores)
tolerance = 1e-3
ixtol = all_tolerances.index(tolerance)

figure()
for method in methods:
    loglog(all_charges, timing[method][:, ixtol, ixcor], 
           label=method, **benchmarks.fmt(method))

title('Complexity')
xlabel('#Charges')
ylabel('Time [s]')
legend()

figure()
for method in methods:
    time_per_particle = timing[method][:, ixtol, ixcor] / all_charges
    loglog(all_charges, time_per_particle,
           label=method, **benchmarks.fmt(method))

title('Complexity')
xlabel('#Charges')
ylabel('Time/Charge [s]')
legend()
        
show()
