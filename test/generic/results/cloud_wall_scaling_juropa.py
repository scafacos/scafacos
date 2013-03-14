import benchmarks
from matplotlib.pylab import *
from numpy import *

testcase = 'cloud_wall_scaling_juropa.xml'

all_charges, all_tolerances, all_cores, timing = \
    benchmarks.read(testcase)

def plot_timing(charges, tolerance):
    ixcha = all_charges.index(charges)
    ixtol = all_tolerances.index(tolerance)

    methods = timing.keys()
    methods.sort()

    for method in methods:
      data = timing[method]
      if not all(isnan(data[ixcha,ixtol,:])):
        plt.loglog(all_cores, data[ixcha,ixtol,:], 
                   label=method, **benchmarks.fmt(method))

    plt.legend()
    plt.xlim((1, charges/200.))
    plt.xlabel('#cores')
    plt.ylabel('Time [s]')

def plot_efficiency(charges, tolerance):
    ixcha = all_charges.index(charges)
    ixtol = all_tolerances.index(tolerance)

    methods = timing.keys()
    methods.sort()

    minimal_time = None
    # look for the minimal timing
    for method in methods:
        methodtimes = timing[method][ixcha,ixtol,:]
        if all(isnan(methodtimes)): continue
        i = 0
        while (isnan(methodtimes[i])): i += 1
        serial_time = all_cores[i] * methodtimes[i]
        if minimal_time is None or serial_time < minimal_time:
            minimal_time = serial_time

    if minimal_time is None:
        print "No data for charges={} tolerance={}".format(charges, tolerance)
        return
            
    for method in methods:
      data = timing[method]
      if not all(isnan(data[ixcha,ixtol,:])):
        plt.semilogx(all_cores, minimal_time/(data[ixcha,ixtol,:]*all_cores), 
                   label=method, **benchmarks.fmt(method))

        #    plt.legend()
    plt.xlim((1, charges/200.))
    plt.xlabel('#cores')
    plt.ylabel('Efficiency')

figure()
# subplot(221, title='8100 charges, tol 1e-3')
# plot_timing(8100, 1e-3)
# subplot(222, title='8100 charges, tol 1e-5')
# plot_timing(8100, 1e-5)

# subplot(231, title='8100 charges, tol 1e-3')
# plot_efficiency(8100, 1e-3)
# plt.legend()
# subplot(234, title='8100 charges, tol 1e-5')
# plot_efficiency(8100, 1e-5)

subplot(221, title='1012500 charges, tol 1e-3')
plt.loglog(all_cores, 10./array(all_cores), 'k--')
plt.loglog(all_cores, 500./array(all_cores), 'k--')
plot_timing(1012500, 1e-3)
plt.legend()

subplot(222, title='8100 charges, tol 1e-3')
plot_efficiency(8100, 1e-3)

subplot(223, title='102900 charges, tol 1e-3')
plot_efficiency(102900, 1e-3)

subplot(224, title='1012500 charges, tol 1e-3')
plot_efficiency(1012500, 1e-3)




#figure()
# subplot(221, title='102900 charges, tol 1e-3')
# plot_timing(102900, 1e-3)
# subplot(222, title='102900 charges, tol 1e-5')
# plot_timing(102900, 1e-5)

#figure()
# subplot(221, title='1012500 charges, tol 1e-3')
# plot_timing(1012500, 1e-3)
# subplot(222, title='1012500 charges, tol 1e-5')
# plot_timing(1012500, 1e-5)


show()
