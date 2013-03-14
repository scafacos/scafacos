import benchmarks
from matplotlib.pyplot import *
from numpy import *

def plot_timing(charges, tolerance):
    ixcha = all_charges.index(charges)
    ixtol = all_tolerances.index(tolerance)

    methods = timing.keys()
    methods.sort()

    for method in methods:
      data = timing[method]
      if not all(isnan(data[ixcha,ixtol,:])):
        loglog(all_cores, data[ixcha,ixtol,:], 
                   label=method, **benchmarks.fmt(method))

    xlim((1, charges/200.))
    xlabel('#cores')
    ylabel('Time [s]')

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
        semilogx(all_cores, minimal_time/(data[ixcha,ixtol,:]*all_cores), 
                   label=method, **benchmarks.fmt(method))

        #    legend()
    xlim((1, charges/200.))
    xlabel('#cores')
    ylabel('Efficiency')

all_charges, all_tolerances, all_cores, timing = \
    benchmarks.read('silica_melt_scaling_juropa.xml')

sp = 1
for charges in all_charges:
    subplot(2, 3, sp)
    sp += 1
    title('Juropa, {} charges'.format(charges))
    plot_efficiency(charges, 1e-3)

# create legend
methods = []
artists = []
for method in benchmarks.method2format:
    methods.append(method)
    artists.append(
        Line2D((0,0),(0,1),
               **benchmarks.method2format[method]))
    
subplot(234)
legend(artists, methods, loc='lower left',
       bbox_to_anchor=(-0.5, -0.3), ncol=8)
    
show()
