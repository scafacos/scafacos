from matplotlib.pyplot import *
from numpy import *
import benchmarks

matplotlib.rc('font', family='serif', size=18)
matplotlib.rc('text', usetex=True)
matplotlib.rc('axes.formatter', limits=(-3,+3))
matplotlib.rc('legend', fontsize='small')

def plot_timing(charges, tolerance):
    ixcha = all_charges.index(charges)
    ixtol = all_tolerances.index(tolerance)

    # look for the minimal timing
    minimal_time = None
    for method in all_methods:
        methodtimes = timing[method][ixcha,ixtol,:]
        if all(isnan(methodtimes)): continue
        i = 0
        while (isnan(methodtimes[i])): i += 1
        serial_time = all_cores[i] * methodtimes[i]
        if minimal_time is None or serial_time < minimal_time:
            minimal_time = serial_time

    loglog(all_cores, minimal_time/array(all_cores), 'k:')
    loglog(all_cores, 10.*minimal_time/array(all_cores), 'k:')
    loglog(all_cores, 100.*minimal_time/array(all_cores), 'k:')

    for method in all_methods:
      data = timing[method]
      if not all(isnan(data[ixcha,ixtol,:])):
        loglog(all_cores, data[ixcha,ixtol,:], 
                   label=method, **benchmarks.fmt(method))

    xlim((1, charges/200.))
    xscale('log', basex=2)
    xlabel('\#Cores $P$')
    ylabel('Time $t$ [s]')

def plot_efficiency(charges, tolerance):
    ixcha = all_charges.index(charges)
    ixtol = all_tolerances.index(tolerance)
    
    # look for the minimal timing
    minimal_time = None
    for method in all_methods:
        methodtimes = timing[method][ixcha,ixtol,:]
        if all(isnan(methodtimes)): continue
        i = 0
        while (isnan(methodtimes[i])): i += 1
        serial_time = all_cores[i] * methodtimes[i]
        if minimal_time is None or serial_time < minimal_time:
            minimal_time = serial_time

    if minimal_time is None:
        print "No data for charges={} tolerance={}"\
          .format(charges, tolerance)
        return
            
    for method in all_methods:
      data = timing[method]
      if not all(isnan(data[ixcha,ixtol,:])):
        semilogx(all_cores, minimal_time/(data[ixcha,ixtol,:]*all_cores), 
                   label=method, **benchmarks.fmt(method))

        #    legend()
    xlim((1, charges/200.))
    xscale('log', basex=2)
    xlabel(r'\#Cores $P$')
    ylabel(r'Relative Parallel Efficiency $e(P)$')

def create_legend():
  # create legend
  methods = []
  artists = []
  for method in benchmarks.method2format:
    methods.append(method)
    artists.append(
        Line2D((0,0),(0,1),
               **benchmarks.method2format[method]))
  legend()
    

####################
# COMPLEXITY
####################

data = benchmarks.read('complexity.xml')
all_methods = data['methods']
all_cores = data['cores']
all_tolerances = data['tolerances']
all_charges = data['charges']
timing = data['timing']

cores = 1
ixcor = all_cores.index(cores)
tolerance = 1e-3
ixtol = all_tolerances.index(tolerance)

figure('Complexity')
for method in data['methods']:
    time_per_particle = \
      timing[method][:, ixtol, ixcor] / all_charges
    loglog(all_charges, time_per_particle,
           label=method, **benchmarks.fmt(method))

xlabel(r'\#Charges')
ylabel(r'Time $t$/\#Charges [s]')
create_legend()
savefig('complexity.pdf')

####################
# ACCURACY
####################
data = benchmarks.read('accuracy.xml')
all_methods = data['methods']
all_cores = data['cores']
all_tolerances = data['tolerances']
all_charges = data['charges']
timing = data['timing']

figure('Accuracy')
cores = 1
ixcor = all_cores.index(cores)
for charges in data['charges']:
    ixcha = data['charges'].index(charges)
    for method in all_methods:
        loglog(all_tolerances, 
               timing[method][ixcha, :, ixcor]/charges, 
               label=method, **benchmarks.fmt(method))

xlabel(r'Relative RMS potential error $\varepsilon_\mathrm{pot}$')
ylabel(r'Time $t$/\#Charges [s]')
gca().xaxis.set_ticks([1e-13,1e-11,1e-9,1e-7,1e-5,1e-3,1e-1])
create_legend()
savefig('accuracy.pdf')

####################
# CLOUD WALL SCALING JUROPA
####################
data = benchmarks.read('cloud_wall_scaling_juropa.xml')
all_methods = data['methods']
all_cores = data['cores']
all_tolerances = data['tolerances']
all_charges = data['charges']
timing = data['timing']

figure('Timing Juropa 1012500 charges')
plot_timing(1012500, 1e-3)
create_legend()
savefig('timing-juropa-1012500.pdf')

figure('Juropa 8100 charges')
plot_efficiency(8100, 1e-3)
create_legend()
savefig('eff-juropa-8100.pdf')

figure('Juropa 102900 charges')
plot_efficiency(102900, 1e-3)
create_legend()
savefig('eff-juropa-102900.pdf')

figure('Juropa 1012500 charges')
plot_efficiency(1012500, 1e-3)
create_legend()
savefig('eff-juropa-1012500.pdf')

####################
# CLOUD WALL SCALING JUGENE
####################
data = benchmarks.read('cloud_wall_scaling_jugene.xml')
all_methods = data['methods']
all_cores = data['cores']
all_tolerances = data['tolerances']
all_charges = data['charges']
timing = data['timing']

figure('JuGene 1012500 charges')
plot_efficiency(1012500, 1e-3)
create_legend()
savefig('eff-jugene-1012500.pdf')

figure('JuGene 9830400 charges')
plot_efficiency(9830400, 1e-3)
create_legend()
savefig('eff-jugene-9830400.pdf')

figure('JuGene 1012500000 charges')
plot_efficiency(1012500000, 1e-3)
create_legend()
savefig('eff-jugene-1012500000.pdf')

#tight_layout()
                        
show()
