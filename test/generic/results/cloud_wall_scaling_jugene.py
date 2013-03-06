import benchmarks
from matplotlib.pylab import *
from numpy import *

testcase = 'cloud_wall_scaling_jugene.xml'

all_charges, all_tolerances, all_cores, timing, _, _ = benchmarks.read(testcase)

def plot_against_cores(testcase, charges, tolerance): 
    ixcha = all_charges.index(charges)
    ixtol = all_tolerances.index(tolerance)

    methods = timing.keys()
    methods.sort()

    rawdata      = zeros([len(all_cores), len(methods)+1])
    rawdata[:,0] = all_cores

    for method in methods:
        data = timing[method]
        if not all(isnan(data[ixcha,ixtol,:])):
            plt.loglog(all_cores, data[ixcha,ixtol,:], 
                       label=method, **benchmarks.fmt(method))
    plt.legend()
    plt.xlabel('#cores')
    plt.ylabel('Time [s]')

figure('8100 charges')
subplot(121, title='tol 1e-3')
plot_against_cores(testcase, 8100, 1e-3)
subplot(122, title='tol 1e-5')
plot_against_cores(testcase, 8100, 1e-5)
tight_layout()

figure('102900 charges')
subplot(121, title='tol 1e-3')
plot_against_cores(testcase, 102900, 1e-3)
subplot(122, title='tol 1e-5')
plot_against_cores(testcase, 102900, 1e-5)
tight_layout()

figure('1012500 charges')
subplot(121, title='tol 1e-3')
plot_against_cores(testcase, 1012500, 1e-3)
subplot(122, title='tol 1e-5')
plot_against_cores(testcase, 1012500, 1e-5)
tight_layout()

figure('9830400 charges')
subplot(121, title='tol 1e-3')
plot_against_cores(testcase, 9830400, 1e-3)
subplot(122, title='tol 1e-5')
plot_against_cores(testcase, 9830400, 1e-5)
tight_layout()

figure('1012500000 charges')
subplot(121, title='tol 1e-3')
plot_against_cores(testcase, 1012500000, 1e-3)
subplot(122, title='tol 1e-5')
plot_against_cores(testcase, 1012500000, 1e-5)
tight_layout()
                        
plt.show()
