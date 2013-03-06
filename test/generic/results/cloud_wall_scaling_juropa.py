import benchmarks
from matplotlib.pylab import *

testcase = 'cloud_wall_scaling_juropa.xml'


figure()

subplot(321, title='8100 charges, tol 1e-3')
benchmarks.plot_against_cores(testcase, 8100, 1e-3, plot_efficiency=False)
subplot(322, title='8100 charges, tol 1e-5')
benchmarks.plot_against_cores(testcase, 8100, 1e-5, plot_efficiency=False)

subplot(323, title='102900 charges, tol 1e-3')
benchmarks.plot_against_cores(testcase, 102900, 1e-3, plot_efficiency=False)
subplot(324, title='102900 charges, tol 1e-5')
benchmarks.plot_against_cores(testcase, 102900, 1e-5, plot_efficiency=False)

subplot(325, title='1012500 charges, tol 1e-3')
benchmarks.plot_against_cores(testcase, 1012500, 1e-3, plot_efficiency=False)
subplot(326, title='1012500 charges, tol 1e-5')
benchmarks.plot_against_cores(testcase, 1012500, 1e-5, plot_efficiency=False)

plt.show()
