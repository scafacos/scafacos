import benchmarks
import matplotlib.pylab as plt

for testcase in ['cloud_wall_scaling_juropa.xml']:
    for charges in [8100, 102900, 1012500, 9830400]:
        for tolerance in [1e-3, 1e-5]: 
            benchmarks.plot_against_cores(testcase, charges, tolerance, 
                                          plot_efficiency=False)

plt.show()
