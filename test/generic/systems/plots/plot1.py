import benchmarks
import matplotlib.pylab as plt

for testcase in ['3d-periodic/bm_cloud_wall_juropa.xml']:
    for charges in [8100, 102900, 1012500, 9830400]:
        for tolerance in [1e-3]:
            benchmarks.plot_against_cores(testcase, charges, tolerance, False, "")

plt.show()
