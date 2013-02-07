import benchmarks
import matplotlib.pylab as plt

for testcase in ['../3d-periodic/bm_icc_gcc.xml']:
    for charges in [102900]:
        for tolerance in [1e-3]:
            benchmarks.plot_against_cores(testcase, charges, tolerance)

plt.show()
