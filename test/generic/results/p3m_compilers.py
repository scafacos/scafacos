import benchmarks
import matplotlib.pylab as plt

for testcase in ['p3m_compilers.xml']:
    for charges in [102900]:
        for tolerance in [1e-3]:
            benchmarks.plot_against_cores(testcase, charges, tolerance)

plt.show()
