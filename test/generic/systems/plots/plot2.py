import benchmarks
import matplotlib.pylab as plt

for testcase in ['3d-periodic/bm_cloud_wall_jugene.xml']:
    for cores in [32, 64]:
        for tolerance in [1e-3, 1e-5]:
            benchmarks.plot_against_charges(testcase, cores, tolerance, False, "")

plt.show()
