#!/usr/bin/python

import benchmarks
import matplotlib.pylab as plt
import numpy as np
import os as os

#              [filename                              ,  [ particlenumbers                                     ]
testcases = [  ['3d-periodic/bm_cloud_wall_jugene.xml',  [8100, 102900, 1012500, 9830400, 102900000, 1012500000]]
              ,['3d-periodic/bm_cloud_wall_juropa.xml',  [8100, 102900, 1012500, 9830400, 102900000, 1012500000]]
              ,['3d-periodic/bm_silica_melt_jugene.xml', [12960, 103680, 829440, 9447840, 103680000, 1030410720]]
              ,['3d-periodic/bm_silica_melt_juropa.xml', [12960, 103680, 829440, 9447840, 103680000, 1030410720]]
              ,['nonperiodic/bm_cloud_wall_jugene.xml',  [8100, 102900, 1012500, 9830400, 102900000, 1012500000]]
              ,['nonperiodic/bm_cloud_wall_juropa.xml',  [8100, 102900, 1012500, 9830400, 102900000, 1012500000]]
              ,['nonperiodic/bm_silica_melt_jugene.xml', [12960, 103680, 829440, 9447840, 103680000, 1030410720]]
              ,['nonperiodic/bm_silica_melt_juropa.xml', [12960, 103680, 829440, 9447840, 103680000, 1030410720]]
              ,['nonperiodic/bm_galaxies_jugene.xml',    [213671, 1709368, 13674944, 109399552]]
              ,['nonperiodic/bm_galaxies_juropa.xml',    [213671]]
              ,['1d-periodic/bm_cloud_wall_juropa.xml',  [1012500]]
              ,['2d-periodic/bm_cloud_wall_juropa.xml',  [1012500]]
            ]

subdir     = "./plot_all/"
corelist   = 2**np.array(range(0,13))
tolerances = [1e-3, 1e-5]

# create target directory if it does not exist
if not os.path.exists(subdir):
    os.makedirs(subdir)

# plot runtime depending on total particles for different processor numbers
for testcase in testcases:
    case       = testcase[0]
    charges    = testcase[1]

    for charge in charges:
        for tolerance in tolerances:
	    print "\n\n Testcase: %s, Charges: %d, Tolerance: %e" % (case, charge, tolerance)
            benchmarks.plot_against_cores(case, charge, tolerance, True, subdir)

# plot runtime depending on number of processors for different particle numbers
for testcase in testcases:
    case       = testcase[0]

    for cores in corelist:
        for tolerance in tolerances:
	    print "\n\n Testcase: %s, Cores: %d, Tolerance: %e" % (case, cores, tolerance)
            benchmarks.plot_against_charges(case, cores, tolerance, True, subdir)

#plt.show()

print "Done. Created pdf/png-files of all plots in directory %s." % subdir
