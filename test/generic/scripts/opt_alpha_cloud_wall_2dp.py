#!/usr/bin/python

import os
import subprocess
from scipy.optimize import minimize_scalar

def compute_error(alpha):
    """Compute the fast summation error for given alpha."""

    tmp_file="tmp_output.txt"

    line="./scafacos_test p2nfft systems/2d-periodic/" + testcase + ".xml.gz -c pnfft_window_name," + window + ",p2nfft_r_cut," + str(rcut) + ",pnfft_N," + N + ",pnfft_n," + n +",pnfft_m," + m + ",p2nfft_ignore_tolerance,1,p2nfft_alpha," + str(alpha) + ",pnfft_interlaced," + interlaced + ",pnfft_intpol_order," + str(pnfft_intpol_order) + ",p2nfft_intpol_order," + str(p2nfft_intpol_order) + ",pfft_patience_name," + patience + ",tolerance_field," + tolerance + ",p2nfft_p," + str(p2nfft_p) + ",p2nfft_epsB," + str(p2nfft_epsB)

    print(line)

#     with open(tmp_file, 'r+a') as f:
#       subprocess.check_call(line, stdout=f, stderr=subprocess.STDOUT, shell=True)
#       output = f.read()

    output = subprocess.check_output(line, stderr=subprocess.STDOUT, shell=True)
#     print("output = ")
#     print(output)
    
    matched_lines = [line for line in output.split('\n') if "rel_rms_potential_error" in line]
#     print("matched_lines = ")
#     print(matched_lines)
    error_pot = float(str(matched_lines[0]).split('=')[1])
    print("alpha = " + str(alpha) + ", error_pot = " + str(error_pot))

    matched_lines = [line for line in output.split('\n') if "rel_rms_field_error" in line]
    error_field = float(str(matched_lines[0]).split('=    ')[1])
    print("alpha = " + str(alpha) + ", error_field = " + str(error_field))
   
    with open(pot_file, 'a') as f:
      f.write(str(alpha) + "  " + str(error_pot) + "\n")

    with open(field_file, 'a') as f:
      f.write(str(alpha) + "  " + str(error_field) + "\n")


    return error_pot


testcase="cloud_wall_300"
# testcase="cloud_wall_8100"
# testcase="cloud_wall_102900"
# N0=str(128)
N0=str(64)
N1=str(64)
# p2nfft_epsB="0.015625"
# p2nfft_epsB="0.003125"
p2nfft_epsB="0.125"
# p2nfft_epsB="0.0625"
N=N0+","+N1+","+N1
n=str(2*int(N0))+","+N1+","+N1
# n=str(1.125*int(N0))+","+N1+","+N1
m="4"
rcut="2.5"
p2nfft_p="8"

interlaced="0"
p2nfft_intpol_order="-1"
pnfft_intpol_order="-1"
patience="estimate"
tolerance="1e-10"

# windows=['gaussian']
windows=['bspline']

for window in windows:
    pot_file="data_pot_" + window + ".txt"
    field_file="data_field_" + window + ".txt"
    
    ofiles=[pot_file, field_file]
    for fname in ofiles:
      with open(fname, 'w') as f:
        f.write('# ' + window + '\n')

    res = minimize_scalar(compute_error, bracket=(0, 10), method='brent', tol=1e-3)
    print[res.x]




## use gnuplot with the following lines to plot the data
# plot 'data_pot_bspline.txt', 'data_pot_kaiser.txt', 'data_pot_bessel_i0.txt', 'data_pot_gaussian.txt'
# plot 'data_field_bspline.txt', 'data_field_kaiser.txt', 'data_field_bessel_i0.txt', 'data_field_gaussian.txt'


