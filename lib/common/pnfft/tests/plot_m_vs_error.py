#!/usr/bin/python

import subprocess

def compute_error(window, m, fname):
    """Compute the relative NFFT error for given window and cutoff parameter 'm'."""

    line="mpirun -np 8 pnfft_check_vs_pfft -pnfft_m " + str(m) + " -pnfft_window " + str(window)
#     print(line)
    print m

    output = subprocess.check_output(line, stderr=subprocess.STDOUT, shell=True)
    
    matched_lines = [line for line in output.split('\n') if "relative maximum error" in line]
    error_rel = float(str(matched_lines[0]).split('=')[1])
#     print("error_rel = " + str(error_rel))

    with open(fname, 'a') as f:
      f.write(str(m) + "  " + str(error_rel) + "\n")

    return error_rel


windows=['PNFFT_WINDOW_KAISER_BESSEL', 'PNFFT_WINDOW_GAUSSIAN', 'PNFFT_WINDOW_BSPLINE', 'PNFFT_WINDOW_SINC_POWER', 'PNFFT_WINDOW_BESSEL_I0']

for i in range(5):
    fname="data_" + windows[i] + "_vs_error.txt"
    
    with open(fname, 'w') as f:
      f.write('# ' + windows[i] + '\n')

    print "Compute errors for window: ", windows[i]
    print "m="
#     for m in range(1,21):
    for m in range(1,21):
      err = compute_error(i, m, fname)



## use gnuplot with the following lines to plot the data
# gnuplot
# set logscale y
# plot 'data_PNFFT_WINDOW_GAUSSIAN_vs_error.txt', 'data_PNFFT_WINDOW_BSPLINE_vs_error.txt', 'data_PNFFT_WINDOW_SINC_POWER_vs_error.txt', 'data_PNFFT_WINDOW_KAISER_BESSEL_vs_error.txt', 'data_PNFFT_WINDOW_BESSEL_I0_vs_error.txt'

## use this to create pdf in gnuplot:
# set terminal postscript enhanced color
# set output '| ps2pdf - plot.pdf'
# plot ...
