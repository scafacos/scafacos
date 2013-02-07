#!/usr/bin/python

import os
import subprocess

def compute_error(window, alpha):
    """Compute the fast summation error for given alpha."""

    tmp_file="tmp_output.txt"

    line="./scafacos_test p2nfft systems/3d-periodic/" + testcase + ".xml.gz -c pnfft_window_name," + window + ",p2nfft_r_cut," + str(rcut) + ",pnfft_N," + N + ",pnfft_n," + n +",pnfft_m," + m + ",p2nfft_ignore_tolerance,1,p2nfft_alpha," + str(alpha) + ",pnfft_interlaced," + interlaced + ",pnfft_intpol_order," + str(pnfft_intpol_order) + ",p2nfft_intpol_order," + str(p2nfft_intpol_order) + ",pfft_patience_name," + patience
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


testcase="cloud_wall_8100"
N="32,32,32"
n="32,32,32"
m="2"
rcut="3"
interlaced="1"
p2nfft_intpol_order="-1"
pnfft_intpol_order="-1"
patience="estimate"

# windows=['gaussian']
windows=['bspline']

for window in windows:
    pot_file="data_pot_" + window + ".txt"
    field_file="data_field_" + window + ".txt"
    
    ofiles=[pot_file, field_file]
    for fname in ofiles:
      with open(fname, 'w') as f:
        f.write('# ' + window + '\n')

    al = 0.0
    au = 10.0
    fal = compute_error(window, al)
    fau = compute_error(window, au)

    for i in range(20):
      # calculate ml=(2*al+au)/3
      ml = (2.0*al+au) / 3.0
      fml = compute_error(window, ml);
    
      # calculate mu=(al+2*au)/3
      mu = (al+2.0*au) / 3.0
      fmu = compute_error(window, mu)

#       print("al = " + str(al) + ", fal = " + str(fal))
#       print("ml = " + str(ml) + ", fml = " + str(fml))
#       print("mu = " + str(mu) + ", fmu = " + str(fmu))
#       print("au = " + str(au) + ", fau = " + str(fau))

      if fml < fmu :
        au=mu
        fau=fmu
      elif fml > fmu :
        al = ml
        fal = fml
      else :
        au=mu
        fau=fmu
        al=ml
        fal=fml

#       print("new al = " + str(al) + ", fal = " + str(fal))
#       print("new ml = " + str(ml) + ", fml = " + str(fml))
#       print("new mu = " + str(mu) + ", fmu = " + str(fmu))
#       print("new au = " + str(au) + ", fau = " + str(fau))

    print("al = " + str(al) + ", fal = " + str(fal))
    print("ml = " + str(ml) + ", fml = " + str(fml))
    print("mu = " + str(mu) + ", fmu = " + str(fmu))
    print("au = " + str(au) + ", fau = " + str(fau))

    
#    error = compute_error(window, 3)
#    print("Das ist der Fehler:", error)




## use gnuplot with the following lines to plot the data
# plot 'data_pot_bspline.txt', 'data_pot_kaiser.txt', 'data_pot_bessel_i0.txt', 'data_pot_gaussian.txt'
# plot 'data_field_bspline.txt', 'data_field_kaiser.txt', 'data_field_bessel_i0.txt', 'data_field_gaussian.txt'


