#!/usr/bin/python

import sys, getopt
import os
import subprocess
from scipy.special import lambertw
from scipy.optimize import minimize_scalar
from math import *
import argparse

# from math import sqrt,ceil,pi

# extrapolate linear dependency
def estimate_time_N(t0, N0, N):
  return (t0 * N) / N0

# extrapolate linear times log dependency
def estimate_time_NlogN(t0, N0, N):
  return (t0 * N * log(N)) / (N0 * log(N0))

# minimize_scalar needs a scalar function,
# therefore, we set up all the other variables as globals
def estimate_time_onearg(rc):
  return estimate_time(t0_near, t0_far_M, t0_far_MlogM, rc_start, M0_start, rc)

def estimate_time(t0_near, t0_far_M, t0_far_MlogM, rc_start, M0_start, rc):
  # set M in dependence of r to fulfill the error bound
  alpha, kc, M, M0, M1, M2, m0, m1, m2, h, epsB = parameter_heuristic(
      rc, L, Q2, N, tol, tol_type, split_tol, P, m, osp, osn)

  t_near = estimate_time_N(t0_near, pow(rc_start,3.0), pow(rc,3.0))
  t_far = estimate_time_N(t0_far_M, pow(M0_start,3.0), pow(M0,3.0)) + estimate_time_NlogN(t0_far_MlogM, pow(M0_start,3.0), pow(M0,3.0))
  t_sum = t_near + t_far

  print("Estimator: "
      + "rc = " + str('%.4e, '%rc)
      + "M0 = " + str('%4d, '%M0)
      + "t_near = " + str('%.4e, '%t_near)
      + "t_far = " + str('%.4e, '% t_far)
      + "t_sum = " + str('%.4e, '% t_sum))
  return t_sum


def filter_p2nfft_runtimes(output):
  rt_tune = search_val(output, "Time:", ":", "")
  rt_all  = search_val(output, "Average time:", ":", "")
  rt_near = search_val(output, "Near field computation takes", "takes", "s")
  rt_far  = search_val(output, "Far field computation takes", "takes", "s")
  return rt_tune, rt_all, rt_near, rt_far 


def filter_pnfft_runtimes(output):
  idx = str(int(search_val(output, "np_pnfft(", "", ",")))
  rt_D  = search_val(output, "pnfft_trf_matrix_D("+idx+")   =", "", ";")
  rt_D += search_val(output, "pnfft_adj_matrix_D("+idx+")   =", "", ";")
  rt_F  = search_val(output, "pnfft_trf_matrix_F("+idx+")   =", "", ";")
  rt_F += search_val(output, "pnfft_adj_matrix_F("+idx+")   =", "", ";")
  rt_C  = search_val(output, "pnfft_trf_loop_B("+idx+")     =", "", ";")
  rt_C += search_val(output, "pnfft_adj_loop_B("+idx+")     =", "", ";")
  rt_G  = search_val(output, "pnfft_trf_gcells("+idx+")     =", "", ";")
  rt_G += search_val(output, "pnfft_adj_gcells("+idx+")     =", "", ";")
  return rt_D, rt_F, rt_C, rt_G 



def filter_pfft_runtimes(output):
  idx = str(int(search_val(output, "np_pfft(", "", ",")))

  iter_forw =  int(search_val(output, "pfft_forw_iter("+idx+")    =", "", ";"))

  # PFFT is executed several times for ik-differentiation 
  iter_forw =  int(search_val(output, "pfft_forw_iter("+idx+")    =", "", ";"))
  iter_back =  int(search_val(output, "pfft_back_iter("+idx+")    =", "", ";"))

  # FFTs of size M0 * M1 * M2 -> M0 * M1 * m2
  t0  = iter_back * search_val(output, "pfft_back_trafo1("+idx+", 2)   =", "", ";")
  t0 += iter_forw * search_val(output, "pfft_forw_trafo6("+idx+", 2)   =", "", ";")

  # FFTs of size M0 * M1 * m2 -> M0 * m1 * m2
  t1  = iter_back * search_val(output, "pfft_back_trafo2("+idx+", 2)   =", "", ";")
  t1 += iter_forw * search_val(output, "pfft_forw_trafo5("+idx+", 2)   =", "", ";")

  # FFTs of size M0 * m1 * m2 -> m0 * m1 * m2
  t2  = iter_back * search_val(output, "pfft_back_trafo3("+idx+", 2)   =", "", ";")
  t2 += iter_forw * search_val(output, "pfft_forw_trafo4("+idx+", 2)   =", "", ";")

  # altough twiddling operates on different array sizes, the scaling factor remains the same
  t3  = iter_back * search_val(output, "pfft_back_itwiddle("+idx+")   =", "", ";")
  t3 += iter_back * search_val(output, "pfft_back_otwiddle("+idx+")   =", "", ";")
  t3 += iter_forw * search_val(output, "pfft_forw_itwiddle("+idx+")   =", "", ";")
  t3 += iter_forw * search_val(output, "pfft_forw_otwiddle("+idx+")   =", "", ";")

  t4  = iter_back * search_val(output, "pfft_back_remap1("+idx+", 1) =", "", ";")
  t4 += iter_back * search_val(output, "pfft_back_remap2("+idx+", 1) =", "", ";")
  t4 += iter_forw * search_val(output, "pfft_forw_remap3("+idx+", 1) =", "", ";")
  t4 += iter_forw * search_val(output, "pfft_forw_remap4("+idx+", 1) =", "", ";")

  t5  = iter_back * search_val(output, "pfft_back_remap_3dto2d("+idx+", 2)   =", "", ";")
  t5 += iter_forw * search_val(output, "pfft_forw_remap_2dto3d("+idx+", 2)   =", "", ";")

  tsum = t0 + t1 + t2 + t3 + t4 + t5

  return t0, t1, t2, t3, t4, t5, tsum


def parameter_heuristic(rc, L, Q2, N, tol, tol_type, split_tol, P, m, osp, osn):
  # invoke Kolafa/Perram parameter tuning
  alpha = inv_kolper_alpha(rc, L, Q2, N, tol * split_tol, tol_type)
  kc    = inv_kolper_kc(alpha, L, Q2, N, tol * split_tol, tol_type)
  M     = kc2M(kc)
#   print("Kolafa-Perram Force = " + str(kolper_near_force(rc, L, Q2, alpha, N)))
#   print("alpha = "+ str('%.2f' % alpha))
#   print("kc = "+ str('%.2f' % kc))

  # M = 64
  # alpha = 0.7
  # alpha = 0.748120

  kc0 = kc
  if nd == "2d":
    kc0 = inv_kolper_kc(alpha, 2*L, Q2, N, tol * split_tol, tol_type)
  if nd == "1d":
    kc0 = inv_kolper_kc(alpha, sqrt(2)*2*L, Q2, N, tol * split_tol, tol_type)
  if nd == "0d":
    kc0 = inv_kolper_kc(alpha, sqrt(3)*2*L, Q2, N, tol * split_tol, tol_type)

  M0 = M1 = M2 = M
  m0 = m1 = m2 = oversampled_gridsize(M, m, osp)
  if nd == "2d":
    M0 = even( kc2M(kc0) + P )
    m0 = oversampled_gridsize(M0, m, osn)
  if nd == "1d":
    M0 = M1 = even( kc2M(kc0) + P )
    m0 = m1 = oversampled_gridsize(M0, m, osn)
  if nd == "0d":
    M0 = M1 = M2 = even( kc2M(kc0) + P )
    m0 = m1 = m2 = oversampled_gridsize(M0, m, osn)

  h = L*M0/M

  epsB = 0.5*P/M0 if nd != "3d" else 0.0

  # Compare correctness of near- and farfield error formulae.
  # Therefore, choose either large rc or large kc.
  # M = 2*M; kc = M/2.0-1.0
  # rc = 1.5*rc

  return alpha, kc, M, M0, M1, M2, m0, m1, m2, h, epsB 



split_tol = 0.5
# split_tol = 1.0/sqrt(2.0)

def is_in_text(text, name):
  if text == "":
    return 0
  
  matched_lines = [line for line in text if name in line]
  if matched_lines==[]:
    return 0
  else:
    return 1


def search_val(text, name, pre_delim, post_delim):
  if text == "":
    return 0.0
  else:
    matched_lines = [line for line in text if name in line]
    if pre_delim == "":
      pre_delim = name
    if post_delim == "":
      return float(str(matched_lines[0]).split(pre_delim)[1])
    else:
      return float(str(matched_lines[0]).split(pre_delim)[1].split(post_delim)[0])

def even(N):
  return int(2 * ceil(N / 2.0))

def oversampled_gridsize(M, m, flag):
  # if flag is a float it is interpreted as oversampling factor
  if str(flag) != str(int(flag)):
    return even(M*flag)

  # if flag is an integer it determines the way of constant oversampling
  if flag == -1:
    return M + 2*m+2
  if flag == -2:
    return 2*M
  else:
    return M + flag

# overwrite default parameters, if command line arguments are given
def set_str(args, nr, default):
  if len(args) > nr:
    if args[nr] >= 0:
      return str(args[nr])
  return default

def set_int(args, nr, default):
  if len(args) > nr:
    if args[nr] >= 0:
      return int(args[nr])
  return default

def set_float(args, nr, default):
  if len(args) > nr:
    if args[nr] > 0:
      return float(args[nr])
  return default

def set_number(string):
  try:
    n = float(string) if '.' in string else int(string)
  except:
    print("Error in set_number: string is neither float nor int")
    sys.exit()
  return n

def inv_kolper_potential_alpha(rc, L, Q2, tol):
#   if rc < 0.0:
#     return 1e20
  lamb = lambertw( 1.0 / tol * sqrt(Q2 * rc / L**3) )
  return 1.0/rc * sqrt(lamb.real)
def inv_kolper_field_alpha(rc, L, Q2, tol):
  l = log( 2.0 / tol * sqrt( Q2 / rc / L**3) )
  return 1.0/rc * sqrt(l)

def inv_kolper_potential_kc(alpha, L, Q2, tol):
  lamb = lambertw( 4.0 / 3.0 / L**2 * pow( sqrt(2.0 * Q2 / pi / alpha) / tol, 4.0/3.0) )
  return sqrt(3.0) / 2.0 * alpha * L / pi * sqrt(lamb.real)
def inv_kolper_field_kc(alpha, L, Q2, tol):
  lamb = lambertw(pow( sqrt(Q2 * alpha) / sqrt(pi * L**3) * 8.0 / tol , 4.0))
  return 0.5 * alpha * L / pi * sqrt(lamb.real)

def inv_kolper_alpha(rc, L, Q2, N, tol, tol_type):
  if tol_type == "phi":
    return inv_kolper_potential_alpha(rc, L, Q2, tol)
  elif tol_type == "U":
    return inv_kolper_potential_alpha(rc, L, Q2, tol * sqrt(N/Q2))
  elif tol_type == "E":
    return inv_kolper_field_alpha(rc, L, Q2, tol)
  elif tol_type == "F":
    return inv_kolper_field_alpha(rc, L, Q2, tol * sqrt(N/Q2))
  else:
    print("ERROR: tol-type out of range")
    sys.exit()

def inv_kolper_kc(alpha, L, Q2, N, tol, tol_type):
  if tol_type == "phi":
    return inv_kolper_potential_kc(alpha, L, Q2, tol)
  elif tol_type == "U":
    return inv_kolper_potential_kc(alpha, L, Q2, tol * sqrt(N/Q2))
  elif tol_type == "E":
    return inv_kolper_field_kc(alpha, L, Q2, tol)
  elif tol_type == "F":
    return inv_kolper_field_kc(alpha, L, Q2, tol * sqrt(N/Q2))
  else:
    print("ERROR: tol-type out of range")
    sys.exit()

def kolper_near_potential(rc, L, Q2, alpha):
  arg = pow(alpha*rc, 2.0)
  return sqrt(Q2 * rc / L**3) * exp(-arg) / arg
def kolper_near_field(rc, L, Q2, alpha):
  arg = pow(alpha*rc, 2.0)
  return 2.0 * sqrt(Q2) / sqrt(rc * pow(L,3.0)) * exp(-arg)

def kolper_far_potential(kc, L, Q2, alpha):
  arg = pow(pi * kc / (alpha * L), 2.0)
  return 2.0 * sqrt(Q2) * sqrt(0.5 / pow(kc,3.0)) * alpha/(pi*pi) * exp(-arg)
#   return 2.0 * sqrt(Q2/N * (alpha*L/kc)**2.0 + Q2 / 2.0 / pow(kc,3.0)) * alpha/(pi*pi) * exp(-arg)
def kolper_far_field(kc, L, Q2, alpha):
  arg = pow(pi * kc / (alpha * L), 2.0)
  # estimate by Franzisk Nestler
#   return 2.0 * alpha / L * sqrt(Q2) / sqrt(pi * kc) * exp(-arg)
  # estimate by Kolafa/Perram
  return 2.0 * alpha / (pi * L) * sqrt(8 * Q2 / kc) * exp(-arg)

def kolper_near_energy(rc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_near_potential(rc, L, Q2, alpha)
def kolper_far_energy(kc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_far_potential(kc, L, Q2, alpha)
def kolper_near_force(rc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_near_field(rc, L, Q2, alpha)
def kolper_far_force(kc, L, Q2, alpha, N):
  return sqrt(Q2/N) * kolper_far_field(kc, L, Q2, alpha)

def kolper_near_total_energy(rc, L, Q2, alpha):
#   if rc < L:
#     const = sqrt(pi) / L**4.0 / alpha**3.0 * Q2 * exp(-(alpha*L)**2.0)
#   else:
#     const = sqrt(pi) / rc / L**3.0 / alpha**3.0 * Q2 * exp(-(alpha*rc)**2.0)
#   cor = Q2 * energy_correction_near(rc, L, alpha)
#   print("Kolafa-Perram NEAR correction = "+ str('%.2e' % const))
#   print("Numerical NEAR correction     = "+ str('%.2e' % cor))
  return sqrt(Q2/2.0) * kolper_near_potential(rc, L, Q2, alpha) / 2.0
def kolper_far_total_energy(kc, L, Q2, alpha):
#   arg = pow(pi * kc / (alpha * L), 2.0)
#   const = Q2 * alpha**2.0 * L / pi**2.0 / kc * exp(-arg)
#   cor   = Q2 * energy_correction_far(kc, L, alpha)
#   print("Kolafa-Perram FAR correction = "+ str('%.2e' % const))
#   print("Numerical FAR correction     = "+ str('%.2e' % cor))
#   
#   rf20 = energy_correction_far(kc, L, alpha) 
#   arg = pow(pi * kc / (alpha * L), 2.0)
#   rEf2 = alpha/(pi*pi) * pow(kc, -1.5) * exp(-arg)
# 
#   f20 = rf20**2
#   Ef2 = rEf2**2
# 
#   Ef2 = alpha**2 / (2.0 * kc**3 * pi**4) * exp(-2.0*arg)
# 
#   h1 = (Q2*Q2-Q4) * f20
#   h2 =  Q2*Q2 * Ef2 
# 
#   print("h1 = "+ str("%.2e" % h1) + ", h2 = "+ str("%.e" % h2))
#   U = sqrt( h1 + h2 )

  Unc  = sqrt(Q2/2.0) * kolper_far_potential(kc, L, Q2, alpha) / 2.0
#   Uold = cor + Unc
# 
#   print("total energy non-corrected  : "+ str("%.2e" % Unc) )
#   print("total energy corrected (old): "+ str("%.2e" % Uold) )
#   print("total energy corrected (new): "+ str("%.2e" % U) )
  return Unc 

def energy_correction_near(rc, L, alpha):
  R = 40
  cor = 0.0
  for k0 in range(-R,R):
    for k1 in range(-R,R):
      for k2 in range(-R,R):
        Lk = L * sqrt(k0*k0 + k1*k1 + k2*k2)
        if Lk > rc:
          cor = cor + erfc(alpha * Lk) / (Lk)
  return 0.5 * cor

def energy_correction_far(kc, L, alpha):
  R = 2 * int(ceil(kc))
  cor = 0.0
  for k0 in range(-R,R):
    for k1 in range(-R,R):
      for k2 in range(-R,R):
        k = k0*k0 + k1*k1 + k2*k2
        if k > kc**2:
#         if min(k0, k1, k2) < -kc or max(k0, k1, k2) >= kc:
          cor = cor + exp(-(pi/alpha/L)**2.0 * k) / k
  return 1.0 / (2.0*pi*L) * cor

def kc2M(kc):
  return int( 2 * ceil(kc) + 2 )

def err_tot(eps_near, eps_far):
  if split_tol == 0.5:
    return eps_near + eps_far
  else:
    return sqrt(eps_near*eps_near + eps_far*eps_far)

def print_errors(kc):
  print("\t\tnear\t\tfar\t\tboth\t\tnumerical")
  eps_near = kolper_near_potential(rc, L, Q2, alpha)
  eps_far  = kolper_far_potential(kc, L, Q2, alpha)
  print("abs Potential:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_potential))

  eps_near = kolper_near_energy(rc, L, Q2, alpha, N)
  eps_far  = kolper_far_energy(kc, L, Q2, alpha, N)
  print("abs Energy:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_energy))

  eps_near = kolper_near_field(rc, L, Q2, alpha)
  eps_far  = kolper_far_field(kc, L, Q2, alpha)
  print("abs Field:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_field))

  eps_near = kolper_near_force(rc, L, Q2, alpha, N)
  eps_far  = kolper_far_force(kc, L, Q2, alpha, N)
  print("abs Force:\t"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_force))

  eps_near = kolper_near_total_energy(rc, L, Q2, alpha)
  eps_far  = kolper_far_total_energy(kc, L, Q2, alpha)
  print("abs tot. Energy:"+ str('%.2e' % eps_near) +"\t"+ str('%.2e' % eps_far) +"\t"+ str('%.2e' % err_tot(eps_near, eps_far)) +"\t"+ str('%.2e' % abs_error_tenergy))

  print("\n\t\tnear\t\tfar\t\tboth\t\tnumerical")
  print("rel Potential:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_potential))
  print("rel Energy:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_energy))
  print("rel Field:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_field))
  print("rel Force:\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_force))
  print("rel tot. Energy:"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % 0) +"\t"+ str('%.2e' % rel_error_tenergy))


def unsigned(value):
  try:
    ivalue = int(value)
  except:
    raise argparse.ArgumentTypeError("invalid unsigned int value: %s" % value)
  if ivalue < 0:
    raise argparse.ArgumentTypeError("invalid unsigned int value: %s" % value)
  return ivalue

def posint(value):
  try:
    ivalue = int(value)
  except:
    raise argparse.ArgumentTypeError("invalid positive int value: %s" % value)
  if ivalue <= 0:
    raise argparse.ArgumentTypeError("invalid positive int value: %s" % value)
  return ivalue

def posfloat(value):
  try:
    fvalue = float(value)
  except:
    raise argparse.ArgumentTypeError("invalid positive float value: %s" % value)
  if fvalue <= 0:
    raise argparse.ArgumentTypeError("invalid positive float value: %s" % value)
  return fvalue


parser = argparse.ArgumentParser()

parser.add_argument("nd",
    help="number of periodic dimensions",
    choices=['0d', '1d', '2d', '3d'],
    type=str)
parser.add_argument("testcase",
    help="choose a test case",
    choices=['random', 'nacl', 'cloud-wall', 'cloud-ball', 'silica'],
    type=str)
parser.add_argument("tolerance",
    help="set tolerance",
    type=posfloat)

parser.add_argument("-s", "--testsize",
    help="all testcases support different number of particles (in most cases TESTSIZE is the multiplier of each box length)",
    default=1,
    type=posint)
parser.add_argument("-t", "--tolerance-type",
    help="set (absolute) tolerance type to potential (phi), energy (U), field (E), or force (F)",
    choices=['phi', 'U', 'E', 'F'],
    default='phi',
    type=str)
parser.add_argument("-n", "--num-procs",
    help="start parallel test runs with: mpirun -np NUM_PROCS",
    default=1,
    type=posint)

parser.add_argument("-r", "--rspace-cutoff",
    help="set real space cutoff",
    default=-1.0,
    type=float)
parser.add_argument("-d", "--direct", "--ndft",
    help="enable NDFT instead of NFFT",
    action="store_true")
parser.add_argument("-w", "--window",
    help="choose NFFT window function",
    choices=['gaussian', 'bspline', 'sinc', 'kaiser', 'bessel_i0'],
    default='bspline',
    type=str)
parser.add_argument("-m", "--window-cutoff",
    help="set NFFT window cutoff (ignored for ndft)",
    default=6,
    type=posint)
parser.add_argument("-o", "--oversampling",
    help="set oversampling '-2' (double: 2*M), '-1' (constant: M +2*m+2), integer 'value' >= 0 (add: M + value), float 'value' >= 0 (multiply: value*M)",
    default='1.0',
    type=str)
parser.add_argument("--oversampling-periodic", "--op",
    help="set periodic oversampling '-2' (double: 2*M), '-1' (constant: M +2*m+2), integer 'value' >= 0 (add: M + value), float 'value' >= 0 (multiply: value*M)",
    default='-3',
    type=str)
parser.add_argument("--oversampling-nonperiodic", "--on",
    help="set nonperiodic oversampling '-2' (double: 2*M), '-1' (constant: M +2*m+2), integer 'value' >= 0 (add: M + value), float 'value' >= 0 (multiply: value*M)",
    default='-3',
    type=str)
parser.add_argument("-p", "--reg-smoothness",
    help="set degree of smoothness for regularization (ignored for 3dp)",
#     metavar='p',
    type=unsigned)
parser.add_argument("-P", "--reg-gridpoints",
    help="set number of gridpoints for regularization (ignored for 3dp)",
    type=unsigned)
parser.add_argument("--kspace-ball",
    help="enable radial cutoff scheme in Fourier space",
    action="store_true")
parser.add_argument("--ignore-field",
    help="omit field computation",
    action="store_true")

args = parser.parse_args()

print args

nd = args.nd
tc = args.testcase
ts = args.testsize
tol = args.tolerance
tol_type = args.tolerance_type
m = args.window_cutoff
p = args.reg_smoothness
P = args.reg_gridpoints
window = args.window
osp = set_number(args.oversampling)
if args.oversampling_periodic != "-3":
  osp = set_number(args.oversampling_periodic)
osn = set_number(args.oversampling)
if args.oversampling_nonperiodic != "-3":
  osn = set_number(args.oversampling_nonperiodic)
nprocs = args.num_procs
rc = args.rspace_cutoff
kc_type = "ball" if args.kspace_ball else "box"




cwd = os.getcwd()

test_dir = "nonperiodic" if nd == "0d" else nd +"-periodic"

# switch test case
if tc=="cloud-wall":
  # cloud_wall
  scale = 2**(ts-1)
  N  = 300 * scale**3
  L  = 10.0 * scale
  rc = 4 if rc<=0 else rc
  Q2 = N
  Q4 = N
  tn = "cloud_wall"
  testname = "systems/"+ test_dir +"/cloud_wall_" + str(N) + ".xml.gz"
elif tc=="silica":
  # silica melt
  scale = 2**(ts-1)
  N  = 12960 * scale**3
  L  = 62.05966799605923 * scale
  rc = 8.0 if rc<=0 else rc
  Q2 = 3.73248000e+04 * scale**3
  Q4 = 1.61243136e+05 * scale**3
  tn = "silica_melt"
  testname = "systems/"+ test_dir +"/silica_melt_" + str(N) + ".xml.gz"
elif tc=="random":
  # random distribution
  N  = 11000 *(ts-1) + 1000
  L  = pow(N/1000.0, 1.0/3.0)
  rc = 0.62 if rc<=0 else rc
#   rc = 1.2
  Q2 = N
  Q4 = N
  tn = "random_dist"
  testname = "numerics/accuracy-tests-increased-box/rand_"+ nd +"p_" + str(N) + "_fmm.xml.gz"
  if not os.path.isfile(testname):
    print "ERROR: Test case is not part of ScaFaCoS! Did you create a proper softlink? If not execute"
    print "ERROR:   ln -s $PATH_TO_POTTS_AG_SVN/potts_ag/paper/2d-periodic-Ewald/numerics ."
    sys.exit()
elif tc=="cloud-ball":
  # inhomogeneous testcase of Prof. Winkler
  N = 50898
  L = 400
  rc = 32.379540 if rc<=0 else rc
  Q2 = 50.898075320917272
  Q4 = 0.050898150641900279
  tn = "cloud_ball"
  testname = "systems/"+ test_dir +"/cloud_ball_" + str(N) + ".xml.gz"
else:
  # nacl crystal
  N = 10648
  L = 22.0
  rc = 10.0 if rc<=0 else rc
  Q2 = N
  Q4 = N
  tn = "nacl"
  testname = "systems/"+ test_dir +"/nacl_" + str(N) + ".xml.gz"

print("Tune for tol = "+ str('%.2e' % tol))
print("L  = "+ str('%.2f' % L))
print("rc = "+ str('%.2f' % rc))
print("Q2 = %.8e" % Q2)
# print("Q4 = %.8e" % Q4)

alpha, kc, M, M0, M1, M2, m0, m1, m2, h, epsB = parameter_heuristic(
    rc, L, Q2, N, tol, tol_type, split_tol, P, m, osp, osn)

conf = " -c "
conf += "p2nfft_verbose_tuning,1,"
conf += "pfft_patience_name,estimate,"
conf += "p2nfft_ignore_tolerance,1,"
conf += "tolerance_field,1.000000e-20,"
conf += "p2nfft_r_cut,"+ str(rc) +","
conf += "p2nfft_alpha,"+ str(alpha) +","
conf += "pnfft_window_name,"+ window +","
conf += "pnfft_N,"+ str(M0) +","+ str(M1) +","+ str(M2) +","
conf += "pnfft_n,"+ str(m0) +","+ str(m1) +","+ str(m2) +","
conf += "pnfft_m,"+ str(m) +"," if m != None else ""
conf += "p2nfft_epsB,"+ str('%.4e' % epsB) +","
conf += "p2nfft_p,"+ str(p) +"," if p != None else ""
conf += "pnfft_direct,1," if args.direct else "pnfft_direct,0,"

intpol_order = "3"
conf += "p2nfft_intpol_order,"+ intpol_order +","
conf += "pnfft_intpol_order,"+ intpol_order +","
conf += "pnfft_grad_ik,0,"
conf += "p2nfft_k_cut,"+ str(kc) +"," if args.kspace_ball else ""
conf += "p2nfft_ignore_field,1," if args.ignore_field else ""

conf += "p2nfft_reg_kernel,0," if nd == "0d" else ""

nofield=""
# nofield+=" -u nofield"

startmpi = "mpirun -np "+ str(args.num_procs) +" " if args.num_procs > 1 else ""

line = startmpi + cwd + "/scafacos_test p2nfft " + testname + nofield + conf
print(line)

output = subprocess.check_output(line, stderr=subprocess.STDOUT, shell=True)
output = output.split('\n')
# output= ""

abs_error_potential = search_val(output, "abs_rms_potential_error", "=", "")
abs_error_energy    = search_val(output, "abs_rms_energy_error",    "=", "")
abs_error_field     = search_val(output, "abs_rms_field_error",     "=", "")
abs_error_force     = search_val(output, "abs_rms_force_error",     "=", "")
abs_error_tenergy   = search_val(output, "abs_total_energy_error",  "=", "")

rel_error_potential = search_val(output, "rel_rms_potential_error", "=", "")
rel_error_energy    = search_val(output, "rel_rms_energy_error",    "=", "")
rel_error_field     = search_val(output, "rel_rms_field_error",     "=", "")
rel_error_force     = search_val(output, "rel_rms_force_error",     "=", "")
rel_error_tenergy   = search_val(output, "rel_total_energy_error",  "=", "")

Q2 = search_val(output, "Q^2 =", "=", "")
Q4 = search_val(output, "Q^4 =", "=", "")
print("Test computed Q2 = %.8e" % Q2)
print("Test computed Q4 = %.8e" % Q4)

print("Error estimates with kc = "+ str('%.2f' % kc) +" predict following errors:")
print_errors(kc)

# print("Error estimates with kc = M/2 = "+ str('%.2f' % (M/2.0)) +" predict following errors:")
# print_errors(M/2.0)

print("Error estimates with kc = M/2-1 = "+ str('%.2f' % (M/2.0 - 1.0)) +" predict following errors:")
print_errors(M/2.0-1.0)

# print("Error estimates with kc = M/2-1 = "+ str('%.2e' % (M/2.0-1)) +" predict following errors:")
# print_errors(M/2.0-1.0)

if not args.direct:
  if is_in_text(output, "pnfft_trf_matrix_D("):
    rt_D, rt_F, rt_C, rt_G = filter_pnfft_runtimes(output)
    print("\nPNFFT runtimes:")
    print("D\t\tF\t\tC\t\tG")
    print(
        str('%.4e' % rt_D) +"\t"+
        str('%.4e' % rt_F) +"\t"+
        str('%.4e' % rt_C) +"\t"+
        str('%.4e' % rt_G)
        )

    rt_F_0, rt_F_1, rt_F_2, rt_F_twiddle, rt_F_remap, rt_F_3dto2d, rt_F_sum = filter_pfft_runtimes(output)
    headline = "trafo0\t\ttrafo1\t\ttrafo2\t\tremap\t\ttwiddle\t\t3dto2d\tsum"
    line = \
        str('%.4e' % rt_F_0) +"\t" \
      + str('%.4e' % rt_F_1) +"\t" \
      + str('%.4e' % rt_F_2) +"\t" \
      + str('%.4e' % rt_F_remap) +"\t" \
      + str('%.4e' % rt_F_twiddle) +"\t" \
      + str('%.4e' % rt_F_3dto2d) +"\t" \
      + str('%.4e' % rt_F_sum) +"\t"


    print("\nPFFT runtimes:")
    print(headline)
    print(line)

    time_file = "pfft_runtimes.txt"
    if not os.path.isfile(time_file):
      with open(time_file, 'a') as f:
        f.write("  M0  " + "  m0  " + headline + "\n")

    with open(time_file, 'a') as f:
      f.write(
          str('%4d  ' % M0)
        + str('%4d  ' % m0)
        + line + "\n") 


    # set up global variables for (scalar) objective function
    t0_near = search_val(output, "Near field computation takes", "takes", "s")
    t0_far_M = rt_F_twiddle + rt_F_3dto2d + rt_F_remap + rt_D + rt_C + rt_G
    t0_far_MlogM = rt_F_0 + rt_F_1 + rt_F_2
    rc_start = rc
    M0_start = M0

  #   t = estimate_time(t0_near, t0_far_M, t0_far_MlogM, rc, M0, 0.7)
  #   print("t = "+ str(t))

    res = minimize_scalar(estimate_time_onearg, bracket=(rc, 1.5*rc), method='brent', tol=1e-2)

    print("\noptimal rc = "+ str(res.x) + " gives estimated runtime of " + str('%.4e' % res.fun))

  else:
    print("\n!!! PNFFT runtimes not available, configure with --enable-fcs-timing=pnfft\n")

if is_in_text(output, "Near field computation takes"):
  rt_tune, rt_all, rt_near, rt_far = filter_p2nfft_runtimes(output)
  print("\nP2NFFT runtimes:")
  print("all\t\tnear\t\tfar\t\ttune")
  print(str('%.4e' % rt_all) +"\t"+  str('%.4e' % rt_near) +"\t"+ str('%.4e' % rt_far) +"\t"+ str('%.4e' % rt_tune))

  outfile="parameters_"+ tn +"_"+ nd +"p_"+ str('%.2e'%tol)
  outfile += "_ndft" if args.direct else "_nfft"
  outfile += "_kc" if args.kspace_ball else ""
  outfile += ".txt"
  if not os.path.isfile(outfile):
    with open(outfile, 'a') as f:
      f.write("# generated with:\n# \n" \
          "N          "+\
          "Q^2       "+\
          "L         "+\
          "rc        "+\
          "alpha     "+\
          "  M0  "+\
          "  M1  "+\
          "  M2  "+\
          "  m0  "+\
          "  m1  "+\
          "  m2  "+\
          "epsB      "+\
          "h         "+\
          " m  "+\
          "   P  "+\
          " p  "+\
          "abs_energy  "+\
          "abs_force   "+\
          "abs_tenergy  "+\
          "rel_energy  "+\
          "rel_force   "+\
          "rel_tenergy  "+\
          "rt-run    "+\
          "rt-near   "+\
          "rt-far    "+\
          "rt-tune   "+\
          "nfft  "+\
          "kc        "+\
          "nprocs  "+\
          "\n")

  with open(outfile, 'a') as f:
    txt  = str('%9d  ' % N)
    txt += str('%.2e  ' % Q2)
    txt += str('%.2e  ' % L)  
    txt += str('%f  ' % rc) 
    txt += str('%f  ' % alpha) 
    txt += str('%4d  ' % M0) 
    txt += str('%4d  ' % M1) 
    txt += str('%4d  ' % M2) 
    txt += str('%4d  ' % m0) 
    txt += str('%4d  ' % m1) 
    txt += str('%4d  ' % m2) 
    txt += str('%.2e  ' % epsB) 
    txt += str('%.2e  ' % h)    
    txt += str('%2d  ' % m) if m != None else "No"
    txt += str('%4d  ' % P) if P != None else "None" 
    txt += str('%4d  ' % p) if p != None else "None" 
    txt += str('%.2e    ' % abs_error_energy) 
    txt += str('%.2e    ' % abs_error_force)  
    txt += str('%.2e     ' % abs_error_tenergy)
    txt += str('%.2e    ' % rel_error_energy) 
    txt += str('%.2e    ' % rel_error_force)  
    txt += str('%.2e     ' % rel_error_tenergy)
    txt += str('%.2e  ' % rt_all)  
    txt += str('%.2e  ' % rt_near) 
    txt += str('%.2e  ' % rt_far)  
    txt += str('%.2e  ' % rt_tune) 
    txt += str(not args.direct) + "  " 
    txt += str('%.2e  ' % kc) if kc_type == "ball" else "box       "
    txt += str('%6d  ' % args.num_procs) 
    txt += "\n"
    f.write( txt.expandtabs(2) )
else:
  print("\n!!! P2NFFT runtimes not available, configure with --enable-fcs-timing\n")

