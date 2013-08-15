/*
  Copyright (C) 2012-2013 Michael Pippig

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.

  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

/* define some shortcuts in order to define setter, getter and set_tune at once */
#include "fcs_p2nfft_shortcuts.h"

/************************************************************
 *     Setter and Getter functions for p2nfft parameters
 ************************************************************/
/* setters/getters for tolerance */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_2(tolerance, tolerance, fcs_int, tolerance_type, fcs_float, tolerance_value)

/* Getters and Setters for near field cutoff radius */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(r_cut, r_cut, fcs_float, r_cut)

/* Getters and Setters for scaled near field cutoff radius */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(epsI, epsI, fcs_float, eps_I)

/* Getters and Setters for scaled far field regularization border */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(epsB, epsB, fcs_float, eps_I)

/* Getter and Setter for far field continuation value c used by taylor2p */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(c, c, fcs_float, c)

/* Getters and Setters for Ewald splitting parameter */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(alpha, alpha, fcs_float, alpha)

/* Getters and Setters for charge assignment order, this wrapper is compliant to P3M interface */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(cao, cao, fcs_int, cao)

/* Getters and Setters for order of P2NFFT near and far field interpolation */
FCS_P2NFFT_SET_GET_WRAPPER_1(interpolation_order, interpolation_order, fcs_int, intpol_order)

/* Getters and Setters for P2NFFT near field regularization flag */
FCS_P2NFFT_SET_GET_WRAPPER_1(reg_near, reg_near, fcs_int, reg)
FCS_P2NFFT_INTERFACE_WRAPPER_1(set_reg_near_by_name, set_reg_near_by_name, char*, reg_name)

/* Getters and Setters for P2NFFT far field regularization flag */
FCS_P2NFFT_SET_GET_WRAPPER_1(reg_far, reg_far, fcs_int, reg)
FCS_P2NFFT_INTERFACE_WRAPPER_1(set_reg_far_by_name, set_reg_far_by_name, char*, reg_name)

/* Getters and Setters for polynomial degree of near field regularization */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(p, p, fcs_int, p)

/* setters/getters for virial */
FCS_P2NFFT_INTERFACE_WRAPPER_1(require_virial,   require_virial,   fcs_int,    require_virial)
FCS_P2NFFT_INTERFACE_WRAPPER_1(get_virial,       get_virial,       fcs_float*, virial)
FCS_P2NFFT_INTERFACE_WRAPPER_1(virial_is_active, virial_is_active, fcs_int*,   yes_or_no)

/* Getters and Setters for ignore tolerance flag */
FCS_P2NFFT_SET_GET_WRAPPER_1(ignore_tolerance, ignore_tolerance, fcs_int, set_ignore_tolerance)

/************************************************************
 *     Setter and Getter functions for pnfft parameters
 ************************************************************/
/* Getters and Setters for FFT grid size, 2nd wrapper is compliant to P3M interface */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_3(pnfft_N,  pnfft_N, fcs_int,  N0, fcs_int,  N1, fcs_int,  N2)
#if FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_3(grid,     pnfft_N, fcs_int,  N0, fcs_int,  N1, fcs_int,  N2)
#endif

/* Getters and Setters for oversampled FFT grid size, 2nd wrapper is compliant to P3M interface */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_3(pnfft_n,          pnfft_n, fcs_int,  n0, fcs_int,  n1, fcs_int,  n2)
#if FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_3(oversampled_grid, pnfft_n, fcs_int,  n0, fcs_int,  n1, fcs_int,  n2)
#endif

/* Getters and Setters for PNFFT window function */
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_window, pnfft_window, fcs_int, window)
FCS_P2NFFT_INTERFACE_WRAPPER_1(set_pnfft_window_by_name, set_pnfft_window_by_name, char*, window_name)

/* Getters and Setters for charge assignment order */
FCS_P2NFFT_SET_GET_TUNE_WRAPPER_1(pnfft_m, pnfft_m, fcs_int, m)

/* Getters and Setters for PNFFT window interpolation order */
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_interpolation_order, pnfft_interpolation_order, fcs_int, intpol_order)

/* Getters and Setters for PNFFT flags */
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_pre_phi_hat,     pnfft_pre_phi_hat,     fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_fg_psi,          pnfft_fg_psi,          fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_fft_in_place,    pnfft_fft_in_place,    fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_sort_nodes,      pnfft_sort_nodes,      fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_interlaced,      pnfft_interlaced,      fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_grad_ik,         pnfft_grad_ik,         fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_pre_psi,         pnfft_pre_psi,         fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_pre_fg_psi,      pnfft_pre_fg_psi,      fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_pre_full_psi,    pnfft_pre_full_psi,    fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_pre_full_fg_psi, pnfft_pre_full_fg_psi, fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pnfft_real_f,          pnfft_real_f,          fcs_int, yes_or_no)


/************************************************************
 *     Setter and Getter functions for pfft parameters
 ************************************************************/

/* Getters and Setters for PFFT patience */
FCS_P2NFFT_SET_GET_WRAPPER_1(pfft_patience, pfft_patience, fcs_int, pfft_patience_flag)
FCS_P2NFFT_INTERFACE_WRAPPER_1(set_pfft_patience_by_name, set_pfft_patience_by_name, char*, patience_name)

/* Getters and Setters for PFFT flags */
FCS_P2NFFT_SET_GET_WRAPPER_1(pfft_preserve_input, pfft_preserve_input, fcs_int, yes_or_no)
FCS_P2NFFT_SET_GET_WRAPPER_1(pfft_tune,           pfft_tune,           fcs_int, yes_or_no)

