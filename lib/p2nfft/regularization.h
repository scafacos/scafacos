/*
 * Copyright (C) 2011-2013 Michael Pippig
 *
 * This file is part of ScaFaCoS.
 * 
 * ScaFaCoS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ScaFaCoS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *	
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _P2NFFT_REGULARIZATION_H
#define _P2NFFT_REGULARIZATION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "FCSDefinitions.h"
#include "kernels.h"

fcs_float ifcs_p2nfft_reg_far_rad_sym(
    ifcs_p2nfft_kernel k, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float epsI, fcs_float epsB);
fcs_float ifcs_p2nfft_reg_far_rad_expl_cont(
    ifcs_p2nfft_kernel k, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float epsI, fcs_float epsB, fcs_float c);
fcs_float ifcs_p2nfft_reg_far_rad_expl_cont_noncubic(
    ifcs_p2nfft_kernel k, fcs_float x2norm, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float r_cut, fcs_float eps_B, fcs_float c);
fcs_float ifcs_p2nfft_reg_far_rad_impl_cont(
    ifcs_p2nfft_kernel k, fcs_float xsnorm,
    fcs_int p, const fcs_float *param,
    fcs_float epsI, fcs_float epsB);

fcs_float ifcs_p2nfft_reg_far_rect_sym(
    fcs_float *x, fcs_float *h,
    fcs_int p, fcs_float epsB);
fcs_float ifcs_p2nfft_reg_far_rect_sym_version2(
    fcs_float *x, fcs_float *h,
    fcs_int p, fcs_float epsB);
fcs_float ifcs_p2nfft_reg_far_rect_expl_cont(
    fcs_float *x, fcs_float *h,
    fcs_int p, fcs_float epsB, fcs_float c);
fcs_float ifcs_p2nfft_reg_far_rect_impl_cont(
    fcs_float *x, fcs_float *h,
    fcs_int p, fcs_float epsB);

fcs_float ifcs_p2nfft_reg_far_no_singularity(
    ifcs_p2nfft_kernel k, fcs_float xx, fcs_int p,
    const fcs_float *param, fcs_float epsB);

fcs_float ifcs_regkern3_2ptaylor_x_inv(fcs_float *x, fcs_int p, fcs_float a, fcs_float b, fcs_float *n);

fcs_float ifcs_p2nfft_interpolation(
    fcs_float x, fcs_float one_over_epsI,
    fcs_int order, fcs_int num_nodes,
    const fcs_float *table);

fcs_float ifcs_p2nfft_nearfield_correction_taylor2p(
    fcs_float xx, fcs_int p, const fcs_float *param);
fcs_float ifcs_p2nfft_nearfield_correction_taylor2p_derive(
    fcs_float xx, fcs_int p, const fcs_float *param);

fcs_int ifcs_p2nfft_load_taylor2p_coefficients(
   fcs_int p,
   fcs_float *param);
fcs_int ifcs_p2nfft_load_taylor2p_derive_coefficients(
   fcs_int p,
   fcs_float *param);

#endif

