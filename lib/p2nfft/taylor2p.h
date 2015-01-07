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

#ifndef _P2NFFT_TAYLOR2P_H
#define _P2NFFT_TAYLOR2P_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "FCSDefinitions.h"
#include "kernels.h"

fcs_float ifcs_p2nfft_interpolate_symmetric(
    ifcs_p2nfft_kernel k, const fcs_float *param,
    fcs_int p, fcs_float x0, fcs_float x1, fcs_float x);
fcs_float ifcs_p2nfft_interpolate_ec(
    ifcs_p2nfft_kernel k, const fcs_float *param, fcs_float c,
    fcs_int p, fcs_float x0, fcs_float x1, fcs_float x);
fcs_float ifcs_p2nfft_interpolate_ic(
    ifcs_p2nfft_kernel k, const fcs_float *param,
    fcs_int p, fcs_float x0, fcs_float x1, fcs_float x);

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

