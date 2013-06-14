/*
 * Copyright (C) 2013 Michael Pippig
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

#ifndef _P2NFFT_BESSEL_K_H_
#define _P2NFFT_BESSEL_K_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "FCSDefinitions.h"

fcs_float ifcs_p2nfft_bessel_k(
    fcs_float nu, fcs_float x);
fcs_float ifcs_p2nfft_inc_lower_bessel_k(
    fcs_float nu, fcs_float x, fcs_float y, fcs_float eps);
fcs_float ifcs_p2nfft_inc_upper_bessel_k(
    fcs_float nu, fcs_float x, fcs_float y, fcs_float eps);
void ifcs_p2nfft_plot_slavinsky_safouhi_table(
    void);

#endif
