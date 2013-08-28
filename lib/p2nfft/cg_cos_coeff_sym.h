/*
 * Copyright (C) 2012,2013 Michael Pippig
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

#ifndef _P2NFFT_CG_COS_COEFF_SYM_H_
#define _P2NFFT_CG_COS_COEFF_SYM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"

fcs_int ifcs_p2nfft_load_cg_cos_coeff_sym(
    fcs_int N, fcs_int log2_eps,
    fcs_int *m, fcs_int *p, fcs_float *f_hat);

#endif


