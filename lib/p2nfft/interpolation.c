/*
 * Copyright (C) 2011-2014 Michael Pippig
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef HAVE_COMPLEX_H
#include <complex.h>
#endif

#include "interpolation.h"


fcs_float ifcs_p2nfft_interpolation_near(
    fcs_float x, fcs_float one_over_epsI,
    fcs_int order, fcs_int num_nodes,
    const fcs_float *table
    )
{
  switch(order){
    case 0: return ifcs_p2nfft_intpol_even_const(x, table, num_nodes, one_over_epsI);
    case 1: return ifcs_p2nfft_intpol_even_lin(x, table, num_nodes, one_over_epsI);
    case 2: return ifcs_p2nfft_intpol_even_quad(x, table, num_nodes, one_over_epsI);
    default: return ifcs_p2nfft_intpol_even_cub(x, table, num_nodes, one_over_epsI);
  }
}

fcs_float ifcs_p2nfft_interpolation_far(
    fcs_float x, fcs_float one_over_epsB,
    fcs_int order, fcs_int num_nodes,
    const fcs_float *table
    )
{
  switch(order){
    case 0: return ifcs_p2nfft_intpol_const(x, table, num_nodes, one_over_epsB);
    case 1: return ifcs_p2nfft_intpol_lin(x, table, num_nodes, one_over_epsB);
    case 2: return ifcs_p2nfft_intpol_quad(x, table, num_nodes, one_over_epsB);
    default: return ifcs_p2nfft_intpol_cub(x, table, num_nodes, one_over_epsB);
  }
}

