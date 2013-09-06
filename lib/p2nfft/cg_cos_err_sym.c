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


#include "cg_cos_err_sym.h"


fcs_float ifcs_p2nfft_get_cg_cos_err_sym(
    fcs_int N, fcs_int log2_eps,
    fcs_int *m, fcs_int *p
    )
{
  /* N_exp=32, N_cos=16, p=6, eps_I=0.125000, eps_B=0.125000, err_1d=4.162614e-06, err_3d=-1.000000e+00 */
  if((N == 16) && (log2_eps == 3)){
    *m = 5; *p = 6; return -1.000000e+00;
  }


  /* Configuration not found */
  return -1.0;
}

