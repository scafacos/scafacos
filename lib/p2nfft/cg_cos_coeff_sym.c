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

#include "cg_cos_coeff_sym.h"


fcs_int ifcs_p2nfft_load_cg_cos_coeff_sym(
    fcs_int N, fcs_int log2_eps,
    fcs_int *m, fcs_int *p, fcs_float *f_hat
    )
{
  /* N_exp=16, N_cos=8, p=16, eps_I=0.125000, eps_B=0.125000, err_1d=3.408852e-03, err_3d=-1.000000e+00 */
  if((N == 8) && (log2_eps == 3)){
    *m = -1; *p = 16;
    f_hat[0] = 5.9794406382675973e+00;    f_hat[1] = 5.4394969049218114e+00;    f_hat[2] = 2.5252186038430637e+00;    f_hat[3] = 1.3954098329782483e+00;
    f_hat[4] = 6.4854991600274414e-01;    f_hat[5] = 3.1167550276815259e-01;    f_hat[6] = 1.0561198531900384e-01;    f_hat[7] = 3.7601840103756001e-02;
    return 0;
  }

  /* N_exp=32, N_cos=16, p=6, eps_I=0.125000, eps_B=0.125000, err_1d=4.162614e-06, err_3d=-1.000000e+00 */
  if((N == 16) && (log2_eps == 3)){
    *m = 4; *p = 6;
    f_hat[0] = 6.6363313034308007e+00;    f_hat[1] = 6.7015017957535532e+00;    f_hat[2] = 3.6961006927906621e+00;    f_hat[3] = 2.3945683485233360e+00;
    f_hat[4] = 1.4711679594787628e+00;    f_hat[5] = 9.2772861106985027e-01;    f_hat[6] = 5.5006897143005318e-01;    f_hat[7] = 3.2211274071481200e-01;
    f_hat[8] = 1.7603504753832683e-01;    f_hat[9] = 9.2167504589197272e-02;    f_hat[10] = 4.4501557780113228e-02;    f_hat[11] = 1.9734983628827810e-02;
    f_hat[12] = 7.8669799979926135e-03;    f_hat[13] = 2.6462246429786178e-03;    f_hat[14] = 7.3031189488311345e-04;    f_hat[15] = 1.2912420373359756e-04;
    return 0;
  }

  /* Configuration not found */
  return 1;
}

