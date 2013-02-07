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

/** 24.02.2012: N = 8,16,32,64; eps_I = 1/8, 1/16, 1/32, 1/64 */
/** 24.02.2012: N = 128,256; eps_I = 1/8, 1/16, 1/32, 1/64 */
/** 24.02.2012: N = 8,...,256; eps_I = 1/128 */
/** 27.07.2012: N = 8,...,128; eps_I = 1/256, 1/512 */

/** 26.04.2012: N = 512,...,4096; eps_I = 1/64, 1/128, 1/256, 1/512 (no 3d error check) */
/** 27.07.2012: N = 512,...,4096; eps_I = 1/1024, 1/2048 (no 3d error check) */


#include "cg_cos_err.h"


fcs_float ifcs_p2nfft_get_cg_cos_err(
    fcs_int N, fcs_int log2_eps,
    fcs_int *m, fcs_int *p
    )
{
/** 24.02.2012: N = 8,16,32,64; eps_I = 1/8, 1/16, 1/32, 1/64 */
  /* N_exp=8, N_cos=4, p=16, eps_I=0.125000, eps_B=0.125000, err_1d=1.719732e-01, err_3d=3.996114e-02 */
  if((N == 4) && (log2_eps == 3)){
    *m = 2; *p = 16; return 3.996114e-02;
  }

  /* N_exp=16, N_cos=8, p=5, eps_I=0.125000, eps_B=0.125000, err_1d=7.428177e-03, err_3d=8.191734e-04 */
  if((N == 8) && (log2_eps == 3)){
    *m = 3; *p = 5; return 8.191734e-04;
  }

  /* N_exp=32, N_cos=16, p=6, eps_I=0.125000, eps_B=0.125000, err_1d=1.585061e-04, err_3d=1.005167e-04 */
  if((N == 16) && (log2_eps == 3)){
    *m = 4; *p = 6; return 1.005167e-04;
  }

  /* N_exp=64, N_cos=32, p=8, eps_I=0.125000, eps_B=0.125000, err_1d=3.277484e-07, err_3d=1.891323e-06 */
  if((N == 32) && (log2_eps == 3)){
    *m = 6; *p = 8; return 1.891323e-06;
  }

  /* N_exp=8, N_cos=4, p=11, eps_I=0.062500, eps_B=0.062500, err_1d=1.627212e+00, err_3d=1.438792e-01 */
  if((N == 4) && (log2_eps == 4)){
    *m = 2; *p = 11; return 1.438792e-01;
  }

  /* N_exp=16, N_cos=8, p=14, eps_I=0.062500, eps_B=0.062500, err_1d=2.330418e-01, err_3d=3.446089e-02 */
  if((N == 8) && (log2_eps == 4)){
    *m = 2; *p = 14; return 3.446089e-02;
  }

  /* N_exp=32, N_cos=16, p=8, eps_I=0.062500, eps_B=0.062500, err_1d=9.350762e-03, err_3d=2.913352e-04 */
  if((N == 16) && (log2_eps == 4)){
    *m = 4; *p = 8; return 2.913352e-04;
  }

  /* N_exp=64, N_cos=32, p=6, eps_I=0.062500, eps_B=0.062500, err_1d=1.174461e-04, err_3d=4.282328e-05 */
  if((N == 32) && (log2_eps == 4)){
    *m = 4; *p = 6; return 4.282328e-05;
  }

  /* N_exp=8, N_cos=4, p=8, eps_I=0.031250, eps_B=0.031250, err_1d=7.229135e+00, err_3d=1.708188e-01 */
  if((N == 4) && (log2_eps == 5)){
    *m = 2; *p = 8; return 1.708188e-01;
  }

  /* N_exp=16, N_cos=8, p=7, eps_I=0.031250, eps_B=0.031250, err_1d=2.553071e+00, err_3d=1.294960e-01 */
  if((N == 8) && (log2_eps == 5)){
    *m = 2; *p = 7; return 1.294960e-01;
  }

  /* N_exp=32, N_cos=16, p=9, eps_I=0.031250, eps_B=0.031250, err_1d=4.371732e-01, err_3d=3.131242e-02 */
  if((N == 16) && (log2_eps == 5)){
    *m = 3; *p = 9; return 3.131242e-02;
  }

  /* N_exp=64, N_cos=32, p=9, eps_I=0.031250, eps_B=0.031250, err_1d=1.537654e-02, err_3d=3.156110e-04 */
  if((N == 32) && (log2_eps == 5)){
    *m = 4; *p = 9; return 3.156110e-04;
  }

  /* N_exp=8, N_cos=4, p=6, eps_I=0.015625, eps_B=0.015625, err_1d=2.402003e+01, err_3d=1.370532e-01 */
  if((N == 4) && (log2_eps == 6)){
    *m = 2; *p = 6; return 1.370532e-01;
  }

  /* N_exp=16, N_cos=8, p=6, eps_I=0.015625, eps_B=0.015625, err_1d=1.309021e+01, err_3d=4.124942e-01 */
  if((N == 8) && (log2_eps == 6)){
    *m = 2; *p = 6; return 4.124942e-01;
  }

  /* N_exp=32, N_cos=16, p=2, eps_I=0.015625, eps_B=0.015625, err_1d=5.011635e+00, err_3d=6.648415e-02 */
  if((N == 16) && (log2_eps == 6)){
    *m = 2; *p = 2; return 6.648415e-02;
  }

  /* N_exp=64, N_cos=32, p=3, eps_I=0.015625, eps_B=0.015625, err_1d=8.308264e-01, err_3d=2.876863e-02 */
  if((N == 32) && (log2_eps == 6)){
    *m = 3; *p = 3; return 2.876863e-02;
  }

/** 24.02.2012: N = 128,256; eps_I = 1/8, 1/16, 1/32, 1/64 */

  /* N_exp=128, N_cos=64, p=15, eps_I=0.125000, eps_B=0.125000, err_1d=3.833787e-10, err_3d=7.530696e-10 */
  if((N == 64) && (log2_eps == 3)){
    *m = 6; *p = 15; return 7.530696e-10;
  }

  /* N_exp=256, N_cos=128, p=7, eps_I=0.125000, eps_B=0.125000, err_1d=1.296310e-09, err_3d=5.436340e-10 */
  if((N == 128) && (log2_eps == 3)){
    *m = 6; *p = 7; return 5.436340e-10;
  }

  /* N_exp=128, N_cos=64, p=8, eps_I=0.062500, eps_B=0.062500, err_1d=7.704587e-07, err_3d=4.767166e-07 */
  if((N == 64) && (log2_eps == 4)){
    *m = 6; *p = 8; return 4.767166e-07;
  }

  /* N_exp=256, N_cos=128, p=16, eps_I=0.062500, eps_B=0.062500, err_1d=9.239169e-10, err_3d=2.718172e-10 */
  if((N == 128) && (log2_eps == 4)){
    *m = 7; *p = 16; return 2.718172e-10;
  }

  /* N_exp=128, N_cos=64, p=6, eps_I=0.031250, eps_B=0.031250, err_1d=5.354156e-05, err_3d=2.180233e-05 */
  if((N == 64) && (log2_eps == 5)){
    *m = 4; *p = 6; return 2.180233e-05;
  }

  /* N_exp=256, N_cos=128, p=8, eps_I=0.031250, eps_B=0.031250, err_1d=4.425121e-07, err_3d=2.512948e-07 */
  if((N == 128) && (log2_eps == 5)){
    *m = 6; *p = 8; return 2.512948e-07;
  }

  /* N_exp=128, N_cos=64, p=16, eps_I=0.015625, eps_B=0.015625, err_1d=2.855394e-02, err_3d=8.541932e-04 */
  if((N == 64) && (log2_eps == 6)){
    *m = 4; *p = 16; return 8.541932e-04;
  }

  /* N_exp=256, N_cos=128, p=5, eps_I=0.015625, eps_B=0.015625, err_1d=3.989503e-05, err_3d=1.276327e-05 */
  if((N == 128) && (log2_eps == 6)){
    *m = 5; *p = 5; return 1.276327e-05;
  }

/** 24.02.2012: N = 8,...,256; eps_I = 1/128 */

  /* N_exp=8, N_cos=4, p=3, eps_I=0.007812, eps_B=0.007812, err_1d=5.060053e+01, err_3d=7.380878e-01 */
  if((N == 4) && (log2_eps == 7)){
    *m = 2; *p = 3; return 7.380878e-01;
  }

  /* N_exp=16, N_cos=8, p=11, eps_I=0.007812, eps_B=0.007812, err_1d=4.836725e+01, err_3d=1.582866e-01 */
  if((N == 8) && (log2_eps == 7)){
    *m = 2; *p = 11; return 1.582866e-01;
  }

  /* N_exp=32, N_cos=16, p=3, eps_I=0.007812, eps_B=0.007812, err_1d=2.347648e+01, err_3d=4.796375e-01 */
  if((N == 16) && (log2_eps == 7)){
    *m = 2; *p = 3; return 4.796375e-01;
  }

  /* N_exp=64, N_cos=32, p=10, eps_I=0.007812, eps_B=0.007812, err_1d=1.018307e+01, err_3d=5.488739e-02 */
  if((N == 32) && (log2_eps == 7)){
    *m = 3; *p = 10; return 5.488739e-02;
  }

  /* N_exp=128, N_cos=64, p=16, eps_I=0.007812, eps_B=0.007812, err_1d=1.428162e+00, err_3d=2.859214e-02 */
  if((N == 64) && (log2_eps == 7)){
    *m = 3; *p = 16; return 2.859214e-02;
  }

  /* N_exp=256, N_cos=128, p=16, eps_I=0.007812, eps_B=0.007812, err_1d=5.392272e-02, err_3d=6.143414e-04 */
  if((N == 128) && (log2_eps == 7)){
    *m = 4; *p = 16; return 6.143414e-04;
  }

/** 27.07.2012: N = 8,...,256; eps_I = 1/256, 1/512 */
  /* N_exp=8, N_cos=4, p=2, eps_I=0.003906, eps_B=0.003906, err_1d=4.591131e+01, err_3d=1.861977e+00 */
  if((N == 4) && (log2_eps == 8)){
    *m = 2; *p = 2; return 1.861977e+00;
  }

  /* N_exp=16, N_cos=8, p=2, eps_I=0.003906, eps_B=0.003906, err_1d=5.208141e+01, err_3d=2.034139e+00 */
  if((N == 8) && (log2_eps == 8)){
    *m = 2; *p = 2; return 2.034139e+00;
  }

  /* N_exp=32, N_cos=16, p=9, eps_I=0.003906, eps_B=0.003906, err_1d=5.055835e+01, err_3d=1.420772e+00 */
  if((N == 16) && (log2_eps == 8)){
    *m = 2; *p = 9; return 1.420772e+00;
  }

  /* N_exp=64, N_cos=32, p=9, eps_I=0.003906, eps_B=0.003906, err_1d=3.390603e+01, err_3d=3.865937e-01 */
  if((N == 32) && (log2_eps == 8)){
    *m = 2; *p = 9; return 3.865937e-01;
  }

  /* N_exp=128, N_cos=64, p=3, eps_I=0.003906, eps_B=0.003906, err_1d=2.008646e+01, err_3d=1.368964e-01 */
  if((N == 64) && (log2_eps == 8)){
    *m = 3; *p = 3; return 1.368964e-01;
  }

  /* N_exp=256, N_cos=128, p=16, eps_I=0.003906, eps_B=0.003906, err_1d=3.412298e+00, err_3d=3.775531e-02 */
  if((N == 128) && (log2_eps == 8)){
    *m = 3; *p = 16; return 3.775531e-02;
  }

  /* N_exp=32, N_cos=16, p=2, eps_I=0.001953, eps_B=0.001953, err_1d=1.091123e+02, err_3d=2.777744e+00 */
  if((N == 16) && (log2_eps == 9)){
    *m = 2; *p = 2; return 2.777744e+00;
  }

  /* N_exp=64, N_cos=32, p=15, eps_I=0.001953, eps_B=0.001953, err_1d=8.615697e+01, err_3d=5.884714e-01 */
  if((N == 32) && (log2_eps == 9)){
    *m = 2; *p = 15; return 5.884714e-01;
  }

  /* N_exp=128, N_cos=64, p=11, eps_I=0.001953, eps_B=0.001953, err_1d=6.579775e+01, err_3d=8.598914e-01 */
  if((N == 64) && (log2_eps == 9)){
    *m = 2; *p = 11; return 8.598914e-01;
  }

  /* N_exp=256, N_cos=128, p=12, eps_I=0.001953, eps_B=0.001953, err_1d=3.609815e+01, err_3d=1.567691e-01 */
  if((N == 128) && (log2_eps == 9)){
    *m = 3; *p = 12; return 1.567691e-01;
  }

/** 26.04.2012: N = 512,...,4096; eps_I = 1/64, 1/128, 1/256, 1/512 (no 3d error check) */

  /* N_exp=512, N_cos=256, p=7, eps_I=0.015625, eps_B=0.015625, err_1d=9.703791e-08, err_3d=-1.000000e+00 */
  if((N == 256) && (log2_eps == 6)){
    *m = -1; *p = 7; return -1.000000e+00;
  }

  /* N_exp=1024, N_cos=512, p=15, eps_I=0.015625, eps_B=0.015625, err_1d=1.260005e-10, err_3d=-1.000000e+00 */
  if((N == 512) && (log2_eps == 6)){
    *m = -1; *p = 15; return -1.000000e+00;
  }

  /* N_exp=2048, N_cos=1024, p=11, eps_I=0.015625, eps_B=0.015625, err_1d=7.321432e-11, err_3d=-1.000000e+00 */
  if((N == 1024) && (log2_eps == 6)){
    *m = -1; *p = 11; return -1.000000e+00;
  }

  /* N_exp=4096, N_cos=2048, p=11, eps_I=0.015625, eps_B=0.015625, err_1d=1.895728e-11, err_3d=-1.000000e+00 */
  if((N == 2048) && (log2_eps == 6)){
    *m = -1; *p = 11; return -1.000000e+00;
  }

  /* N_exp=512, N_cos=256, p=8, eps_I=0.007812, eps_B=0.007812, err_1d=7.144070e-05, err_3d=-1.000000e+00 */
  if((N == 256) && (log2_eps == 7)){
    *m = -1; *p = 8; return -1.000000e+00;
  }

  /* N_exp=1024, N_cos=512, p=7, eps_I=0.007812, eps_B=0.007812, err_1d=4.788143e-08, err_3d=-1.000000e+00 */
  if((N == 512) && (log2_eps == 7)){
    *m = -1; *p = 7; return -1.000000e+00;
  }

  /* N_exp=2048, N_cos=1024, p=14, eps_I=0.007812, eps_B=0.007812, err_1d=3.452971e-11, err_3d=-1.000000e+00 */
  if((N == 1024) && (log2_eps == 7)){
    *m = -1; *p = 14; return -1.000000e+00;
  }

  /* N_exp=4096, N_cos=2048, p=12, eps_I=0.007812, eps_B=0.007812, err_1d=1.224976e-11, err_3d=-1.000000e+00 */
  if((N == 2048) && (log2_eps == 7)){
    *m = -1; *p = 12; return -1.000000e+00;
  }

  /* N_exp=512, N_cos=256, p=16, eps_I=0.003906, eps_B=0.003906, err_1d=8.970433e-02, err_3d=-1.000000e+00 */
  if((N == 256) && (log2_eps == 8)){
    *m = -1; *p = 16; return -1.000000e+00;
  }

  /* N_exp=1024, N_cos=512, p=7, eps_I=0.003906, eps_B=0.003906, err_1d=1.444791e-04, err_3d=-1.000000e+00 */
  if((N == 512) && (log2_eps == 8)){
    *m = -1; *p = 7; return -1.000000e+00;
  }

  /* N_exp=2048, N_cos=1024, p=7, eps_I=0.003906, eps_B=0.003906, err_1d=2.352167e-08, err_3d=-1.000000e+00 */
  if((N == 1024) && (log2_eps == 8)){
    *m = -1; *p = 7; return -1.000000e+00;
  }

  /* N_exp=4096, N_cos=2048, p=14, eps_I=0.003906, eps_B=0.003906, err_1d=1.844569e-11, err_3d=-1.000000e+00 */
  if((N == 2048) && (log2_eps == 8)){
    *m = -1; *p = 14; return -1.000000e+00;
  }

  /* N_exp=512, N_cos=256, p=4, eps_I=0.001953, eps_B=0.001953, err_1d=5.337899e+00, err_3d=-1.000000e+00 */
  if((N == 256) && (log2_eps == 9)){
    *m = -1; *p = 4; return -1.000000e+00;
  }

  /* N_exp=1024, N_cos=512, p=5, eps_I=0.001953, eps_B=0.001953, err_1d=1.725060e-01, err_3d=-1.000000e+00 */
  if((N == 512) && (log2_eps == 9)){
    *m = -1; *p = 5; return -1.000000e+00;
  }

  /* N_exp=2048, N_cos=1024, p=6, eps_I=0.001953, eps_B=0.001953, err_1d=2.745911e-04, err_3d=-1.000000e+00 */
  if((N == 1024) && (log2_eps == 9)){
    *m = -1; *p = 6; return -1.000000e+00;
  }

  /* N_exp=4096, N_cos=2048, p=7, eps_I=0.001953, eps_B=0.001953, err_1d=1.115607e-08, err_3d=-1.000000e+00 */
  if((N == 2048) && (log2_eps == 9)){
    *m = -1; *p = 7; return -1.000000e+00;
  }

/** 27.07.2012: N = 512,...,4096; eps_I = 1/1024, 1/2048 (no 3d error check) */

  /* N_exp=512, N_cos=256, p=9, eps_I=0.000977, eps_B=0.000977, err_1d=7.129043e+01, err_3d=-1.000000e+00 */
  if((N == 256) && (log2_eps == 10)){
    *m = -1; *p = 9; return -1.000000e+00;
  }

  /* N_exp=1024, N_cos=512, p=7, eps_I=0.000977, eps_B=0.000977, err_1d=1.081330e+01, err_3d=-1.000000e+00 */
  if((N == 512) && (log2_eps == 10)){
    *m = -1; *p = 7; return -1.000000e+00;
  }

  /* N_exp=2048, N_cos=1024, p=3, eps_I=0.000977, eps_B=0.000977, err_1d=3.530704e-01, err_3d=-1.000000e+00 */
  if((N == 1024) && (log2_eps == 10)){
    *m = -1; *p = 3; return -1.000000e+00;
  }

  /* N_exp=4096, N_cos=2048, p=5, eps_I=0.000977, eps_B=0.000977, err_1d=5.101622e-04, err_3d=-1.000000e+00 */
  if((N == 2048) && (log2_eps == 10)){
    *m = -1; *p = 5; return -1.000000e+00;
  }

  /* N_exp=512, N_cos=256, p=15, eps_I=0.000488, eps_B=0.000488, err_1d=2.071658e+02, err_3d=-1.000000e+00 */
  if((N == 256) && (log2_eps == 11)){
    *m = -1; *p = 15; return -1.000000e+00;
  }

  /* N_exp=1024, N_cos=512, p=6, eps_I=0.000488, eps_B=0.000488, err_1d=1.336593e+02, err_3d=-1.000000e+00 */
  if((N == 512) && (log2_eps == 11)){
    *m = -1; *p = 6; return -1.000000e+00;
  }

  /* N_exp=2048, N_cos=1024, p=7, eps_I=0.000488, eps_B=0.000488, err_1d=2.067662e+01, err_3d=-1.000000e+00 */
  if((N == 1024) && (log2_eps == 11)){
    *m = -1; *p = 7; return -1.000000e+00;
  }

  /* N_exp=4096, N_cos=2048, p=11, eps_I=0.000488, eps_B=0.000488, err_1d=6.886330e-01, err_3d=-1.000000e+00 */
  if((N == 2048) && (log2_eps == 11)){
    *m = -1; *p = 11; return -1.000000e+00;
  }




  /* Configuration not found */
  return -1.0;
}

