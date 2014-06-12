/*
 * Copyright (c) 2013 Benedikt Morbach
 *
 * This file is part of PNFFT.
 *
 * PNFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PNFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PNFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <complex.h>
#include "pnfft.h"
#include "ipnfft.h"

void PNX(debug_sum_print_strides)(
    R *data, INT max, int strides, int is_complex, const char *msg
    )
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  for(int c=0; c<strides; c++) {
    if(is_complex) {
      C sum = 0.0, gsum = 0.0;
      C *cdata = (C*)data;
  
      for(INT t=0; t<max; t++)
        sum += pnfft_fabs(pnfft_creal(cdata[strides*t+c])) + _Complex_I * pnfft_fabs(pnfft_cimag(cdata[strides*t+c]));
  
      MPI_Reduce(&sum, &gsum, 2, PNFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  
      if(!myrank) {
        fprintf(stderr, msg, strides);
        fprintf(stderr, ": %e + I* %e\n", pnfft_creal(gsum), pnfft_cimag(gsum));
      }
    } else {
      R sum = 0.0, gsum = 0.0;
  
      for(INT t=0; t<max; t++)
        sum += pnfft_fabs(data[strides*t+c]);
  
      MPI_Reduce(&sum, &gsum, 1, PNFFT_MPI_REAL_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  
      if(!myrank) {
        fprintf(stderr, msg, strides);
        fprintf(stderr, ": %e\n", gsum);
      }
    }
  }
}

void PNX(debug_sum_print)(
    R *data, INT max, int is_complex, const char *msg
    )
{
  PNX(debug_sum_print_strides)(data, max, 1, is_complex, msg);
}

