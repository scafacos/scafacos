/*
 * Copyright (c) 2011-2013 Michael Pippig
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

#include "pnfft.h"
#include "ipnfft.h"

static int local_nodes_out_of_range(
    const R *lo, const R *up, const R *x_max);
static R random_number_less_than_one(
    void);
static void default_fft_size(
    const INT *N,
    INT *n, R *x_max);
static int default_m(
    void);


void PNX(init_x_3d_adv)(
    const R *lo, const R *up,
    const R *x_max, INT loc_M,
    R *x
    )
{
  R tmp;
  
  if(local_nodes_out_of_range(lo, up, x_max))
    return;
  
  for (INT j=0; j<loc_M; j++){
    for(int t=0; t<3; t++){
      do{
        tmp = random_number_less_than_one();
        tmp = (up[t]-lo[t]) * tmp + lo[t];
      }
      while( (-x_max[t] < tmp) || (tmp <= x_max[t]) );
      x[3*j+t] = tmp;
    }
  }
}


static int local_nodes_out_of_range(
    const R *lo, const R *up, const R *x_max
    )
{
  int rtnval = 0;
  
  for(int t=0; t<3; t++)
    if( (lo[t] < -x_max[t]) || (x_max[t] < up[t]) )
      rtnval = 1;
  
  return rtnval;
}


static R random_number_less_than_one(
    void
    )
{
  int itmp;
  R tmp;
  
  do {
    itmp = rand();
    tmp = (R) itmp/RAND_MAX;
  } while(tmp - 1.0 < PNFFT_EPSILON);
  
  return tmp;
}
    
    

void PNX(local_size_adv)(
    int d, const INT *N, MPI_Comm comm_cart,
    unsigned pnfft_flags,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    )
{
  int m;
  INT n[3];
  R x_max[3];

  m = default_m();
  default_fft_size(N,
      n, x_max);

  PNX(local_size_guru)(
      d, N, n, x_max, m, comm_cart, pnfft_flags,
      local_N, local_N_start, lower_border, upper_border);
}


PNX(plan) PNX(init_adv)(
    int d, const INT *N,
    INT local_M,
    unsigned pnfft_flags, unsigned pfft_flags,
    MPI_Comm comm_cart
    )
{
  int m;
  INT n[3];
  R x_max[3];

  m = default_m();
  default_fft_size(N,
      n, x_max);

  return PNX(init_guru)(
      d, N, n, x_max, local_M, m,
      pnfft_flags, pfft_flags, comm_cart);
}


static void default_fft_size(
    const INT *N,
    INT *n, R *x_max
    )
{
  for(int t=0; t<3; t++){
    n[t] = 2*N[t];
    x_max[t] = 0.5;
  }
}


/* We choose the same default real space cutoff 'm' for all
 * window functions. Therefore, we do not need the pnfft_flags for 
 * local_size_adv */
static int default_m(
    void
    )
{
  return 6;
}
