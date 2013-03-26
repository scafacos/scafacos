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

static unsigned extract_pfft_opt_flags(
    unsigned pfft_flags);
static void fft_output_size(
    const INT *n, const R *x_max, int m,
    INT *no);


void PNX(local_size_guru)(
    int d, const INT *N, const INT *n, const R *x_max, int m,
    MPI_Comm comm_cart, unsigned pnfft_flags,
    INT *local_N, INT *local_N_start,
    R *lower_border, R *upper_border
    )
{
  INT no[3], local_no[3], local_no_start[3];

  if(d != 3){
    PX(fprintf)(comm_cart, stderr, "!!! Error in PNFFT: d != 3 not yet implemented !!!\n");
    return;
  }

  fft_output_size(n, x_max, m,
      no);

  PNX(local_size_internal)(N, n, no, comm_cart, pnfft_flags,
      local_N, local_N_start, local_no, local_no_start);

  PNX(node_borders)(n, local_no, local_no_start, x_max,
      lower_border, upper_border);
}


PNX(plan) PNX(init_guru)(
    int d, const INT *N, const INT *n, const R *x_max,
    INT local_M, int m,
    unsigned pnfft_flags, unsigned pfft_flags,
    MPI_Comm comm_cart
    )
{
  INT no[3];
  PNX(plan) ths;
  unsigned pfft_opt_flags = extract_pfft_opt_flags(pfft_flags);
  
  if(d != 3){
    PX(fprintf)(comm_cart, stderr, "!!! Error in PNFFT: d != 3 not yet implemented !!!\n");
    return NULL;
  }
  
  fft_output_size(n, x_max, m,
    no);

#if PNFFT_DEBUG_USE_KAISER_BESSEL | PNFFT_DEBUG_USE_GAUSSIAN | PNFFT_DEBUG_USE_BSPLINE | PNFFT_DEBUG_USE_SINC_POWER
  int rank=0;
  MPI_Comm_rank(comm_cart, &rank);
#endif

#if PNFFT_DEBUG_USE_KAISER_BESSEL
  if(!rank) fprintf(stderr, "!!! Debugging: Force PNFFT_WINDOW_KAISER_BESSEL !!!\n");
  pnfft_flags &= ~(PNFFT_WINDOW_GAUSSIAN | PNFFT_WINDOW_BSPLINE | PNFFT_WINDOW_SINC_POWER);
#elif PNFFT_DEBUG_USE_GAUSSIAN
  if(!rank) fprintf(stderr, "!!! Debugging: Force PNFFT_WINDOW_GAUSSIAN !!!\n");
  pnfft_flags &= ~(PNFFT_WINDOW_GAUSSIAN | PNFFT_WINDOW_BSPLINE | PNFFT_WINDOW_SINC_POWER);
  pnfft_flags |= PNFFT_WINDOW_GAUSSIAN;
#elif PNFFT_DEBUG_USE_BSPLINE
  if(!rank) fprintf(stderr, "!!! Debugging: Force PNFFT_WINDOW_BSPLINE !!!\n");
  pnfft_flags &= ~(PNFFT_WINDOW_GAUSSIAN | PNFFT_WINDOW_BSPLINE | PNFFT_WINDOW_SINC_POWER);
  pnfft_flags |= PNFFT_WINDOW_BSPLINE;
#elif PNFFT_DEBUG_USE_SINC_POWER
  if(!rank) fprintf(stderr, "!!! Debugging: Force PNFFT_WINDOW_SINC_POWER !!!\n");
  pnfft_flags &= ~(PNFFT_WINDOW_GAUSSIAN | PNFFT_WINDOW_BSPLINE | PNFFT_WINDOW_SINC_POWER);
  pnfft_flags |= PNFFT_WINDOW_SINC_POWER;
#endif

  ths = PNX(init_internal)(d, N, n, no, local_M, m, pnfft_flags, pfft_opt_flags, comm_cart);

  /* Quick fix to save x_max in PNFFT plan */
  for(int t=0; t<d; t++)
    ths->x_max[t] = x_max[t];

  return ths;
}

static unsigned extract_pfft_opt_flags(
    unsigned pfft_flags
    )
{
  /* PFFT default: PFFT_MEASURE| PFFT_NO_TUNE */
  unsigned flags = 0;

  if(pfft_flags & PFFT_ESTIMATE)
    flags = PFFT_ESTIMATE;
  if(pfft_flags & PFFT_PATIENT)
    flags = PFFT_PATIENT;
  if(pfft_flags & PFFT_EXHAUSTIVE)
    flags = PFFT_EXHAUSTIVE;

  if(pfft_flags & PFFT_TUNE)
    flags |= PFFT_TUNE;
  if(pfft_flags & PFFT_PRESERVE_INPUT)
    flags |= PFFT_PRESERVE_INPUT;
  if(pfft_flags & PFFT_DESTROY_INPUT)
    flags |= PFFT_DESTROY_INPUT;

  return flags;
}




static void fft_output_size(
    const INT *n, const R *x_max, int m,
    INT *no
    )
{
  INT c;

  /* For the default case x_max == 0.5, we do not allow nodes with x = +0.5,
   * since they can be mapped to -0.5.
   * Things change, if we use a smaller x_max. Then we need to save a further grid point to enable
   * the calculation of nodes with x == +x_max. This explains the +1 in calculation of no. */
  /* Calculation of no:
   * Project the right most node into the grid: cr = floor(n*x_max) (lrint converts this into integer).
   * Project the left most node into the grid: cl = floor(-n*x_max) (lrint converts this into integer).
   * From these indices, we sum up the next 'm' grid points to the left and the next 'm+1' grid points to the right.
   * Therefore we need all coefficients g_l with cl - m <= l <= cr + m+1.
   * This is also fullfiled, if we choose -cr-m-1 <= l <= cr+m+1 (since we have either cl=-cr or cl=-cr-1).
   * Since FFT gives back g_l with -no/2 <= l <= no/2-1, we choose no = 2*(cr+m+1 +1) = 2*(cr+m+2)  */
  
  for(int t=0; t<3; t++){
    c = pnfft_lrint(pnfft_floor(n[t]*x_max[t]));
    no[t] = PNFFT_MIN(n[t], 2*(c+m+2));
  }
}


