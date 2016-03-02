/*
  Copyright (C) 2011,2012 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/


#ifndef __EWALD_H__
#define __EWALD_H__


#define FCS_EWALD_USE_ERFC_APPROXIMATION 0

/* Mathematical constants, from gcc's math.h */
#ifndef M_PI
#define M_E             2.7182818284590452353602874713526625L  /* e */
#define M_LOG2E         1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10E        0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2           0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10          2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PI            3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2          1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4          0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PI          0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PI          0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPI      1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2         1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2       0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif

#define EWALD_DEFAULT_TOLERANCE_FIELD 1.0e-6

/* precision of alpha */
#ifdef FCS_FLOAT_IS_FLOAT
#define ALPHA_OPT_PREC 1.0e-5
#else
#define ALPHA_OPT_PREC 1.0e-10
#endif

/* Debug macros */
#ifdef FCS_ENABLE_DEBUG
#define ADDITIONAL_CHECKS
#define FCS_TRACE(cmd)      do { if (d->comm_rank == 0) { cmd; } } while (0)
#define FCS_DEBUG(cmd)      do { if (d->comm_rank == 0) { cmd; } } while (0)
#define FCS_DEBUG_ALL(cmd)  do { cmd; } while (0)
#else
#define FCS_TRACE(cmd)      do { } while (0)
#define FCS_DEBUG(cmd)      do { } while (0)
#define FCS_DEBUG_ALL(cmd)  do { } while (0)
#endif

#ifdef FCS_ENABLE_INFO
#define FCS_INFO(cmd)       do { if (d->comm_rank == 0) { cmd; } } while (0)
#else
#define FCS_INFO(cmd)       do { } while (0)
#endif

#define MAXKMAX_DEFAULT 100


/*static inline int on_root()
{
  static int rank = -1;
  if (rank < 0)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
}*/


inline static fcs_float SQR(fcs_float x) { return x*x; }


inline static fcs_int linindex(fcs_int ix, fcs_int iy, fcs_int iz, fcs_int kmax) {
  fcs_int km = kmax+1;
  return ix*km*km + iy*km + iz;
}


typedef struct {
  /* MPI communicator */
  MPI_Comm comm;
  int comm_size;
  int comm_rank;

  /* Cartesian communicator (for near field solver) */
  MPI_Comm comm_cart;

  /* System size in x,y,z */
  fcs_float box_l[3];

  /** Required accuracy. */
  fcs_float tolerance_field;

  /** Splitting parameter */
  fcs_float alpha;

  /** Near-field cutoff */
  fcs_float r_cut;

  /** Kspace cutoff */
  fcs_int kmax;

  /** maximal Kspace cutoff used by tuning */
  fcs_int maxkmax;

  /* influence function */
  fcs_float* G;

  /* Whether or not it is necessary to retune the method before running it. */
  int needs_retune;
  /** Whether or not rcut is to be automatically tuned. */
  int tune_r_cut;
  /** Whether or not kmax is to be automatically tuned. */
  int tune_kmax;
  /** Whether or not alpha is to be automatically tuned. */
  int tune_alpha;

  /* The components of the fields and potentials */
  fcs_float* far_fields;
  fcs_float* near_fields;
  fcs_float* far_potentials;
  fcs_float* near_potentials;
} ewald_data_struct;


void ewald_tune_alpha(fcs_int N, fcs_float sum_q2, fcs_float box_l[3], fcs_float r_cut, fcs_int kmax, fcs_float *alpha, fcs_float *error, int tune);
void ewald_compute_kspace(ewald_data_struct* d, fcs_int num_particles, fcs_float *positions, fcs_float *charges, fcs_float *fields, fcs_float *potentials);
void ewald_compute_rspace(ewald_data_struct* d, fcs_int num_particles, fcs_int max_num_particles, fcs_float *positions, fcs_float *charges, fcs_float *fields, fcs_float *potentials);


#endif /* __EWALD_H__ */
