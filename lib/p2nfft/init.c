/*
 * Copyright (C) 2011-2013 Michael Pippig
 * Copyright (C) 2011 Sebastian Banert
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
#include "init.h"
#include "types.h"
#include "utils.h"

/* P2NFFT naturally needs a two-dimensional procmesh.
 * If we use a three-dimensional procmesh, PNFFT performs an expensive 3d to 2d remap.
 * Nevertheless, the near field module needs a three-dimensional procmesh.
 * Therefore, create a three-dimensional procmesh with last dimension equal to 1.*/
#define FCS_P2NFFT_USE_3D_PROCMESH 1 

static void calc_grid_sizes(
    fcs_int nprocs, fcs_int ndim, int *np);
static ifcs_p2nfft_data_struct* mkplan_p2nfft(
    void);
static int comm_is_cart_3d(
    MPI_Comm comm);
static void comm_get_periodicity(
    MPI_Comm comm, fcs_int *periodicity);

#if FCS_P2NFFT_USE_3D_PROCMESH 
static void comm_create_cart_3d(
    MPI_Comm comm, MPI_Comm *comm_cart_3d, int *np);
#else
static void comm_create_cart_2d(
    MPI_Comm comm, MPI_Comm *comm_cart_2d, MPI_Comm *comm_cart_3d, int *np);
#endif



FCSResult ifcs_p2nfft_init(
    void **rd, MPI_Comm comm
    )
{
  const char *fnc_name = "ifcs_p2nfft_init";
  ifcs_p2nfft_data_struct *d;

  /* return error if method context is already allocated */
  if (*rd != NULL) 
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, "Multiple init of method context without finalize.");
  
  /* Initialize the PNFFT library */
  FCS_PNFFT(init)();
 
  /* Create data structure */ 
  d = mkplan_p2nfft();

  /* return error if allocation failed */
  if (d == NULL)
    return fcsResult_create(FCS_ALLOC_FAILED, fnc_name, "Allocation of the p2nfft data structure failed.");

#if FCS_P2NFFT_USE_3D_PROCMESH 
  /* Create a three-dimensional cartesian comm
     from the given (possibly non-cartesian) one. */
  if( !comm_is_cart_3d(comm) )
    comm_create_cart_3d(comm, &d->cart_comm_3d, d->np);
  else {
    int periods[3], coords[3];
    MPI_Cart_get(comm, 3,
        d->np, periods, coords);

    if( periods[0] && periods[1] && periods[2] )
      MPI_Comm_dup(comm, &d->cart_comm_3d);
    else {
      for(int t=0; t<3; t++)
        periods[t] = 1;
      MPI_Cart_create(comm, 3, d->np, periods, 0, &d->cart_comm_3d);
    }
  }
  MPI_Comm_dup(d->cart_comm_3d, &d->cart_comm_pnfft);
#else  
  /* create 2d cart procmesh for PNFFT and its 3d counterpart */
  comm_create_cart_2d(comm, &d->cart_comm_pnfft, &d->cart_comm_3d, d->np);
#endif

  /* Set the default values */
  d->needs_retune = 1;
  d->tune_alpha = 1;
  d->tune_r_cut = 1;
  d->tune_epsI = 1;
  d->tune_epsB = 1;
  d->tune_N = 1;
  d->tune_n = 1;
  d->tune_m = 1;
  d->tune_p = 1;
  d->tune_c = 1;
  d->flags = FCS_P2NFFT_CHECK_TOLERANCE; /* 1: continue even if accuracy estimation fails */

  d->pnfft_flags = PNFFT_MALLOC_F_HAT| PNFFT_PRE_PHI_HAT | PNFFT_FFT_OUT_OF_PLACE | PNFFT_TRANSPOSED_F_HAT;
  d->pnfft_interpolation_order = 3;
  d->pnfft_window = FCS_P2NFFT_DEFAULT_PNFFT_WINDOW;
  d->pfft_flags = PFFT_NO_TUNE | PFFT_DESTROY_INPUT;
  d->pfft_patience = FCS_P2NFFT_DEFAULT_PFFT_PATIENCE;

  /* We do not know the default tolerance type at this point, since periodicity
   * may be changed via fcs_set_periodicity after fcs_init. */
  d->tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
  d->tolerance = -1.0;
  
  d->N[0] = d->N[1] = d->N[2] = 16;
  d->m = 4;
  d->p = 8;
  d->c = 0.0;

  /* init to same nonsense on all processes */
  d->alpha = -1.0;
  d->r_cut = -1.0;
  d->one_over_r_cut = -1.0;
  d->epsI = -1.0;
  d->epsB = -1.0;
  d->num_nodes = -1;
  d->sum_qpart = -1;
  d->sum_q2 = -1.0;
  d->sum_q = 0.0;
  d->bg_charge = 0.0;
  for(int t=0; t<3; t++){
    d->box_l[t] = -1.0; 
    d->box_scales[t] = 1.0;
    d->box_shifts[t] = 0.0;
  }
  
  comm_get_periodicity(comm, d->periodicity);

  d->short_range_flag = -1;
  d->reg_near = FCS_P2NFFT_REG_NEAR_DEFAULT;
  d->reg_far  = FCS_P2NFFT_REG_FAR_DEFAULT;

  /* init local data distribution of PNFFT:
   * local_N, local_N_start, lower_border, upper_border */
  for(int t=0; t<3; t++){
    d->local_N[t] = -1;
    d->local_N_start[t] = -1;
    d->lower_border[t] = -1;
    d->upper_border[t] = -1;
  }

  d->regkern_hat = NULL;  

  d->max_particle_move = -1;
  d->resort = d->local_num_particles = 0;
  d->gridsort_resort = FCS_GRIDSORT_RESORT_NULL;
  d->gridsort_cache = FCS_GRIDSORT_CACHE_NULL;

  *rd = d;

  return NULL;
}


static ifcs_p2nfft_data_struct* mkplan_p2nfft(
    void
    )
{
  ifcs_p2nfft_data_struct *d; 

  /* allocate the memory for the p2nfft data structure */
  d = (ifcs_p2nfft_data_struct*) malloc(sizeof(ifcs_p2nfft_data_struct));

  /* set PNFFT plan to NULL to avoid senseless references */
  d->pnfft = NULL;

  d->use_ewald = -1;

  /* initialize pointer to interpolation table */
  d->interpolation_order = 3;
  d->near_interpolation_num_nodes = 0;
  d->far_interpolation_num_nodes = 0;
  d->near_interpolation_table_potential = NULL;  
  d->near_interpolation_table_force = NULL;  
  d->far_interpolation_table_potential = NULL;  
  
  d->taylor2p_coeff = NULL;
  d->taylor2p_derive_coeff = NULL;

  d->cg_cos_coeff = NULL;
  d->cg_sin_coeff = NULL;

  /* no virial computed on default */
  d->virial = NULL;

  return d;
}

void ifcs_p2nfft_destroy(
    void *rd
    )
{
  ifcs_p2nfft_data_struct *d = (ifcs_p2nfft_data_struct *) rd;

  if (d == NULL)
    return;
 
  /* finalize PNFFT */
  if (d->pnfft != NULL)
    FCS_PNFFT(finalize)(d->pnfft, PNFFT_FREE_F_HAT|PNFFT_FREE_F|PNFFT_FREE_X|PNFFT_FREE_GRAD_F);

  /* free interpolation tables */
  if(d->near_interpolation_table_potential != NULL)
    free(d->near_interpolation_table_potential);
  if(d->near_interpolation_table_force != NULL)
    free(d->near_interpolation_table_force);
  if(d->far_interpolation_table_potential != NULL)
    free(d->far_interpolation_table_potential);

  /* free horner coefficients of 2 point Taylor polynomials */
  if(d->taylor2p_coeff != NULL)
    free(d->taylor2p_coeff);
  if(d->taylor2p_derive_coeff != NULL)
    free(d->taylor2p_derive_coeff);

  /* free cosine and sine coefficients of CG approximation */
  if(d->cg_cos_coeff != NULL)
    free(d->cg_cos_coeff);
  if(d->cg_sin_coeff != NULL)
    free(d->cg_sin_coeff);

  /* destroy precomputed Fourier coefficients */
  if(d->regkern_hat != NULL)
    FCS_PFFT(free)(d->regkern_hat);
  
  /* destroy virial */
  if(d->virial != NULL)
    free(d->virial);

  /* free Cartesian communicators */
  MPI_Comm_free(&d->cart_comm_pnfft);
  MPI_Comm_free(&d->cart_comm_3d);

  /* free mem of data struct */
  free(d);
}



#if FCS_P2NFFT_USE_3D_PROCMESH 
static void comm_create_cart_3d(
    MPI_Comm comm, MPI_Comm *comm_cart_3d, int *np
    )
{
  int size;

  MPI_Comm_size(comm, &size);
  calc_grid_sizes(size, 3, np);
  FCS_PNFFT(create_procmesh)(3, comm, np, comm_cart_3d);
}
#else
/* create a two-dimensional carteasian communicator and
 * its three-dimensional counterpart (last dimension equals 1) */
static void comm_create_cart_2d(
    MPI_Comm comm, MPI_Comm *comm_cart_2d, MPI_Comm *comm_cart_3d, int *np
    )
{
  int size;

  MPI_Comm_size(comm, &size);
  calc_grid_sizes(size, 2, np);
  np[2] = 1;

  FCS_PNFFT(create_procmesh)(2, comm, np, comm_cart_2d);
  FCS_PNFFT(create_procmesh)(3, comm, np, comm_cart_3d);
}
#endif

static void comm_get_periodicity(
    MPI_Comm comm, fcs_int *periodicity
    )
{
  int dims[3], periods[3], coords[3];
 
  /* default: no periodicity given */ 
  for(int t=0; t<3; t++)
    periodicity[t] = -1;

  if( !comm_is_cart_3d(comm) )
    return;
  
  /* for 3d cart comm use periodcity of comm */ 
  MPI_Cart_get(comm, 3,
      dims, periods, coords);
  for(int t=0; t<3; t++)
    periodicity[t] = periods[t];
}

static int comm_is_cart_3d(
    MPI_Comm comm
    )
{
  int ndims, status;

  MPI_Topo_test(comm, &status);
  if (status != MPI_CART)
    return 0;

  MPI_Cartdim_get(comm, &ndims);
  if (ndims != 3)
    return 0;

  return 1;
}


static void factorize_int(
    int n, int numfac, int *buffer
    )
{
  if(numfac == 1){
    buffer[0] = n;
    return;
  }
 
  for (int i = (int) (pow(n, 1.0/numfac)+0.5); i >= 1; --i){
    if (n % i == 0) {
      buffer[0] = i;
      factorize_int(n/i, numfac-1, buffer+1);
      return;
    }
  }
}

/* This is an stupid n^2 method for sorting n integers from largest to smallest one.
 * Since n ist small (here between 1 and 4) this should be no problem.
 * In one step we search for the minimum and copy it to the last position.
 * Then we do this recursively with the rest of the array. */
static void sort_int_size(
    int n, int *array
    )
{
  int tmp, t, t_min;
 
  if(n<=1)
    return;

  for(t=0, t_min=0; t<n; t++)
    if(array[t] < array[t_min])
      t_min = t;

  /* exchange last element and minimum */
  tmp = array[t_min];
  array[t_min] = array[n-1];
  array[n-1] = tmp;

  sort_int_size(n-1, array);
}





/** @brief Calculates the sizes of a processor grid.
 *  @param[in] nprocs The number of processors to be in the grid.
 *  @param[in] ndim The dimension of the grid.
 *  @param[out] np An array where the sizes of the grid are
 *         stored. Has to be an array of size at least @p ndim.
 */
static void calc_grid_sizes(
    fcs_int nprocs, fcs_int ndim,
    int *np
    )
{
  /* substitute for MPI_Dims_create(nprocs, ndim, np); */
  factorize_int(nprocs, ndim, np);
  sort_int_size(ndim, np);
 
#if FCS_ENABLE_DEBUG || FCS_P2NFFT_DEBUG
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(!myrank){
    fprintf(stderr, "P2NFFT_DEBUG: Procmesh: %d", np[0]);
    for(int t=1; t<ndim; t++)
      fprintf(stderr, " x %d", np[t]);
    fprintf(stderr, "\n");
  }
#endif
}

