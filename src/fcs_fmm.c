/*
  Copyright (C) 2011-2012 Rene Halver
  Copyright (C) 2016 Michael Hofmann

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



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "fcs_fmm.h"
#include "FCSCommon.h"
#include "../lib/fmm/src/fmm_cbindings.h"
#include "../lib/fmm/sl_fmm/mpi_fmm_resort.h"


#define FMM_CHECK_RETURN_RESULT(_h_, _f_)  do { \
  CHECK_HANDLE_RETURN_RESULT(_h_, _f_); \
  CHECK_METHOD_RETURN_RESULT(_h_, _f_, FCS_METHOD_FMM, "fmm"); \
  } while (0)

#define FMM_CHECK_RETURN_VAL(_h_, _f_, _v_)  do { \
  CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_); \
  CHECK_METHOD_RETURN_VAL(_h_, _f_, FCS_METHOD_FMM, "fmm", _v_); \
  } while (0)


/*
typedef struct fmm_internal_parameters_t
{
  long serroranalysis;
  long nerroranalysis;
  long pgd;
  long depth;
  long nmultipoles;
  long parabola;
  long ilinearpotentialsv;
  long iplummerpontialsv;


  _Bool firsterroranalysis;
  _Bool hugep[101];
  double fracdepth;
  double shmonopole;
  double linearodistancesv[3];
  double aoplummersv;
  double hugef[100];
} fmm_internal_parameters_t;
*/

/*
void fcs_fmm_setup_f(void *handle, fcs_int absrel, fcs_float deltaE, fcs_int dipole_correction, fcs_int *return_value)
{
  FCSResult result = fcs_fmm_setup((FCS)handle, absrel, deltaE, dipole_correction,0LL);
  if (result == FCS_RESULT_SUCCESS)
    *return_value = 0;
  else
    *return_value = fcs_result_get_return_code(result);
}
*/

/* combined setter function for all fmm parameters */
FCSResult fcs_fmm_setup(FCS handle, fcs_int absrel, fcs_float tolerance_energy, fcs_int dipole_correction, fcs_int system, fcs_int maxdepth, fcs_int unroll_limit, fcs_int load/*, fcs_int potential, fcs_float radius*/)
{
  FCSResult result;

  result = fcs_fmm_set_absrel(handle,absrel);
  CHECK_RESULT_RETURN(result);

  result = fcs_fmm_set_tolerance_energy(handle,tolerance_energy);
  CHECK_RESULT_RETURN(result);

  result = fcs_fmm_set_dipole_correction(handle, dipole_correction);
  CHECK_RESULT_RETURN(result);

  result = fcs_fmm_set_internal_tuning(handle, system);
  CHECK_RESULT_RETURN(result);

  result = fcs_fmm_set_maxdepth(handle, maxdepth);
  CHECK_RESULT_RETURN(result);

  result = fcs_fmm_set_unroll_limit(handle, unroll_limit);
  CHECK_RESULT_RETURN(result);

  result = fcs_fmm_set_balanceload(handle, load);
  CHECK_RESULT_RETURN(result);

/*
  result = fcs_fmm_set_potential(handle, potential);
  CHECK_RESULT_RETURN(result);

  result = fcs_fmm_set_dipole_correction(handle, dipole_correction);
  CHECK_RESULT_RETURN(result);
*/

  return FCS_RESULT_SUCCESS;
}


/* method to check if fmm parameters are consistent with requirements */
FCSResult fcs_fmm_check(FCS handle, fcs_int local_particles)
{
  const fcs_float *a,*b,*c;
  fcs_float norm[3],period_length;
  const fcs_int *periodicity;
  MPI_Comm comm;
  fcs_int total_particles;
  int comm_size;
  int i,p,ok;

  fcs_int absrel;
  fcs_fmm_get_absrel(handle, &absrel);
  if (absrel == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "fmm: absrel not set");

  fcs_float tolerance_energy;
  fcs_fmm_get_tolerance_energy(handle, &tolerance_energy);
  if (tolerance_energy == -1.0)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "fmm: energy tolerance not set");

  comm = fcs_get_communicator(handle);
  MPI_Comm_size(comm, &comm_size);
  total_particles = fcs_get_total_particles(handle);
  if (total_particles < comm_size)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, "fmm: there have to be at least as much particles as processes");
  
  if (local_particles <= 0)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, "fmm: each process has to receive at least one particle");
  
  a = fcs_get_box_a(handle); 
  norm[0] = fcs_norm(a);
  b = fcs_get_box_b(handle);
  norm[1] = fcs_norm(b);
  c = fcs_get_box_c(handle);
  norm[2] = fcs_norm(c);
  periodicity = fcs_get_periodicity(handle);

  p = 0; ok = 1;
  for (i = 0; i < 3; ++i)
    if (periodicity[i]) 
    {
      if (0==p) 
        period_length = norm[i];
      else 
        ok = ok && fcs_float_is_equal( period_length, norm[i] );
      p++;
    }

  if (p && !(fcs_uses_principal_axes(a,b,c))) 
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
      "fmm: with periodic boundaries, box must be arranged along principle axes");

  if (p && !ok)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
      "fmm: all periodic directions must have equal length");

  if (p && (p<3)) 
  {
    ok = 1;
    for (i = 0; i < 3; ++i)
      if (!periodicity[i]) 
        ok = ok && (norm[i] <= period_length);
    if (!ok)
      return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
        "fmm: axes in non-periodic directions must not be longer than periodic ones");
  }

  return FCS_RESULT_SUCCESS;
}


#if 0

FCSResult fcs_fmm_check(FCS handle)
{
  const fcs_float *a,*b,*c;
  const fcs_int *periodicity;
  int i,p;

  fcs_int absrel;
  fcs_fmm_get_absrel(handle, &absrel);
  if (absrel == -1)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "fmm: absrel not set");

  fcs_float deltaE;
  fcs_fmm_get_deltaE(handle, &deltaE);
  if (deltaE == -1.0)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "fmm: deltaE not set");

  a = fcs_get_box_a(handle);
  b = fcs_get_box_b(handle);
  c = fcs_get_box_c(handle);
  periodicity = fcs_get_periodicity(handle);

  p = 0;
  for (i = 0; i < 3; ++i)
    if (periodicity[i]) p++;

  switch(p)
  {
    case 0:
      return FCS_RESULT_SUCCESS;
    case 1:
      for (i = 0; i < 3; ++i)
        if (periodicity[i]) break;
      switch(i)
      {
        case 0:
          if (!( (fcs_norm(a) >= fcs_norm(b)) && (fcs_norm(a) >= fcs_norm(c)) ))
            return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
              "fmm: longest box vector not in direction of 1D-periodicity direction");


          return FCS_RESULT_SUCCESS;
        case 1:
          if (!( (fcs_norm(b) >= fcs_norm(a)) && (fcs_norm(b) >= fcs_norm(a)) ))
            return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
              "fmm: longest box vector not in direction of 1D-periodicity direction");
          return FCS_RESULT_SUCCESS;
        case 2:
          if (!( (fcs_norm(c) >= fcs_norm(a)) && (fcs_norm(c) >= fcs_norm(b)) ))
            return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
              "fmm: longest box vector not in direction of 1D-periodicity direction");
          return FCS_RESULT_SUCCESS;
      }
    case 2:
      for (i = 0; i < 3; ++i)
        if (!periodicity[i]) break;
      switch(i)
      {
        case 0:
          if (!( (fcs_float_is_equal(fcs_norm(b),fcs_norm(c))) && (fcs_norm(b) >= fcs_norm(a)) ))


            return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
              "fmm: longest box vector not in direction of 2D-periodicity direction");
          return FCS_RESULT_SUCCESS;
        case 1:
          if (!( (fcs_float_is_equal(fcs_norm(a),fcs_norm(c))) && (fcs_norm(a) >= fcs_norm(b)) ))
            return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
              "fmm: longest box vector not in direction of 2D-periodicity direction");
          return FCS_RESULT_SUCCESS;
        case 2:
          if (!( (fcs_float_is_equal(fcs_norm(a),fcs_norm(b))) && (fcs_norm(a) >= fcs_norm(c)) ))
            return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
              "fmm: longest box vector not in direction of 2D-periodicity direction");
          return FCS_RESULT_SUCCESS;
      }
    case 3:
      if (!(fcs_uses_principal_axes(a,b,c))) 
        return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, "fmm: cannot use a non-cubic system in 3D-periodicity");
      return FCS_RESULT_SUCCESS;
  }

  return FCS_RESULT_SUCCESS;
}

#endif

/* initialization function for basic fmm parameters */
FCSResult fcs_fmm_init(FCS handle)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  handle->fmm_param = malloc(sizeof(*handle->fmm_param));
  /* setting fmm parameters to invalid values (or default values, if possible) */
  handle->fmm_param->absrel = FCS_FMM_STANDARD_ERROR;
  handle->fmm_param->tolerance_energy = -1.0;
  handle->fmm_param->dipole_correction = -1;
  handle->fmm_param->potential = -1;
  handle->fmm_param->cusp_radius = -1.0;

  fcs_fmm_set_absrel( handle, FCS_FMM_CUSTOM_RELATIVE );
  fcs_fmm_set_tolerance_energy( handle, 1e-3 );
  fcs_fmm_set_dipole_correction( handle, FCS_FMM_ACTIVE_DIPOLE_CORRECTION );
  fcs_fmm_set_internal_tuning( handle, FCS_FMM_HOMOGENOUS_SYSTEM );
  fcs_fmm_set_balanceload( handle, 1 );
  fcs_fmm_set_define_loadvector( handle, 1 );
  fcs_fmm_set_maxdepth( handle, 20 );
  fcs_fmm_set_unroll_limit( handle, 9 );
  /* FCSResult result; */
  void* ptr;
  ptr = malloc(4096);
  fmm_cinit(ptr);
  fcs_set_method_context( handle, ptr );

  handle->fmm_param->wignersize = 0;
  handle->fmm_param->wignerptr = NULL;

  fcs_fmm_set_max_particle_move(handle, -1);
  fcs_fmm_set_resort(handle, 0);
  handle->fmm_param->fmm_resort = FCS_FMM_RESORT_NULL;

  handle->shift_positions = 0;

  handle->destroy = fcs_fmm_destroy;
  handle->set_tolerance = fcs_fmm_set_tolerance;
/*  handle->get_tolerance = fcs_fmm_get_tolerance;*/
  handle->set_parameter = fcs_fmm_set_parameter;
  handle->print_parameters = fcs_fmm_print_parameters;
  handle->tune = fcs_fmm_tune;
  handle->run = fcs_fmm_run;
  handle->set_compute_virial = fcs_fmm_require_virial;
  handle->get_virial = fcs_fmm_get_virial;

  handle->set_max_particle_move = fcs_fmm_set_max_particle_move;
  handle->set_resort = fcs_fmm_set_resort;
  handle->get_resort = fcs_fmm_get_resort;
  handle->get_resort_availability = fcs_fmm_get_resort_availability;
  handle->get_resort_particles = fcs_fmm_get_resort_particles;
  handle->resort_ints = fcs_fmm_resort_ints;
  handle->resort_floats = fcs_fmm_resort_floats;
  handle->resort_bytes = fcs_fmm_resort_bytes;

  return FCS_RESULT_SUCCESS;
}

/* internal fmm-specific tuning function */
FCSResult fcs_fmm_tune(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges)
{
  FCSResult result;

  fcs_int dotune;
  long long ll_tp;
  long long ll_lp;
  long long ll_absrel;
  long long ll_dip_corr;
  const fcs_int* periodicity;
  long long* ll_periodicity;
  int i;
  fcs_float tolerance_energy;
  fcs_float period_length;
  void* params;
  long long wignersize;
  long long r;
  const fcs_float* box_vector;
  void* loadptr = NULL;

  FMM_CHECK_RETURN_RESULT(handle, __func__);

  result = fcs_fmm_check(handle, local_particles);
  CHECK_RESULT_RETURN(result);

  ll_periodicity = (long long*)malloc(3*sizeof(long long));

  ll_tp = (long long)fcs_get_total_particles(handle);
  ll_lp = (long long)local_particles;
  fcs_int absrel;
  fcs_fmm_get_absrel(handle, &absrel);
  ll_absrel = absrel;
  fcs_fmm_get_tolerance_energy(handle, &tolerance_energy);
  fcs_int dip_corr;
  fcs_fmm_get_dipole_correction(handle, &dip_corr);
  ll_dip_corr = dip_corr;
  periodicity = fcs_get_periodicity(handle);
  for (i = 0; i < 3; i++)
    ll_periodicity[i] = (long long)periodicity[i];
  params = fcs_get_method_context(handle);
  box_vector = fcs_get_box_a(handle);
  period_length = fcs_norm(box_vector);

  long long ll_unroll_limit = handle->fmm_param->limit;
  long long ll_maxdepth = handle->fmm_param->maxdepth;
  long long ll_balance_load = handle->fmm_param->balance;
  long long define_loadvector = handle->fmm_param->define_loadvector;

  fcs_fmm_get_internal_tuning( handle, &dotune );
  if (dotune == FCS_FMM_INHOMOGENOUS_SYSTEM)
  {
    if (define_loadvector == 1)
    {
      handle->fmm_param->define_loadvector = 0;
      long long ll_loadvectorsize = 4*local_particles;
      fcs_float val = 1e0;
      loadptr = malloc(sizeof(ll_loadvectorsize)*ll_loadvectorsize);
      fmm_cinitload(params,loadptr,ll_loadvectorsize);
      fmm_csetload(params,val);
    }

    fmm_ctune(ll_lp,positions,charges,ll_tp,ll_absrel,tolerance_energy,ll_dip_corr,
      ll_periodicity, period_length, ll_maxdepth, ll_unroll_limit, ll_balance_load,params, &wignersize, &r);

  } else if ( dotune == FCS_FMM_HOMOGENOUS_SYSTEM )
  {
    fmm_ctunehomogen(params,&wignersize, &r);

  } else
  {
    result = fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "wrong kind of internal tuning chosen");
    return result;
  }

  if (handle->fmm_param->wignerptr == NULL || handle->fmm_param->wignersize < wignersize)
  {
    if (handle->fmm_param->wignerptr) free(handle->fmm_param->wignerptr);

    handle->fmm_param->wignersize = wignersize;
    handle->fmm_param->wignerptr = malloc(wignersize);
  }

  fmm_ccomputewigner(handle->fmm_param->wignerptr,params,dotune);

  free(ll_periodicity);
  if (loadptr) free(loadptr);

  if (r == 0)
  {
    result = fcs_result_create(FCS_ERROR_FORTRAN_CALL, __func__, "error in fmm_ctune (FORTRAN)");
    return result;
  }
  
  return FCS_RESULT_SUCCESS;
}

int fcs_mpi_fmm_sort_front_part, fcs_mpi_fmm_sort_back_part, fcs_mpi_fmm_sort_front_merge_presorted;

/* internal fmm-specific run function */
FCSResult fcs_fmm_run(FCS handle, fcs_int local_particles,
                      fcs_float *positions, fcs_float *charges, 
                      fcs_float *field, fcs_float *potentials)
{
  FCSResult result;

  long long ll_tp;
  long long ll_lp;
  long long ll_absrel;
  long long ll_dip_corr;
  const fcs_int* periodicity;
  long long* ll_periodicity;
  int i;
  fcs_float tolerance_energy;
  fcs_float period_length;
  void* params;
  fcs_int dotune;
  long long r;
  const fcs_float* box_vector;
  void* loadptr = NULL;

  FMM_CHECK_RETURN_RESULT(handle, __func__);

  result = fcs_fmm_check(handle, local_particles);
  CHECK_RESULT_RETURN(result);

  ll_periodicity = (long long*)malloc(3*sizeof(long long));

  ll_tp = (long long)fcs_get_total_particles(handle);
  ll_lp = (long long)local_particles;

  fcs_int absrel;
  fcs_fmm_get_absrel(handle, &absrel);
  ll_absrel = absrel;

  fcs_fmm_get_tolerance_energy(handle, &tolerance_energy);

  fcs_int dip_corr;
  fcs_fmm_get_dipole_correction(handle, &dip_corr);

  ll_dip_corr = dip_corr;

  periodicity = fcs_get_periodicity(handle);
  for (i = 0; i < 3; i++)
    ll_periodicity[i] = (long long)periodicity[i];
  params = fcs_get_method_context(handle);
  box_vector = fcs_get_box_a(handle);
  period_length = fcs_norm(box_vector);

  fcs_fmm_get_internal_tuning( handle, &dotune );
  if (dotune != FCS_FMM_INHOMOGENOUS_SYSTEM && dotune != FCS_FMM_HOMOGENOUS_SYSTEM)
  {
    result = fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "wrong kind of internal tuning chosen");
    return result;
  }

  long long ll_unroll_limit = handle->fmm_param->limit;
  long long ll_maxdepth = handle->fmm_param->maxdepth;
  long long ll_balance_load = handle->fmm_param->balance;
  long long define_loadvector = handle->fmm_param->define_loadvector;
 
  if (define_loadvector == 1)
  {
    handle->fmm_param->define_loadvector = 0;
    long long ll_loadvectorsize = ll_loadvectorsize = 229074;/*local_particles*4;*/
    fcs_float val = 1e0;
    loadptr = malloc(sizeof(ll_loadvectorsize)*ll_loadvectorsize);
    fmm_cinitload(params,loadptr,ll_loadvectorsize);
    fmm_csetload(params,val);
  }

  int old_fcs_mpi_fmm_sort_front_part = fcs_mpi_fmm_sort_front_part;

  int comm_size;
  MPI_Comm_size(fcs_get_communicator(handle), &comm_size);
  fcs_float max_merge_move = fcs_pow(fcs_norm(fcs_get_box_a(handle)) * fcs_norm(fcs_get_box_b(handle)) * fcs_norm(fcs_get_box_c(handle)) / comm_size, 1.0 / 3.0);

  fcs_float max_particle_move;
  MPI_Allreduce(&handle->fmm_param->max_particle_move, &max_particle_move, 1, FCS_MPI_FLOAT, MPI_MAX, fcs_get_communicator(handle));

  if (max_particle_move >= 0 && max_particle_move < max_merge_move)
  {
/*    fmm_csetpresorted(params, 1);*/
    fcs_mpi_fmm_sort_front_part = 0;
    fcs_mpi_fmm_sort_front_merge_presorted = 1;

  } else
  {
/*    fmm_csetpresorted(params, 0);*/
    fcs_mpi_fmm_sort_front_merge_presorted = 0;
  }

  fcs_fmm_resort_destroy(&handle->fmm_param->fmm_resort);
  if (handle->fmm_param->resort) fcs_fmm_resort_create(&handle->fmm_param->fmm_resort, local_particles, fcs_get_communicator(handle));
  fmm_cinitresort(params, handle->fmm_param->fmm_resort);
  fmm_csetresort(params, (long long) handle->fmm_param->resort);

  fmm_crun(ll_lp,positions,charges,potentials,field,handle->fmm_param->virial,ll_tp,ll_absrel,tolerance_energy,
    ll_dip_corr, ll_periodicity, period_length, dotune, ll_maxdepth,ll_unroll_limit,ll_balance_load,params, &r);

  fcs_mpi_fmm_sort_front_part = old_fcs_mpi_fmm_sort_front_part;

  free(ll_periodicity);
  if (loadptr) free(loadptr);

  if (r == 0)
  {
    result = fcs_result_create(FCS_ERROR_FORTRAN_CALL, __func__, "error in fmm_run (FORTRAN)");
    return result;
  }

  return FCS_RESULT_SUCCESS;
}

/* clean-up function for fmm */
FCSResult fcs_fmm_destroy(FCS handle)
{
  fcs_int dotune;

  FMM_CHECK_RETURN_RESULT(handle, __func__);

  fcs_fmm_get_internal_tuning( handle, &dotune);

  fmm_cfinalize(fcs_get_method_context(handle),dotune);

  if (handle->fmm_param->wignerptr) free(handle->fmm_param->wignerptr);

  free(fcs_get_method_context(handle));

  fcs_fmm_resort_destroy(&handle->fmm_param->fmm_resort);

  free(handle->fmm_param);
  handle->fmm_param = NULL;

  return FCS_RESULT_SUCCESS;
}

/******************************************************************************************************
 *
 *            Setter and Getter functions for fmm parameters
 *
 ******************************************************************************************************/


/* setter function for fmm parameter absrel */
FCSResult fcs_fmm_set_absrel(FCS handle, fcs_int choice)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (choice < 0 || choice > 2)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"error switch has to be between 0 and 2");

  handle->fmm_param->absrel = choice;

  return FCS_RESULT_SUCCESS;
}

/* getter function for fmm parameter absrel */
FCSResult fcs_fmm_get_absrel(FCS handle, fcs_int *absrel)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!absrel)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for absrel");

  *absrel = handle->fmm_param->absrel;

  return FCS_RESULT_SUCCESS;
}

/* setter function for fmm parameter internal tuning */
FCSResult fcs_fmm_set_internal_tuning(FCS handle, fcs_int system)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (system != FCS_FMM_HOMOGENOUS_SYSTEM && system != FCS_FMM_INHOMOGENOUS_SYSTEM)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"unknown system type chosen, use either: FCS_FMM_HOMOGENOUS_SYSTEM or FCS_FMM_INHOMOGENOUS_SYSTEM");

  handle->fmm_param->system = system;

  return FCS_RESULT_SUCCESS;
}

/* getter function for fmm parameter internal tuning */
FCSResult fcs_fmm_get_internal_tuning(FCS handle, fcs_int *system)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!system)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for internal tuning parameter");

  *system = handle->fmm_param->system;

  return FCS_RESULT_SUCCESS;
}

/* setter function for fmm parameter energy tolerance (deltaE) */
FCSResult fcs_fmm_set_tolerance_energy(FCS handle, fcs_float tolerance_energy)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (tolerance_energy <= 0.0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"chosen error not valid, has to be larger than zero");

  handle->fmm_param->tolerance_energy = tolerance_energy;

  return FCS_RESULT_SUCCESS;
}



/* getter function for fmm parameter energy tolerance (deltaE) */
FCSResult fcs_fmm_get_tolerance_energy(FCS handle, fcs_float *tolerance_energy)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!tolerance_energy)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for tolerance_energy");

  *tolerance_energy = handle->fmm_param->tolerance_energy;

  return FCS_RESULT_SUCCESS;
}


/* setter function for fmm parameter dipole correction */
FCSResult fcs_fmm_set_dipole_correction(FCS handle, fcs_int dipole_correction)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!(dipole_correction == FCS_FMM_NO_DIPOLE_CORRECTION || dipole_correction == FCS_FMM_STANDARD_DIPOLE_CORRECTION || dipole_correction == FCS_FMM_ACTIVE_DIPOLE_CORRECTION))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"invalid dipole correction chosen");

  handle->fmm_param->dipole_correction = dipole_correction;

  return FCS_RESULT_SUCCESS;
}

/* getter function for fmm parameter dipole correction */
FCSResult fcs_fmm_get_dipole_correction(FCS handle, fcs_int *dipole_correction)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!dipole_correction)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for dipole_correction");

  *dipole_correction = handle->fmm_param->dipole_correction;

  return FCS_RESULT_SUCCESS;
}

/* setter function for fmm parameter potential */
FCSResult fcs_fmm_set_potential(FCS handle, fcs_int potential)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!(potential == FCS_FMM_COULOMB || potential == FCS_FMM_CUSP))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"invalid (fmm) potential chosen");

  handle->fmm_param->potential = potential;

  return FCS_RESULT_SUCCESS;
}

/* getter function for fmm parameter potential */
FCSResult fcs_fmm_get_potential(FCS handle, fcs_int *potential)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!potential)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for potential");

  *potential = handle->fmm_param->potential;

  return FCS_RESULT_SUCCESS;
}

/* setter function for maximum fmm tree depth */
FCSResult fcs_fmm_set_maxdepth(FCS handle, fcs_int depth)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  handle->fmm_param->maxdepth = depth;

  return FCS_RESULT_SUCCESS;
}

/* getter function for maximum fmm tree depth */
FCSResult fcs_fmm_get_maxdepth(FCS handle, fcs_int *depth)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!depth)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for depth");

  *depth = handle->fmm_param->maxdepth;

  return FCS_RESULT_SUCCESS;
}

/* setter function for maximum fmm unroll limit */
FCSResult fcs_fmm_set_unroll_limit(FCS handle, fcs_int limit)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (limit < 0) limit = 0;
  if (limit > 50) limit = 50;

  handle->fmm_param->limit = limit;

  return FCS_RESULT_SUCCESS;
}

/* getter function for fmm unroll limit */
FCSResult fcs_fmm_get_unroll_limit(FCS handle, fcs_int *limit)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!limit)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for limit");

  *limit = handle->fmm_param->limit;

  return FCS_RESULT_SUCCESS;
}

/* setter function for status of fmm load balancing */
FCSResult fcs_fmm_set_balanceload(FCS handle, fcs_int load)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);
  
  handle->fmm_param->balance = (load != 0)?1:0;

  return FCS_RESULT_SUCCESS;
}

/* setter function for status of fmm loadvector definition */
FCSResult fcs_fmm_set_define_loadvector(FCS handle, fcs_int define_loadvector)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (define_loadvector != 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"invalid (fmm) define_loadvector chosen (only value 1 is current supported)");

  handle->fmm_param->define_loadvector = define_loadvector;

  return FCS_RESULT_SUCCESS;
}

/* getter function for status of fmm define loadvector */
FCSResult fcs_fmm_get_define_loadvector(FCS handle, fcs_int *define_loadvector)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!define_loadvector)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for define_loadvector");

  *define_loadvector = handle->fmm_param->define_loadvector;

  return FCS_RESULT_SUCCESS;
}


/* getter function for status of fmm load balancing */
FCSResult fcs_fmm_get_balanceload(FCS handle, fcs_int *load)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!load)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for load");

  *load = handle->fmm_param->balance;

  return FCS_RESULT_SUCCESS;
}

/* setter function for fmm parameter cusp_radius */
FCSResult fcs_fmm_set_cusp_radius(FCS handle, fcs_float radius)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (radius < 0.0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"cusp radius must be non-negative");

  handle->fmm_param->cusp_radius = radius;

  return FCS_RESULT_SUCCESS;
}

/* getter function for fmm parameter cusp radius */
FCSResult fcs_fmm_get_cusp_radius(FCS handle, fcs_float *cusp_radius)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!cusp_radius)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for cusp_radius");

  *cusp_radius = handle->fmm_param->cusp_radius;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_fmm_require_virial(FCS handle, fcs_int flag) {

  FMM_CHECK_RETURN_RESULT(handle, __func__);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_fmm_get_virial(FCS handle, fcs_float *virial) {

  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (!virial)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,__func__,"null pointer supplied for virial"); 

  fcs_int i;
  for (i=0; i < 9; i++)
    virial[i] = handle->fmm_param->virial[i];

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY)
  {
    fcs_fmm_set_absrel(handle, FCS_FMM_CUSTOM_ABSOLUTE);
    fcs_fmm_set_tolerance_energy(handle, tolerance);

  } else if (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY_REL)
  {
    fcs_fmm_set_absrel(handle, FCS_FMM_CUSTOM_RELATIVE);
    fcs_fmm_set_tolerance_energy(handle, tolerance);

  } else
  {
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, __func__, "Unsupported tolerance type. FMM only supports FCS_TOLERANCE_TYPE_ENERGY and FCS_TOLERANCE_TYPE_ENERGY_REL.");
  }
  
  return FCS_RESULT_SUCCESS;
}


/*FCSResult fcs_fmm_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance)
{
}*/


FCSResult fcs_fmm_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_absrel",            fmm_set_absrel,            FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_tolerance_energy",  fmm_set_tolerance_energy,  FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_dipole_correction", fmm_set_dipole_correction, FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_potential",         fmm_set_potential,         FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_cusp_radius",       fmm_set_cusp_radius,       FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_internal_tuning",   fmm_set_internal_tuning,   FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_maxdepth",          fmm_set_maxdepth,          FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_unroll_limit",      fmm_set_unroll_limit,      FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_balanceload",       fmm_set_balanceload,       FCS_PARSE_VAL(fcs_int));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_print_parameters(FCS handle)
{
  fcs_int absrel;
  fcs_float tolerance_energy;
  fcs_int dipole_correction;
  fcs_int tuning;
  fcs_int maxdepth;
  fcs_int limit;
  fcs_int load;

  FMM_CHECK_RETURN_RESULT(handle, __func__);

  fcs_fmm_get_absrel(handle, &absrel);
  fcs_fmm_get_tolerance_energy(handle, &tolerance_energy);
  fcs_fmm_get_dipole_correction(handle, &dipole_correction);
  fcs_fmm_get_internal_tuning(handle, &tuning);
  fcs_fmm_get_balanceload(handle, &load);
  fcs_fmm_get_maxdepth(handle, &maxdepth);
  fcs_fmm_get_unroll_limit(handle, &limit);

  printf("fmm absrel: %" FCS_LMOD_INT "d\n", absrel);
  printf("fmm tolerance value: %e\n", tolerance_energy);
  printf("fmm dipole correction: %" FCS_LMOD_INT "d\n", dipole_correction);
  printf("fmm internal tuning: %c\n", (tuning)?'T':'F');
  printf("fmm maxdepth: %" FCS_LMOD_INT "d\n", maxdepth);
  printf("fmm unroll limit: %" FCS_LMOD_INT "d\n", limit);
  printf("fmm internal balance load: %c\n", (load)?'T':'F');
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  handle->fmm_param->max_particle_move = max_particle_move;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_set_resort(FCS handle, fcs_int resort)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  handle->fmm_param->resort = resort;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_get_resort(FCS handle, fcs_int *resort)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  *resort = handle->fmm_param->resort;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_get_resort_availability(FCS handle, fcs_int *availability)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  if (handle->fmm_param->fmm_resort != FCS_FMM_RESORT_NULL) *availability = fcs_resort_is_available(handle->fmm_param->fmm_resort->resort);
  else *availability = 0;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  *resort_particles = fcs_resort_get_original_particles(handle->fmm_param->fmm_resort->resort);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  fcs_resort_resort_ints(handle->fmm_param->fmm_resort->resort, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  fcs_resort_resort_floats(handle->fmm_param->fmm_resort->resort, src, dst, n, comm);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_fmm_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  FMM_CHECK_RETURN_RESULT(handle, __func__);

  fcs_resort_resort_bytes(handle->fmm_param->fmm_resort->resort, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}
