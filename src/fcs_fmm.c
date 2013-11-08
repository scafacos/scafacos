/*
  Copyright (C) 2011-2012 Rene Halver

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
  if (NULL == result)
    *return_value = 0;
  else
    *return_value = fcsResult_getReturnCode(result);
}
*/

/* combined setter function for all fmm parameters */
extern FCSResult fcs_fmm_setup(FCS handle, fcs_int absrel, fcs_float tolerance_value, fcs_int dipole_correction, long long system, long long maxdepth, long long unroll_limit, long long load/*, fcs_int potential, fcs_float radius*/)
{
/*  char* fnc_name = "fcs_fmm_setup"; */
  FCSResult result;

  result = fcs_fmm_set_absrel(handle,absrel);
  if (result != NULL)
    return result;
  result = fcs_fmm_set_tolerance_energy(handle,tolerance_value);
  if (result != NULL)
    return result;
  result = fcs_fmm_set_dipole_correction(handle, dipole_correction);
  if (result != NULL)
    return result;
  result = fcs_fmm_set_internal_tuning(handle, system);
  if (result != NULL)
    return result;
  result = fcs_fmm_set_maxdepth(handle, maxdepth);
  if (result != NULL)
    return result;
  result = fcs_fmm_set_unroll_limit(handle, unroll_limit);
  if (result != NULL)
    return result;
  result = fcs_fmm_set_balanceload(handle, load);
  if (result != NULL)
    return result;
  return (FCSResult)NULL;
  
/*
  result = fcs_fmm_set_potential(handle, potential);
  if (result != NULL)
    return result;
  result = fcs_fmm_set_dipole_correction(handle, dipole_correction);
  if (result != NULL)
    return result;
  return NULL;
*/
}


/* method to check if fmm parameters are consistent with requirements */
extern FCSResult fcs_fmm_check(FCS handle, fcs_int local_particles)
{
  char* fnc_name = "fcs_fmm_check";
  fcs_float *a,*b,*c,norm[3],period_length;
  fcs_int *periodicity;
  MPI_Comm comm;
  fcs_int total_particles;
  int comm_size;
  int i,p,ok;

  fcs_int absrel;
  fcs_fmm_get_absrel(handle, &absrel);
  if (absrel == -1)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "fmm: absrel not set");

  fcs_float tolerance_value;
  fcs_fmm_get_tolerance_energy(handle, &tolerance_value);
  if (tolerance_value == -1.0)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "fmm: energy tolerance not set");

  comm = fcs_get_communicator(handle);
  MPI_Comm_size(comm, &comm_size);
  total_particles = fcs_get_total_particles(handle);
  if (total_particles < comm_size)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "fmm: there have to be at least as much particles as processes");
  
  if (local_particles <= 0)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "fmm: each process has to receive at least one particle");
  
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
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
      "fmm: with periodic boundaries, box must be arranged along principle axes");

  if (p && !ok)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
      "fmm: all periodic directions must have equal length");

  if (p && (p<3)) 
  {
    ok = 1;
    for (i = 0; i < 3; ++i)
      if (!periodicity[i]) 
        ok = ok && (norm[i] <= period_length);
    if (!ok)
      return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
        "fmm: axes in non-periodic directions must not be longer than periodic ones");
  }

  return NULL;
}

// /* method to check if fmm parameters are entered into checked FCS */
// extern FCSResult fcs_fmm_check(FCS handle)
// {
//   char* fnc_name = "fcs_fmm_check";
//   fcs_float *a,*b,*c;
//   fcs_int *periodicity;
//   int i,p;
// 
//   fcs_int absrel;
//   fcs_fmm_get_absrel(handle, &absrel);
//   if (absrel == -1)
//     return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "fmm: absrel not set");
// 
//   fcs_float deltaE;
//   fcs_fmm_get_deltaE(handle, &deltaE);
//   if (deltaE == -1.0)
//     return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "fmm: deltaE not set");
//   a = fcs_get_box_a(handle);
//   b = fcs_get_box_b(handle);
// 
// 
//   c = fcs_get_box_c(handle);
//   periodicity = fcs_get_periodicity(handle);
// 
//   p = 0;
//   for (i = 0; i < 3; ++i)
//     if (periodicity[i]) p++;
//   switch(p)
//   {
//     case 0:
//       return NULL;
//     case 1:
//       for (i = 0; i < 3; ++i)
//         if (periodicity[i]) break;
//       switch(i)
//       {
//         case 0:
//           if (!( (fcs_norm(a) >= fcs_norm(b)) && (fcs_norm(a) >= fcs_norm(c)) ))
//             return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
//               "fmm: longest box vector not in direction of 1D-periodicity direction");
// 
// 
//           return NULL;
//         case 1:
//           if (!( (fcs_norm(b) >= fcs_norm(a)) && (fcs_norm(b) >= fcs_norm(a)) ))
//             return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
//               "fmm: longest box vector not in direction of 1D-periodicity direction");
//           return NULL;
//         case 2:
//           if (!( (fcs_norm(c) >= fcs_norm(a)) && (fcs_norm(c) >= fcs_norm(b)) ))
//             return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
//               "fmm: longest box vector not in direction of 1D-periodicity direction");
//           return NULL;
//       }
//     case 2:
//       for (i = 0; i < 3; ++i)
//         if (!periodicity[i]) break;
//       switch(i)
//       {
//         case 0:
//           if (!( (fcs_float_is_equal(fcs_norm(b),fcs_norm(c))) && (fcs_norm(b) >= fcs_norm(a)) ))
// 
// 
//             return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
//               "fmm: longest box vector not in direction of 2D-periodicity direction");
//           return NULL;
//         case 1:
//           if (!( (fcs_float_is_equal(fcs_norm(a),fcs_norm(c))) && (fcs_norm(a) >= fcs_norm(b)) ))
//             return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
//               "fmm: longest box vector not in direction of 2D-periodicity direction");
//           return NULL;
//         case 2:
//           if (!( (fcs_float_is_equal(fcs_norm(a),fcs_norm(b))) && (fcs_norm(a) >= fcs_norm(c)) ))
//             return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
//               "fmm: longest box vector not in direction of 2D-periodicity direction");
//           return NULL;
//       }
//     case 3:
//       if (!(fcs_uses_principal_axes(a,b,c))) 
//         return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "fmm: cannot use a non-cubic system in 3D-periodicity");
//       return NULL;
//   }
// 
//   return NULL;
// }

/* initialization function for basic fmm parameters */
FCSResult fcs_fmm_init(FCS handle)
{
  fcs_fmm_set_absrel( handle, FCS_FMM_CUSTOM_RELATIVE );
  fcs_fmm_set_tolerance_energy( handle, 1e-3 );
  fcs_fmm_set_dipole_correction( handle, FCS_FMM_ACTIVE_DIPOLE_CORRECTION );
  fcs_fmm_set_internal_tuning( handle, FCS_FMM_HOMOGENOUS_SYSTEM );
  fcs_fmm_set_balanceload( handle, 1ll );
  fcs_fmm_set_define_loadvector( handle, 1ll );
  fcs_fmm_set_maxdepth( handle, 20ll );
  fcs_fmm_set_unroll_limit( handle, 9ll );
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

  handle->set_max_particle_move = fcs_fmm_set_max_particle_move;
  handle->set_resort = fcs_fmm_set_resort;
  handle->get_resort = fcs_fmm_get_resort;
  handle->get_resort_availability = fcs_fmm_get_resort_availability;
  handle->get_resort_particles = fcs_fmm_get_resort_particles;
  handle->resort_ints = fcs_fmm_resort_ints;
  handle->resort_floats = fcs_fmm_resort_floats;
  handle->resort_bytes = fcs_fmm_resort_bytes;

  return NULL;
}

/* internal fmm-specific tuning function */
FCSResult fcs_fmm_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges)
{
  char* fnc_name = "fcs_fmm_tune";
  FCSResult result;

  long long dotune;
  long long ll_tp;
  long long ll_lp;
  long long ll_absrel;
  long long ll_dip_corr;
  fcs_int* periodicity;
  long long* ll_periodicity;
  int i;
  fcs_float tolerance_value;
  fcs_float period_length;
  void* params;
  long long wignersize;
  long long r;
  fcs_float* box_vector;
  void* loadptr = NULL;

  result = fcs_fmm_check(handle, local_particles);
  if (result != NULL)
    return result;

  ll_periodicity = (long long*)malloc(3*sizeof(long long));

  ll_tp = (long long)fcs_get_total_particles(handle);
  ll_lp = (long long)local_particles;
  fcs_int absrel;
  fcs_fmm_get_absrel(handle, &absrel);
  ll_absrel = absrel;
  fcs_fmm_get_tolerance_energy(handle, &tolerance_value);
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

    fmm_ctune(ll_lp,positions,charges,ll_tp,ll_absrel,tolerance_value,ll_dip_corr,
      ll_periodicity, period_length, ll_maxdepth, ll_unroll_limit, ll_balance_load,params, &wignersize, &r);

  } else if ( dotune == FCS_FMM_HOMOGENOUS_SYSTEM )
  {
    fmm_ctunehomogen(params,&wignersize, &r);

  } else
  {
    result = fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "wrong kind of internal tuning chosen");
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
    result = fcsResult_create(FCS_FORTRAN_CALL_ERROR, fnc_name, "error in fmm_ctune (FORTRAN)");
    return result;
  }
  
  return NULL;
}

extern int fcs_mpi_fmm_sort_front_part, fcs_mpi_fmm_sort_back_part, fcs_mpi_fmm_sort_front_merge_presorted;

/* internal fmm-specific run function */
FCSResult fcs_fmm_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles, 
                      fcs_float *positions, fcs_float *charges, 
                      fcs_float *field, fcs_float *potentials)
{
  char* fnc_name = "fcs_fmm_run";
  FCSResult result;

  long long ll_tp;
  long long ll_lp;
  long long ll_absrel;
  long long ll_dip_corr;
  fcs_int* periodicity;
  long long* ll_periodicity;
  int i;
  fcs_float tolerance_value;
  fcs_float period_length;
  void* params;
  long long dotune;
  long long r;
  fcs_float* box_vector;
  void* loadptr = NULL;

  result = fcs_fmm_check(handle, local_particles);
  if (result != NULL)
    return result;

  ll_periodicity = (long long*)malloc(3*sizeof(long long));

  ll_tp = (long long)fcs_get_total_particles(handle);
  ll_lp = (long long)local_particles;

  fcs_int absrel;
  fcs_fmm_get_absrel(handle, &absrel);
  ll_absrel = absrel;

  fcs_fmm_get_tolerance_energy(handle, &tolerance_value);

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
    result = fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "wrong kind of internal tuning chosen");
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

  fmm_crun(ll_lp,positions,charges,potentials,field,handle->fmm_param->virial,ll_tp,ll_absrel,tolerance_value,
    ll_dip_corr, ll_periodicity, period_length, dotune, ll_maxdepth,ll_unroll_limit,ll_balance_load,params, &r);

  fcs_mpi_fmm_sort_front_part = old_fcs_mpi_fmm_sort_front_part;

  free(ll_periodicity);
  if (loadptr) free(loadptr);

  if (r == 0)
  {
    result = fcsResult_create(FCS_FORTRAN_CALL_ERROR, fnc_name, "error in fmm_run (FORTRAN)");
    return result;
  }

  return NULL;
}

/* clean-up function for fmm */
extern FCSResult fcs_fmm_destroy(FCS handle)
{
/*  char* fnc_name = "fcs_fmm_destroy"; */
  long long dotune;

  fcs_fmm_get_internal_tuning( handle, &dotune);
  fmm_cfinalize(fcs_get_method_context(handle),dotune);

  if (handle->fmm_param->wignerptr) free(handle->fmm_param->wignerptr);

  free(fcs_get_method_context(handle));

  fcs_fmm_resort_destroy(&handle->fmm_param->fmm_resort);

  return NULL;
}

/******************************************************************************************************
 *
 *            Setter and Getter functions for fmm parameters
 *
 ******************************************************************************************************/


/* setter function for fmm parameter absrel */
extern FCSResult fcs_fmm_set_absrel(FCS handle, fcs_int choice)
{
  char* fnc_name = "fcs_fmm_set_absrel";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");
  if (choice < 0 || choice > 2)
    return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"error switch has to be between 0 and 2");
  else
  {
    handle->fmm_param->absrel = choice;
    return NULL;
  }
}

/* getter function for fmm parameter absrel */
extern FCSResult fcs_fmm_get_absrel(FCS handle, fcs_int *absrel)
{
  char* fnc_name = "fcs_fmm_get_absrel";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!absrel)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for absrel");

  *absrel = handle->fmm_param->absrel;
  return NULL;
}

/* setter function for fmm parameter internal tuning */
extern FCSResult fcs_fmm_set_internal_tuning(FCS handle, long long system)
{
  char* fnc_name = "fcs_fmm_set_internal_tuning";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");
  if (system != FCS_FMM_HOMOGENOUS_SYSTEM && system != FCS_FMM_INHOMOGENOUS_SYSTEM)
    return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"unknown system type chosen, use either: FCS_FMM_HOMOGENOUS_SYSTEM or FCS_FMM_INHOMOGENOUS_SYSTEM");

  handle->fmm_param->system = system;
  return NULL;
}

/* getter function for fmm parameter internal tuning */
extern FCSResult fcs_fmm_get_internal_tuning(FCS handle, long long *system)
{
  char* fnc_name = "fcs_fmm_get_internal_tuning";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!system)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for internal tuning parameter");

  *system = handle->fmm_param->system;
  return NULL;
}

/* setter function for fmm parameter energy tolerance (deltaE) */
extern FCSResult fcs_fmm_set_tolerance_energy(FCS handle, fcs_float tolerance_value)
{
  char* fnc_name = "fcs_fmm_set_tolerance_energy";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
        return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");
  if (tolerance_value <= 0.0)
    return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"chosen error not valid, has to be larger than zero");

  handle->fmm_param->tolerance_value = tolerance_value;
  return NULL;
}



/* getter function for fmm parameter energy tolerance (deltaE) */
extern FCSResult fcs_fmm_get_tolerance_energy(FCS handle, fcs_float *tolerance_value)
{
  char* fnc_name = "fcs_fmm_get_tolerance_energy";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!tolerance_value)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for tolerance_value");

  *tolerance_value = handle->fmm_param->tolerance_value;
  return NULL;
}


/* setter function for fmm parameter dipole correction */
extern FCSResult fcs_fmm_set_dipole_correction(FCS handle, fcs_int dipole_correction)
{
  char* fnc_name = "fcs_fmm_set_dipole_correction";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");
  if (!(dipole_correction == FCS_FMM_NO_DIPOLE_CORRECTION || dipole_correction == FCS_FMM_STANDARD_DIPOLE_CORRECTION || dipole_correction == FCS_FMM_ACTIVE_DIPOLE_CORRECTION))
    return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"invalid dipole correction chosen");

  handle->fmm_param->dipole_correction = dipole_correction;
  return NULL;
}

/* getter function for fmm parameter dipole correction */
extern FCSResult fcs_fmm_get_dipole_correction(FCS handle, fcs_int *dipole_correction)
{
  char* fnc_name = "fcs_fmm_get_dipole_correction";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!dipole_correction)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for dipole_correction");

  *dipole_correction = handle->fmm_param->dipole_correction;
  return NULL;
}

/* setter function for fmm parameter potential */
extern FCSResult fcs_fmm_set_potential(FCS handle, fcs_int potential)
{
  char* fnc_name = "fcs_fmm_set_potential";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
        return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");
  if (!(potential == FCS_FMM_COULOMB || potential == FCS_FMM_CUSP))
    return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"invalid (fmm) potential chosen");

  handle->fmm_param->potential = potential;
  return NULL;
}

/* getter function for fmm parameter potential */
extern FCSResult fcs_fmm_get_potential(FCS handle, fcs_int *potential)
{
  char* fnc_name = "fcs_fmm_get_potential";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!potential)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for potential");

  *potential = handle->fmm_param->potential;
  return NULL;
}

/* setter function for maximum fmm tree depth */
extern FCSResult fcs_fmm_set_maxdepth(FCS handle, long long depth)
{
  char* fnc_name = "fcs_fmm_set_maxdepth";
  fcs_int int_depth = depth;
  
  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");

  handle->fmm_param->maxdepth = int_depth;
  return NULL;
}

/* getter function for maximum fmm tree depth */
extern FCSResult fcs_fmm_get_maxdepth(FCS handle, long long *depth)
{
  char* fnc_name = "fcs_fmm_get_maxdepth";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!depth)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for depth");

  *depth = handle->fmm_param->maxdepth;
  return NULL;
}

/* setter function for maximum fmm unroll limit */
extern FCSResult fcs_fmm_set_unroll_limit(FCS handle, long long limit)
{
  char* fnc_name = "fcs_fmm_set_unroll_limit";
  fcs_int int_limit = limit;
  
  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");

  int_limit = (0 > limit)?0:limit;
  int_limit = (50 < limit)?50:limit;
  handle->fmm_param->limit = int_limit;
  return NULL;
}

/* getter function for fmm unroll limit */
extern FCSResult fcs_fmm_get_unroll_limit(FCS handle, long long *limit)
{
  char* fnc_name = "fcs_fmm_get_unroll_limit";
  
  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!limit)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for limit");

  *limit = handle->fmm_param->limit;
  return NULL;
}

/* setter function for status of fmm load balancing */
extern FCSResult fcs_fmm_set_balanceload(FCS handle, long long load)
{
  char* fnc_name = "fcs_fmm_set_balanceload";
  fcs_int int_load = load;
  
  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");

  int_load = (load != 0)?1:0;
  handle->fmm_param->balance = int_load;
  return NULL;
}

/* setter function for status of fmm loadvector definition */
extern FCSResult fcs_fmm_set_define_loadvector(FCS handle, long long define_loadvector)
{
  char* fnc_name = "fcs_fmm_set_define_loadvector";
  fcs_int int_define_loadvector = define_loadvector;

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");

  handle->fmm_param->define_loadvector = 1;//int_define_loadvector;
  return NULL;
}

/* getter function for status of fmm define loadvector */
extern FCSResult fcs_fmm_get_define_loadvector(FCS handle, long long *define_loadvector)
{
  char* fnc_name = "fcs_fmm_get_define_loadvector";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!define_loadvector)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for define_loadvector");

  *define_loadvector = (long long) handle->fmm_param->define_loadvector;
  return NULL;
}


/* getter function for status of fmm load balancing */
extern FCSResult fcs_fmm_get_balanceload(FCS handle, long long *load)
{
  char* fnc_name = "fcs_fmm_get_balanceload";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!load)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for load");

  *load = handle->fmm_param->balance;
  return NULL;
}

/* setter function for fmm parameter cusp_radius */
extern FCSResult fcs_fmm_set_cusp_radius(FCS handle, fcs_float radius)
{
  char* fnc_name = "fcs_fmm_set_radius";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_FMM)
        return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"fmm\")");
  if (radius < 0.0)
    return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"cusp radius must be non-negative");

  handle->fmm_param->cusp_radius = radius;
  return NULL;
}

/* getter function for fmm parameter cusp radius */
extern FCSResult fcs_fmm_get_cusp_radius(FCS handle, fcs_float *cusp_radius)
{
  char* fnc_name = "fcs_fmm_get_radius";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!cusp_radius)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for cusp_radius");

  *cusp_radius = handle->fmm_param->cusp_radius;
  return NULL;
}

FCSResult fcs_fmm_require_virial(FCS handle, fcs_int flag) { return NULL; }

FCSResult fcs_fmm_get_virial(FCS handle, fcs_float *virial) {
  char* fnc_name = "fcs_fmm_get_virial";

  if (!handle || !handle->fmm_param)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (!virial)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for virial"); 

  fcs_int i;
  for (i=0; i < 9; i++)
    virial[i] = handle->fmm_param->virial[i];
  return NULL;
}


FCSResult fcs_fmm_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  handle->fmm_param->max_particle_move = max_particle_move;

  return NULL;
}


FCSResult fcs_fmm_set_resort(FCS handle, fcs_int resort)
{
  handle->fmm_param->resort = resort;

  return NULL;
}


FCSResult fcs_fmm_get_resort(FCS handle, fcs_int *resort)
{
  *resort = handle->fmm_param->resort;

  return NULL;
}


FCSResult fcs_fmm_get_resort_availability(FCS handle, fcs_int *availability)
{
  if (handle->fmm_param->fmm_resort != FCS_FMM_RESORT_NULL) *availability = fcs_resort_is_available(handle->fmm_param->fmm_resort->resort);
  else *availability = 0;

  return NULL;
}


FCSResult fcs_fmm_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  *resort_particles = fcs_resort_get_original_particles(handle->fmm_param->fmm_resort->resort);

  return NULL;
}


FCSResult fcs_fmm_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  fcs_resort_resort_ints(handle->fmm_param->fmm_resort->resort, src, dst, n, comm);
  
  return NULL;
}


FCSResult fcs_fmm_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  fcs_resort_resort_floats(handle->fmm_param->fmm_resort->resort, src, dst, n, comm);

  return NULL;
}


FCSResult fcs_fmm_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  fcs_resort_resort_bytes(handle->fmm_param->fmm_resort->resort, src, dst, n, comm);
  
  return NULL;
}
