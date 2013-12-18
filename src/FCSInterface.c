/*
  Copyright (C) 2011, 2012, 2013 Rene Halver, Michael Hofmann

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

#include "FCSCommon.h"
#include "FCSInterface.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>


#define CHECK_HANDLE_RETURN_RESULT(_h_, _f_) do { \
  if (handle == FCS_NULL) \
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, _f_, "null handle supplied"); \
  } while (0)


/**
 * @brief function to check whether all initial parameters are set
 * @param handle FCS-object representing an FCS solver
 * @return whether all initial parameters are set
 */
static fcs_int fcs_init_check(FCS handle)
{
  if (handle == FCS_NULL) return 0;

  return (fcs_get_method(handle) != FCS_METHOD_NONE && fcs_get_communicator(handle) != MPI_COMM_NULL);
}


/**
 * @brief function to check whether all common parameters are set
 * @param handle FCS-object representing an FCS solver
 * @return whether all common parameters are set
 */
static fcs_int fcs_common_check(FCS handle)
{
  if (handle == FCS_NULL) return 0;

  return (handle->near_field_flag != -1)
      && (handle->box_a[0] != 0 || handle->box_a[1] != 0 || handle->box_a[2] != 0)
      && (handle->box_b[0] != 0 || handle->box_b[1] != 0 || handle->box_b[2] != 0)
      && (handle->box_c[0] != 0 || handle->box_c[1] != 0 || handle->box_c[2] != 0)
      && (handle->periodicity[0] != -1 && handle->periodicity[1] != -1 && handle->periodicity[2] != -1)
      && (handle->total_particles != -1);
}


/**
 * @brief function to check whether fcs_tune is ready to be called
 * @return whether fcs_tune is ready to be called
 */
static fcs_int fcs_tune_check(FCS handle)
{
  if (handle == FCS_NULL) return 0;

  return (fcs_init_check(handle) && fcs_common_check(handle));
}


/**
 * @brief function to check whether fcs_run is ready to be called
 * @return whether fcs_run is ready to be called
 */
static fcs_int fcs_run_check(FCS handle)
{
  if (handle == FCS_NULL) return 0;

  return (fcs_init_check(handle) && fcs_common_check(handle));
}


/**
 * initialize an FCS solver
 */
FCSResult fcs_init(FCS *new_handle, const char* method_name, MPI_Comm communicator)
{
  const char *fnc_name = "fcs_init";

  if (new_handle == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  *new_handle = FCS_NULL;

  if (method_name == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as method name");

  FCS handle = malloc(sizeof(FCS_t));

  if (handle == NULL) 
    return fcs_result_create(FCS_ERROR_ALLOC_FAILED, fnc_name, "memory allocation for FCS-object failed");

  handle->communicator = communicator;

  /* unset common parameters */
  handle->near_field_flag = -1;
  handle->box_a[0] = handle->box_a[1] = handle->box_a[2] = 0.0;
  handle->box_b[0] = handle->box_b[1] = handle->box_b[2] = 0.0;
  handle->box_c[0] = handle->box_c[1] = handle->box_c[2] = 0.0;
  handle->box_origin[0] = handle->box_origin[1] = handle->box_origin[2] = 0.0;
  handle->periodicity[0] = handle->periodicity[1] = handle->periodicity[2] = 0.0;
  handle->total_particles = -1;

#ifdef FCS_ENABLE_FMM
  handle->fmm_param = NULL;
#endif
#ifdef FCS_ENABLE_MEMD
  handle->memd_param = NULL;
#endif
#ifdef FCS_ENABLE_P2NFFT
  handle->p2nfft_param = NULL;
#endif
#ifdef FCS_ENABLE_PEPC
  handle->pepc_param = NULL;
#endif
#ifdef FCS_ENABLE_PP3MG
  handle->pp3mg_param = NULL;
#endif
#ifdef FCS_ENABLE_VMG
  handle->vmg_param = NULL;
#endif
#ifdef FCS_ENABLE_WOLF
  handle->wolf_param = NULL;
#endif

  handle->method_context = NULL;

  handle->set_max_particle_move = NULL;
  handle->set_resort = NULL;
  handle->get_resort = NULL;
  handle->get_resort_availability = NULL;
  handle->get_resort_particles = NULL;
  handle->resort_ints = NULL;
  handle->resort_floats = NULL;
  handle->resort_bytes = NULL;

  *new_handle = handle;

  /* call the method-specific init functions */
#define METHOD_COMPARE_INIT_AND_RETURN(_n_, _id_) \
  if (strcmp(method_name, #_n_) == 0) { \
    handle->method = _id_; \
    strncpy(handle->method_name, method_name, FCS_MAX_METHOD_NAME_LENGTH); \
    return fcs_##_n_##_init(handle); \
  }

#ifdef FCS_ENABLE_DIRECT
  METHOD_COMPARE_INIT_AND_RETURN(direct, FCS_METHOD_DIRECT);
#endif
#ifdef FCS_ENABLE_EWALD
  METHOD_COMPARE_INIT_AND_RETURN(ewald, FCS_METHOD_EWALD);
#endif
#ifdef FCS_ENABLE_FMM
  METHOD_COMPARE_INIT_AND_RETURN(fmm, FCS_METHOD_FMM);
#endif
#ifdef FCS_ENABLE_MEMD
  METHOD_COMPARE_INIT_AND_RETURN(memd, FCS_METHOD_MEMD);
#endif
#ifdef FCS_ENABLE_MMM1D
  METHOD_COMPARE_INIT_AND_RETURN(mmm1d, FCS_METHOD_MMM1D);
#endif
#ifdef FCS_ENABLE_MMM2D
  METHOD_COMPARE_INIT_AND_RETURN(mmm2d, FCS_METHOD_MMM2D);
#endif
#ifdef FCS_ENABLE_P2NFFT
  METHOD_COMPARE_INIT_AND_RETURN(p2nfft, FCS_METHOD_P2NFFT);
#endif
#ifdef FCS_ENABLE_P3M
  METHOD_COMPARE_INIT_AND_RETURN(p3m, FCS_METHOD_P3M);
#endif
#ifdef FCS_ENABLE_PEPC
  METHOD_COMPARE_INIT_AND_RETURN(pepc, FCS_METHOD_PEPC);
#endif
#ifdef FCS_ENABLE_PP3MG
  METHOD_COMPARE_INIT_AND_RETURN(pp3mg, FCS_METHOD_PP3MG);
#endif
#ifdef FCS_ENABLE_VMG
  METHOD_COMPARE_INIT_AND_RETURN(vmg, FCS_METHOD_VMG);
#endif
#ifdef FCS_ENABLE_WOLF
  METHOD_COMPARE_INIT_AND_RETURN(wolf, FCS_METHOD_WOLF);
#endif

  /* none of the known methods was chosen */
  return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
}


/**
 * destroy an FCS solver
 */
FCSResult fcs_destroy(FCS handle)
{
  const char *fnc_name = "fcs_destroy";
  FCSResult result;

  if (handle == FCS_NULL) return FCS_RESULT_SUCCESS;

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_DIRECT
    case FCS_METHOD_DIRECT:
      result = fcs_direct_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_EWALD
    case FCS_METHOD_EWALD:
      result = fcs_ewald_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      result = fcs_fmm_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_METHOD_MEMD:
      result = fcs_memd_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_METHOD_MMM1D:
      result = fcs_mmm1d_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_MMM2D
    case FCS_METHOD_MMM2D:
      result = fcs_mmm2d_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_METHOD_PEPC:
      result = fcs_pepc_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      result = fcs_p2nfft_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      result = fcs_p3m_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_METHOD_PP3MG:
      result = fcs_pp3mg_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_METHOD_VMG:
      result = fcs_vmg_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_METHOD_WOLF:
      result = fcs_wolf_destroy(handle);
      if (result != FCS_RESULT_SUCCESS) return result;
      break;
#endif
    default:
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
  }

  free(handle);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the numerical identifier of the solver method
 */
fcs_int fcs_get_method(FCS handle)
{
  if (handle == FCS_NULL) return FCS_METHOD_NONE;

  return handle->method;
}


/**
 * return the name of the solver method
 */
const char *fcs_get_method_name(FCS handle)
{
  if (handle == FCS_NULL) return "none";

  return handle->method_name;
}


/**
 * return the MPI communicator used for the parallel execution
 */
MPI_Comm fcs_get_communicator(FCS handle)
{
  const char *fnc_name = "fcs_get_communicator";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning MPI_COMM_NULL", fnc_name);
    return MPI_COMM_NULL;
  }

  return handle->communicator;
}


/**
 * set all obligatory parameters for an FCS solver
 */
FCSResult fcs_set_common(FCS handle, fcs_int near_field_flag, const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c, const fcs_float *box_origin, const fcs_int *periodicity, fcs_int total_particles)
{
  const char *fnc_name = "fcs_set_common";
  FCSResult result;

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  result = fcs_set_near_field_flag(handle, near_field_flag);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_set_box_a(handle, box_a);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_set_box_b(handle, box_b);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_set_box_c(handle, box_c);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_set_box_origin(handle, box_origin);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_set_periodicity(handle, periodicity);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_set_total_particles(handle, total_particles);
  if (result != FCS_RESULT_SUCCESS) return result;

  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * set the dimensions of the system
 */
FCSResult fcs_set_dimensions(FCS handle, fcs_int dim)
{
  const char *fnc_name = "fcs_set_dimensions";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (dim > 3 || dim < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "dimensions must be between 1 and 3");

  handle->dimensions = dim;

  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the dimensions of the system
 */
fcs_int fcs_get_dimensions(FCS handle)
{
  const char *fnc_name = "fcs_get_dimensions";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning -1", fnc_name);
    return -1;
  }

  return handle->dimensions;
}


/**
 * set the near-field flag
 */
FCSResult fcs_set_near_field_flag(FCS handle, fcs_int near_field_flag)
{
  const char *fnc_name = "fcs_set_near_field_flag";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  handle->near_field_flag = near_field_flag;

  fcs_set_values_changed(handle,1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the near-field flag
 */
fcs_int fcs_get_near_field_flag(FCS handle)
{
  const char *fnc_name = "fcs_get_near_field_flag";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning -1", fnc_name);
    return -1;
  }

  return handle->near_field_flag;
}


/**
 * set the first base vector of the system box
 */
FCSResult fcs_set_box_a(FCS handle, const fcs_float *box_a)
{
  const char *fnc_name = "fcs_set_box_a";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (box_a == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as argument");

  handle->box_a[0] = box_a[0];
  handle->box_a[1] = box_a[1];
  handle->box_a[2] = box_a[2];
  
  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the first base vector of the system box
 */
const fcs_float* fcs_get_box_a(FCS handle)
{
  const char *fnc_name = "fcs_get_box_a";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning NULL", fnc_name);
    return NULL;
  }

  return handle->box_a;
}


/**
 * set the second base vector of the system box
 */
FCSResult fcs_set_box_b(FCS handle, const fcs_float *box_b)
{
  const char *fnc_name = "fcs_set_box_b";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (box_b == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as argument");

  handle->box_b[0] = box_b[0];
  handle->box_b[1] = box_b[1];
  handle->box_b[2] = box_b[2];
  
  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the second base vector of the system box
 */
const fcs_float* fcs_get_box_b(FCS handle)
{
  const char *fnc_name = "fcs_get_box_b";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning NULL", fnc_name);
    return NULL;
  }

  return handle->box_b;
}


/**
 * set the third base vector of the system box
 */
FCSResult fcs_set_box_c(FCS handle, const fcs_float *box_c)
{
  const char *fnc_name = "fcs_set_box_c";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (box_c == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as argument");

  handle->box_c[0] = box_c[0];
  handle->box_c[1] = box_c[1];
  handle->box_c[2] = box_c[2];
  
  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the third base vector of the system box
 */
const fcs_float* fcs_get_box_c(FCS handle)
{
  const char *fnc_name = "fcs_get_box_c";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning NULL", fnc_name);
    return NULL;
  }

  return handle->box_c;
}


/**
 * set the origin vector of the system box
 */
FCSResult fcs_set_box_origin(FCS handle, const fcs_float *box_origin)
{
  const char *fnc_name = "fcs_set_box_origin";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (box_origin == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as argument");

  handle->box_origin[0] = box_origin[0];
  handle->box_origin[1] = box_origin[1];
  handle->box_origin[2] = box_origin[2];
  
  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the origin vector of the system box
 */
const fcs_float *fcs_get_box_origin(FCS handle)
{
  const char *fnc_name = "fcs_get_box_origin";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning NULL", fnc_name);
    return NULL;
  }

  return handle->box_origin;
}


/**
 * set the periodicity of the system
 */
FCSResult fcs_set_periodicity(FCS handle, const fcs_int *periodicity)
{
  const char *fnc_name = "fcs_set_periodicity";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (periodicity == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as argument");

  handle->periodicity[0] = periodicity[0];
  handle->periodicity[1] = periodicity[1];
  handle->periodicity[2] = periodicity[2];

  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the periodicity of the system
 */
const fcs_int* fcs_get_periodicity(FCS handle)
{
  const char *fnc_name = "fcs_get_periodicity";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning NULL", fnc_name);
    return NULL;
  }

  return handle->periodicity;
}


/**
 * set the total number of particles in the system
 */
FCSResult fcs_set_total_particles(FCS handle, fcs_int total_particles)
{
  const char *fnc_name = "fcs_set_total_particles";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (total_particles < 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "total number of particles must be at least 1");
  
  handle->total_particles = total_particles;

  fcs_set_values_changed(handle, 1);

  return FCS_RESULT_SUCCESS;
}


/**
 * return the total number of particles in the system
 */
fcs_int fcs_get_total_particles(FCS handle)
{
  const char *fnc_name = "fcs_get_total_particles";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning -1", fnc_name);
    return -1;
  }

  return handle->total_particles;
}


/**
 * set the method context information
 */
FCSResult fcs_set_method_context(FCS handle, void *method_context)
{
  const char *fnc_name = "fcs_set_method_context";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  handle->method_context = method_context;

  return FCS_RESULT_SUCCESS;
}


/**
 * return the method context information
 */
void *fcs_get_method_context(FCS handle)
{
  const char *fnc_name = "fcs_get_method_context";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning NULL", fnc_name);
    return NULL;
  }

  return handle->method_context;
}


/**
 * set whether parameter values of the FCS solver have changed
 */
FCSResult fcs_set_values_changed(FCS handle, fcs_int values_changed)
{
  const char *fnc_name = "fcs_set_values_changed";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  handle->values_changed = values_changed;

  return FCS_RESULT_SUCCESS;
}


/**
 * return whether parameter values of the FCS solver have changed
 */
fcs_int fcs_get_values_changed(FCS handle)
{
  const char *fnc_name = "fcs_get_box_values_changed";

  if (handle == FCS_NULL)
  {
    fprintf(stderr, "%s: null handle supplied, returning -1", fnc_name);
    return -1;
  }

  return handle->values_changed;
}


/**
 * set the error tolerance of the FCS solver
 */
FCSResult fcs_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance)
{
  const char *fnc_name = "fcs_set_tolerance";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_EWALD
    case FCS_METHOD_EWALD:
      if (tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
      {
        fcs_ewald_set_tolerance_field(handle, tolerance);
        return FCS_RESULT_SUCCESS;
      } else return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "Unsupported tolerance type. EWALD only supports FCS_TOLERANCE_TYPE_FIELD.");
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      if (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY)
      {
        fcs_fmm_set_absrel(handle, FCS_FMM_CUSTOM_ABSOLUTE);
        fcs_fmm_set_tolerance_energy(handle, tolerance);
        return FCS_RESULT_SUCCESS;
      } else if (tolerance_type == FCS_TOLERANCE_TYPE_ENERGY_REL)
      {
        fcs_fmm_set_absrel(handle, FCS_FMM_CUSTOM_RELATIVE);
        fcs_fmm_set_tolerance_energy(handle, tolerance);
        return FCS_RESULT_SUCCESS;
      } else return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "Unsupported tolerance type. FMM only supports FCS_TOLERANCE_TYPE_ENERGY and FCS_TOLERANCE_TYPE_ENERGY_REL.");
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_set_tolerance(handle, tolerance_type, tolerance);
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      if (tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
      {
        fcs_p3m_set_tolerance_field(handle, tolerance);
        return FCS_RESULT_SUCCESS;
      } else return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "Unsupported tolerance type. P3M only supports FCS_TOLERANCE_TYPE_FIELD.");
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Setting tolerance not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * return the error tolerance of the FCS solver
 */
FCSResult fcs_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance)
{
  const char *fnc_name = "fcs_get_tolerance";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
  *tolerance = -1.0; 

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_EWALD
    case FCS_METHOD_EWALD:
      *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
      return fcs_ewald_get_tolerance_field(handle, tolerance);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_get_tolerance(handle, tolerance_type, tolerance);
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
      fcs_p3m_get_tolerance_field(handle, tolerance);
      return FCS_RESULT_SUCCESS;
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Return tolerance not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * set a user-defined cutoff radius for the near-field
 */
FCSResult fcs_set_r_cut(FCS handle, fcs_float r_cut)
{
  const char *fnc_name = "fcs_set_r_cut";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      return fcs_p3m_set_r_cut(handle, r_cut);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_set_r_cut(handle, r_cut);
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Setting a user-defined cutoff radius for the near-field not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * disable a user-defined cutoff radius for the near-field
 */
FCSResult fcs_unset_r_cut(FCS handle)
{
  const char *fnc_name = "fcs_unset_r_cut";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      return fcs_p3m_set_r_cut_tune(handle);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_set_r_cut_tune(handle);
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Disabling a user-defined cutoff radius for the near-field not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * return the user-defined cutoff radius for the near-field
 */
FCSResult fcs_get_r_cut(FCS handle, fcs_float *r_cut)
{
  const char *fnc_name = "fcs_get_r_cut";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      return fcs_p3m_get_r_cut(handle, r_cut);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_get_r_cut(handle, r_cut);
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Returning a user-defined cutoff radius for the near-field not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * set the parameters of the FCS solver based on a parameter string
 */
FCSResult fcs_set_parameters(FCS handle, const char *parameters, fcs_bool continue_on_errors)
{
  const char *fnc_name = "fcs_set_parameters";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  FCSResult r = FCS_RESULT_SUCCESS;

  char *cur;
  char *params, *param;
  fcs_int params_strlen, matched;

  params_strlen = strlen(parameters) + 1;
  params = malloc(params_strlen * sizeof(char));
  strncpy(params, parameters, params_strlen);

  cur = params;

  while (cur)
  {
    param = cur;

    cur = strchr(cur, ',');

    if (cur)
    {
      *cur = 0;
      ++cur;
    }

/*    printf("param: %s\n", param);
    printf("cur: %s\n", cur);*/

    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_a",                   set_box_a,           FCS_PARSE_SEQ(fcs_float, 3));
    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_b",                   set_box_b,           FCS_PARSE_SEQ(fcs_float, 3));
    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_c",                   set_box_c,           FCS_PARSE_SEQ(fcs_float, 3));
    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("offset",                  set_box_origin,      FCS_PARSE_SEQ(fcs_float, 3));
    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("periodicity",             set_periodicity,     FCS_PARSE_SEQ(fcs_int, 3));
    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("near_field_flag",         set_near_field_flag, FCS_PARSE_VAL(fcs_int));
    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("total_particles",         set_total_particles, FCS_PARSE_VAL(fcs_int));
    FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("require_virial",          set_compute_virial,  FCS_PARSE_VAL(fcs_int));
    FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT("",                        set_tolerance,       FCS_PARSE_VAL(fcs_int),                                     FCS_PARSE_VAL(fcs_float));
    FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_energy",        set_tolerance,       FCS_PARSE_CONST(fcs_int, FCS_TOLERANCE_TYPE_ENERGY),        FCS_PARSE_VAL(fcs_float));
    FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_energy_rel",    set_tolerance,       FCS_PARSE_CONST(fcs_int, FCS_TOLERANCE_TYPE_ENERGY_REL),    FCS_PARSE_VAL(fcs_float));
    FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_potential",     set_tolerance,       FCS_PARSE_CONST(fcs_int, FCS_TOLERANCE_TYPE_POTENTIAL),     FCS_PARSE_VAL(fcs_float));
    FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_potential_rel", set_tolerance,       FCS_PARSE_CONST(fcs_int, FCS_TOLERANCE_TYPE_POTENTIAL_REL), FCS_PARSE_VAL(fcs_float));
    FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_field",         set_tolerance,       FCS_PARSE_CONST(fcs_int, FCS_TOLERANCE_TYPE_FIELD),         FCS_PARSE_VAL(fcs_float));
    FCS_PARSE_IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_field_rel",     set_tolerance,       FCS_PARSE_CONST(fcs_int, FCS_TOLERANCE_TYPE_FIELD_REL),     FCS_PARSE_VAL(fcs_float));
#ifdef FCS_ENABLE_DIRECT
    r = fcs_direct_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_EWALD
    r = fcs_ewald_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_FMM
    r = fcs_fmm_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_MEMD
    r = fcs_memd_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_MMM1D
    r = fcs_mmm1d_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_MMM2D
    r = fcs_mmm2d_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_P2NFFT
    r = fcs_p2nfft_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_P3M
    r = fcs_p3m_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_PEPC
    r = fcs_pepc_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_PP3MG
    r = fcs_pp3mg_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_VMG
    r = fcs_vmg_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif
#ifdef FCS_ENABLE_WOLF
    r = fcs_wolf_set_parameter(handle, continue_on_errors, &param, &cur, &matched);
    if (matched) goto next_param;
#endif

    FCS_PARSE_DUMMY_REFERENCE_TO_STATIC_FUNCTIONS(param);

    if (r == FCS_RESULT_SUCCESS)
      r = fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "interface (parser): error in parameter string at '%s'!", param); 

    if (FCS_IS_FALSE(continue_on_errors)) break;

next_param:
    ;
  }

  free(params);

  return r;
}


/**
 * print the parameters of an FCS solver to stdout
 */
FCSResult fcs_print_parameters(FCS handle)
{
  const char *fnc_name = "fcs_print_parameters";
  FCSResult result;

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  printf("chosen method: %s\n", fcs_get_method_name(handle));

  printf("near field computations done by solver: %c\n", (fcs_get_near_field_flag(handle)?'T':'F'));
  printf("box vectors: [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f]\n",
    fcs_get_box_a(handle)[0], fcs_get_box_a(handle)[1], fcs_get_box_a(handle)[2],
    fcs_get_box_b(handle)[0], fcs_get_box_b(handle)[1], fcs_get_box_b(handle)[2],
    fcs_get_box_c(handle)[0], fcs_get_box_c(handle)[1], fcs_get_box_c(handle)[2]);
  printf("box origin: [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f]\n",
    fcs_get_box_origin(handle)[0], fcs_get_box_origin(handle)[1], fcs_get_box_origin(handle)[2]);
  printf("periodicity: %c %c %c\n", ((fcs_get_periodicity(handle)[0] == 1)?'T':'F'), ((fcs_get_periodicity(handle)[1] == 1)?'T':'F'),((fcs_get_periodicity(handle)[2] == 1)?'T':'F'));
  printf("total particles: %" FCS_LMOD_INT "d\n", fcs_get_total_particles(handle));
  printf("------------------------");
  printf("solver specific data:\n");

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_DIRECT
    case FCS_METHOD_DIRECT:
      result = fcs_direct_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_EWALD
    case FCS_METHOD_EWALD:
      result = fcs_ewald_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      result = fcs_fmm_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_METHOD_MEMD:
      result = fcs_memd_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_METHOD_MMM1D:
      result = fcs_mmm1d_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_MMM2D
  case FCS_METHOD_MMM2D:
      result = fcs_mmm2d_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_METHOD_PEPC:
      result = fcs_pepc_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      result = fcs_p2nfft_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      result = fcs_p3m_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_METHOD_PP3MG:
      result = fcs_pp3mg_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_METHOD_VMG:
      result = fcs_vmg_print_parameters(handle);
      break;
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_METHOD_WOLF:
      result = fcs_wolf_print_parameters(handle);
      break;
#endif
  }

  if (result != FCS_RESULT_SUCCESS) fcs_result_print_result(result);

  return FCS_RESULT_SUCCESS;
}


/**
 * tune method specific parameters depending on the particles
 */
FCSResult fcs_tune(FCS handle, fcs_int local_particles, fcs_int max_local_particles,
  fcs_float *positions, fcs_float *charges)
{
  const char *fnc_name = "fcs_tune";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (local_particles < 0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "number of local particles must be non negative");
  if (max_local_particles < 0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "maximum number of local particles must be non negative");
  if (local_particles > max_local_particles)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "maximum number of local particles should be greater or equal to current number of particles");

  if (!fcs_init_check(handle) || !fcs_tune_check(handle))
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, fnc_name, "not all needed data has been inserted into the given handle");

  fcs_set_values_changed(handle,0);
    
  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_DIRECT
    case FCS_METHOD_DIRECT:
      return fcs_direct_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_EWALD
    case FCS_METHOD_EWALD:
      return fcs_ewald_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      return fcs_fmm_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_METHOD_MEMD:
      return fcs_memd_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_METHOD_MMM1D:
      return fcs_mmm1d_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_MMM2D
    case FCS_METHOD_MMM2D:
      return fcs_mmm2d_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      return fcs_p3m_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_METHOD_PEPC:
      return fcs_pepc_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_METHOD_PP3MG:
      return fcs_pp3mg_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_METHOD_VMG:
      return fcs_vmg_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_METHOD_WOLF:
      return fcs_wolf_tune(handle, local_particles, max_local_particles, positions, charges);
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Tuning tune method specific parameters not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * run the solver method
 */
FCSResult fcs_run(FCS handle, fcs_int local_particles, fcs_int max_local_particles,
  fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  const char *fnc_name = "fcs_run";
  FCSResult result;

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (local_particles < 0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "number of local particles must be non negative");
  if (max_local_particles < 0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "maximum number of local particles must be non negative");
  if (local_particles > max_local_particles)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "maximum number of local particles should be greater or equal to current number of particles");

  if (fcs_get_values_changed(handle))
  {
    result = fcs_tune(handle, local_particles, max_local_particles, positions, charges);
    if (result != FCS_RESULT_SUCCESS) return result;
  }

  if (!fcs_init_check(handle) || !fcs_run_check(handle))
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, fnc_name, "not all needed data has been inserted into the given handle");

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_DIRECT
    case FCS_METHOD_DIRECT:
      return fcs_direct_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_EWALD
    case FCS_METHOD_EWALD:
      return fcs_ewald_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      return fcs_fmm_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_METHOD_MEMD:
      return fcs_memd_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_METHOD_MMM1D:
      return fcs_mmm1d_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_MMM2D
    case FCS_METHOD_MMM2D:
      return fcs_mmm2d_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_METHOD_PEPC:
      return fcs_pepc_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      return fcs_p3m_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_METHOD_PP3MG:
      return fcs_pp3mg_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_METHOD_VMG:
      return fcs_vmg_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_METHOD_WOLF:
      return fcs_wolf_run(handle, local_particles, max_local_particles, positions, charges, field, potentials);
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Running solver method '%s' not implemented", fcs_get_method_name(handle));
}


/**
 * compute the correction to the field and total energy 
 */
FCSResult fcs_compute_dipole_correction(FCS handle, fcs_int local_particles,
  fcs_float* positions, fcs_float *charges, fcs_float epsilon,
  fcs_float *field_correction, fcs_float *energy_correction)
{
  const char *fnc_name = "fcs_compute_dipole_correction";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  /* Local dipole moment */
  fcs_float local_dipole_moment[3] = {0.0, 0.0, 0.0};
  /* Global dipole moment */
  fcs_float dipole_moment[3];

  fcs_int pid;
  fcs_int dim;

  if (fcs_float_is_zero(epsilon) || epsilon > 0.0) {
    /* Compute the global dipole moment */
    for (pid = 0; pid < local_particles; pid++)
      for (dim = 0; dim < 3; dim++)
        local_dipole_moment[dim] += charges[pid]*positions[pid*3+dim];
    MPI_Allreduce(local_dipole_moment, dipole_moment, 3, FCS_MPI_FLOAT,
      MPI_SUM, handle->communicator);

    const fcs_float *a = fcs_get_box_a(handle);
    const fcs_float *b = fcs_get_box_b(handle);
    const fcs_float *c = fcs_get_box_c(handle);
    
    /* Volume of the parallelepiped */
    fcs_float volume = 
        a[0] * (b[1]*c[2] - b[2]*c[1]) 
      + a[1] * (b[2]*c[0] - b[0]*c[2])
      + a[2] * (b[0]*c[1] - b[1]*c[0]);

    fcs_float pref = 4.0*3.14159265358979323846264338328 
      / (3.0*volume*(epsilon + 1.0));

    if (energy_correction)
      *energy_correction = 0.5*pref*(dipole_moment[0]*dipole_moment[0]
        + dipole_moment[1]*dipole_moment[1]
        + dipole_moment[2]*dipole_moment[2]);
    
    if (field_correction) {
      field_correction[0] = -pref*dipole_moment[0];
      field_correction[1] = -pref*dipole_moment[1];
      field_correction[2] = -pref*dipole_moment[2];
    }
  } else {
    /* metallic BC (epsilon=+infty) */
    if (energy_correction)
      *energy_correction = 0.0;
    if (field_correction) {
      field_correction[0] = 0.0;
      field_correction[1] = 0.0;
      field_correction[2] = 0.0;
    }
  }

  return FCS_RESULT_SUCCESS;
}


/**
 * return whether the solver method supports the delegation of near-field computations
 */
FCSResult fcs_get_near_field_delegation(FCS handle, fcs_int *near_field_delegation)
{
  const char *fnc_name = "fcs_get_near_field_delegation";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  *near_field_delegation = 0;

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      *near_field_delegation = 1;
      break;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      *near_field_delegation = 1;
      break;
#endif
  }

  return FCS_RESULT_SUCCESS;
}


/**
 * compute the near-field components of the potential and the field
 */
FCSResult fcs_compute_near(FCS handle, fcs_float dist, fcs_float *potential, fcs_float *field)
{
  const char *fnc_name = "fcs_compute_near";

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      fcs_p2nfft_compute_near(handle, dist, potential, field);
      return FCS_RESULT_SUCCESS;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      {
        fcs_p3m_near_parameters_t params;
        fcs_p3m_get_near_parameters(handle, &params);
        fcs_p3m_compute_near(params, dist, potential, field);
      }
      return FCS_RESULT_SUCCESS;
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Computing the near-field components of the potential and the field not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * compute the near-field component of the potential
 */
FCSResult fcs_compute_near_potential(FCS handle, fcs_float dist, fcs_float *potential)
{
  const char *fnc_name = "fcs_compute_near_potential";

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      *potential = fcs_p2nfft_compute_near_potential(handle, dist);
      return FCS_RESULT_SUCCESS;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      {
        fcs_p3m_near_parameters_t params;
        fcs_p3m_get_near_parameters(handle, &params);
        *potential = fcs_p3m_compute_near_potential(params, dist);
      }
      return FCS_RESULT_SUCCESS;
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Computing the near-field component of the potential not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * compute the near-field component of the field
 */
FCSResult fcs_compute_near_field(FCS handle, fcs_float dist, fcs_float *field)
{
  const char *fnc_name = "fcs_compute_near_field";

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      *field = fcs_p2nfft_compute_near_field(handle, dist);
      return FCS_RESULT_SUCCESS;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      {
        fcs_p3m_near_parameters_t params;
        fcs_p3m_get_near_parameters(handle, &params);
        *field = fcs_p3m_compute_near_field(params, dist);
      }
      return FCS_RESULT_SUCCESS;
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Computing the near-field component of the field not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * set whether the virial should be computed
 */
FCSResult fcs_set_compute_virial(FCS handle, fcs_int compute_virial)
{
  const char *fnc_name = "fcs_require_virial";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (compute_virial != 0 && compute_virial != 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "parameter compute_virial must be 0 or 1");

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_DIRECT
    case FCS_METHOD_DIRECT:
      return fcs_direct_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_METHOD_PEPC:
      return fcs_pepc_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      return fcs_fmm_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      return fcs_p3m_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_METHOD_MMM1D:
      return fcs_mmm1d_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_MMM2D
    case FCS_METHOD_MMM2D:
      return fcs_mmm2d_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_METHOD_MEMD:
      return fcs_memd_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_METHOD_PP3MG:
      return fcs_pp3mg_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_METHOD_VMG:
      return fcs_vmg_require_virial(handle, compute_virial);
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_METHOD_WOLF:
      return fcs_wolf_require_virial(handle, compute_virial);
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Setting whether the virial should be computed not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * return whether the virial should be computed
 */
FCSResult fcs_get_compute_virial(FCS handle, fcs_int *compute_virial)
{
  const char *fnc_name = "fcs_get_compute_virial";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (compute_virial == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as argument");

  switch (fcs_get_method(handle))
  {
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Returning whether the virial should be computed not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * return the comuputed virial
 */
FCSResult fcs_get_virial(FCS handle, fcs_float *virial)
{
  const char *fnc_name = "fcs_get_virial";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (virial == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as argument");

  switch (fcs_get_method(handle))
  {
#ifdef FCS_ENABLE_DIRECT
    case FCS_METHOD_DIRECT:
      return fcs_direct_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_METHOD_PEPC:
      return fcs_pepc_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METHOD_FMM:
      return fcs_fmm_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METHOD_P3M:
      return fcs_p3m_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_METHOD_MMM1D:
      return fcs_mmm1d_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_MMM2D
    case FCS_METHOD_MMM2D:
      return fcs_mmm2d_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_METHOD_MEMD:
      return fcs_memd_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_METHOD_P2NFFT:
      return fcs_p2nfft_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_METHOD_PP3MG:
      return fcs_pp3mg_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_METHOD_VMG:
      return fcs_vmg_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_METHOD_WOLF:
      return fcs_wolf_get_virial(handle, virial);
#endif
  }

  return fcs_result_create(FCS_ERROR_NOT_IMPLEMENTED, fnc_name, "Returning the computed virial not implemented for solver method '%s'", fcs_get_method_name(handle));
}


/**
 * set the maximum distance the particles have moved since the call of ::fcs_run
 */
FCSResult fcs_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  const char *fnc_name = "fcs_set_max_particle_move";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (handle->set_max_particle_move == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "max. particle move not supported");

  return handle->set_max_particle_move(handle, max_particle_move);
}


/**
 * set whether resort support is requested
 */
FCSResult fcs_set_resort(FCS handle, fcs_int resort)
{
  const char *fnc_name = "fcs_set_resort";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (handle->set_resort == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->set_resort(handle, resort);
}


/**
 * return whether resort support is requested
 */
FCSResult fcs_get_resort(FCS handle, fcs_int *resort)
{
  const char* fnc_name = "fcs_get_resort";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (handle->get_resort == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->get_resort(handle, resort);
}


/**
 * return whether resort support is available
 */
FCSResult fcs_get_resort_availability(FCS handle, fcs_int *availability)
{
  const char *fnc_name = "fcs_get_resort_availability";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  *availability = 0;

  if (handle->get_resort_availability == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->get_resort_availability(handle, availability);
}


/**
 * return the new local number of particles
 */
FCSResult fcs_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  const char *fnc_name = "fcs_get_resort_particles";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (handle->get_resort_particles == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->get_resort_particles(handle, resort_particles);
}


/**
 * sort additional integer particle data
 */
FCSResult fcs_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n)
{
  const char *fnc_name = "fcs_resort_ints";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (handle->resort_ints == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->resort_ints(handle, src, dst, n, fcs_get_communicator(handle));
}


/**
 * sort additional float particle data
 */
FCSResult fcs_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n)
{
  const char* fnc_name = "fcs_resort_floats";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (handle->resort_floats == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->resort_floats(handle, src, dst, n, fcs_get_communicator(handle));
}


/**
 * sort additional byte particle data
 */
FCSResult fcs_resort_bytes(FCS handle, void *src, void *dst, fcs_int n)
{
  const char *fnc_name = "fcs_resort_bytes";

  CHECK_HANDLE_RETURN_RESULT(handle, fnc_name);

  if (handle->resort_bytes == NULL)
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->resort_bytes(handle, src, dst, n, fcs_get_communicator(handle));
}


/**
 * Fortran wrapper function ot initialize an FCS solver method
 */
FCSResult fcs_init_f(FCS *handle, const char *method_name, MPI_Fint communicator)
{
  MPI_Comm c_comm = MPI_Comm_f2c(communicator);
  return fcs_init(handle, method_name, c_comm);
}
