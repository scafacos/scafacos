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


#if defined(FCS_ENABLE_DEBUG) || 0
# define PRINT_PARAM_BEGIN(_f_)     printf("%s: calling " #_f_ "(handle", func_name)
# define PRINT_PARAM_VAL(_f_, _v_)  printf( _f_, _v_)
# define PRINT_PARAM_STR(_str_)     printf("%s", (_str_))
# define PRINT_PARAM_END(_r_)       printf(") -> %p\n", _r_)
#else
# define PRINT_PARAM_BEGIN(_f_)     do {} while (0)
# define PRINT_PARAM_VAL(_f_, _v_)  do {} while (0)
# define PRINT_PARAM_STR(_str_)     do {} while (0)
# define PRINT_PARAM_END(_r_)       do {} while (0)
#endif

typedef long long fcs_long_long_t;
typedef char *fcs_p_char_t;

#define MAKE_TYPE_FUNC(_type_, _atox_, _format_) \
  static inline _type_ *parse_##_type_(char **s, _type_ *v) { \
    if (v == NULL) return NULL; \
    *v = (_type_) _atox_(*s); \
    *s = strchr(*s, ','); \
    if (*s) { **s = 0; *s += 1; } \
    PRINT_PARAM_VAL(_format_, *v); \
    return v; \
  } \
  static inline _type_ *const_##_type_(_type_ c, _type_ *v) { \
    if (v == NULL) return NULL; \
    *v = c; \
    PRINT_PARAM_VAL(_format_, *v); \
    return v; \
  }

static inline fcs_bool atob(const char *nptr)
{
  const char false_str[] = "false";
  if ((strlen(nptr) == 1 && strncmp(nptr, "0", 1) == 0) || (strlen(nptr) == strlen(false_str) && strncasecmp(nptr, false_str, strlen(false_str)))) return FCS_FALSE;
  return FCS_TRUE;
}

MAKE_TYPE_FUNC(fcs_int, atoll, "%" FCS_LMOD_INT "d")
MAKE_TYPE_FUNC(fcs_float, atof, "%" FCS_LMOD_FLOAT "f")
MAKE_TYPE_FUNC(fcs_bool, atob, "%" FCS_LMOD_INT "d")
MAKE_TYPE_FUNC(fcs_long_long_t, atoll, "%lld")
MAKE_TYPE_FUNC(fcs_p_char_t, , "%s")

#define PARSE_SEQ_MAX  3

#define PARAM_SELECTED(_str_, _param_) \
  (strcmp(param, #_param_) == 0 || strcmp(param, _str_) == 0)

#define IF_PARAM_INTRO(_str_, _param_) \
  if (PARAM_SELECTED(_str_, _param_)) { \
    FCSResult _r; \
    struct { \
      void *t; \
      fcs_int v_fcs_int[PARSE_SEQ_MAX]; \
      fcs_float v_fcs_float[PARSE_SEQ_MAX]; \
      fcs_bool v_fcs_bool[PARSE_SEQ_MAX]; \
      fcs_long_long_t v_fcs_long_long_t[PARSE_SEQ_MAX]; \
      fcs_p_char_t v_fcs_p_char_t[PARSE_SEQ_MAX]; \
    } _t; \
    char *_n=NULL; \
    PRINT_PARAM_BEGIN(_param_);

#define IF_PARAM_EXTRO() \
    PRINT_PARAM_END(_r); \
    if (_r != FCS_RESULT_SUCCESS && FCS_IS_FALSE(continue_on_errors)) return _r; \
    goto next_param; \
  }

#define PARSE_VAL(_type_) \
  _t.v_##_type_; \
  if (cur) { \
    PRINT_PARAM_STR(", "); \
    parse_##_type_(&cur, &_t.v_##_type_[0]); \
    _n = cur; cur = NULL; \
  } _type_

#define PARSE_SEQ(_type_, _n_) \
  &_t.v_##_type_; \
  if (cur) { \
    PRINT_PARAM_STR(", ["); \
    for (int _i = 0; _i < (_n_) && _i < PARSE_SEQ_MAX; ++_i) { \
      PRINT_PARAM_STR((_i == 0)?"":", "); \
      parse_##_type_(&cur, &_t.v_##_type_[_i]); \
    } \
    _n = cur; cur = NULL; \
    PRINT_PARAM_STR("]"); \
  } _type_ *

#define CONST_VAL(_type_, _c_) \
  _t.v_##_type_; \
  if (cur) { \
    PRINT_PARAM_STR(", "); \
    const_##_type_(_c_, &_t.v_##_type_[0]); \
    _n = cur; cur = NULL; \
  } _type_

#define IF_PARAM_THEN_FUNC0_GOTO_NEXT(_str_, _param_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _r = fcs_##_param_(handle); \
  IF_PARAM_EXTRO()

#define IF_PARAM_THEN_FUNC1_GOTO_NEXT(_str_, _param_, _p0_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _r = fcs_##_param_(handle, _vv0); \
  IF_PARAM_EXTRO()

#define IF_PARAM_THEN_FUNC2_GOTO_NEXT(_str_, _param_, _p0_, _p1_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _t.t = _p1_ _v1 = *_p1_ _vv1 = _v1; cur = _n; \
    _r = fcs_##_param_(handle, _vv0, _vv1); \
  IF_PARAM_EXTRO()

#define IF_PARAM_THEN_FUNC3_GOTO_NEXT(_str_, _param_, _p0_, _p1_, _p2_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _t.t = _p1_ _v1 = *_p1_ _vv1 = _v1; cur = _n; \
    _t.t = _p2_ _v2 = *_p2_ _vv2 = _v2; cur = _n; \
    _r = fcs_##_param_(handle, _vv0, _vv1, _vv2); \
  IF_PARAM_EXTRO()

#define DUMMY_REFERENCE_TO_STATIC_FUNCTIONS(_str_) \
  if (PARAM_SELECTED(_str_, "EVEN_IF_IT_MATCHES_IT_DOES_NOTHING")) { \
    parse_fcs_int(NULL, NULL); \
    const_fcs_int(0, NULL); \
    parse_fcs_float(NULL, NULL); \
    const_fcs_float(0, NULL); \
    parse_fcs_bool(NULL, NULL); \
    const_fcs_bool(0, NULL); \
    parse_fcs_long_long_t(NULL, NULL); \
    const_fcs_long_long_t(0, NULL); \
    parse_fcs_p_char_t(NULL, NULL); \
    const_fcs_p_char_t(NULL, NULL); \
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
  fcs_int params_strlen;

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

    IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_a",                   set_box_a,           PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_b",                   set_box_b,           PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_c",                   set_box_c,           PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("offset",                  set_box_origin,      PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("periodicity",             set_periodicity,     PARSE_SEQ(fcs_int, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("near_field_flag",         set_near_field_flag, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("total_particles",         set_total_particles, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("require_virial",          set_compute_virial,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("",                        set_tolerance,       PARSE_VAL(fcs_int),                                   PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_energy",        set_tolerance,       CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_ENERGY),        PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_energy_rel",    set_tolerance,       CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_ENERGY_REL),    PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_potential",     set_tolerance,       CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_POTENTIAL),     PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_potential_rel", set_tolerance,       CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_POTENTIAL_REL), PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_field",         set_tolerance,       CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_FIELD),         PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_field_rel",     set_tolerance,       CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_FIELD_REL),     PARSE_VAL(fcs_float));
#ifdef FCS_ENABLE_DIRECT
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_periodic_images",  direct_set_periodic_images,  PARSE_SEQ(fcs_int, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_cutoff",           direct_set_cutoff,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_cutoff_with_near", direct_set_cutoff_with_near, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_EWALD
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_maxkmax", ewald_set_maxkmax, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_kmax",    ewald_set_kmax,    PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_r_cut",   ewald_set_r_cut,   PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_alpha",   ewald_set_alpha,   PARSE_VAL(fcs_float));
#endif
#ifdef FCS_ENABLE_FMM
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_absrel",            fmm_set_absrel,            PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_tolerance_energy",  fmm_set_tolerance_energy,  PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_dipole_correction", fmm_set_dipole_correction, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_potential",         fmm_set_potential,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_cusp_radius",       fmm_set_cusp_radius,       PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_internal_tuning",   fmm_set_internal_tuning,   PARSE_VAL(fcs_long_long_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_maxdepth",          fmm_set_maxdepth,          PARSE_VAL(fcs_long_long_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_unroll_limit",      fmm_set_unroll_limit,      PARSE_VAL(fcs_long_long_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_balanceload",       fmm_set_balanceload,       PARSE_VAL(fcs_long_long_t));
#endif
#ifdef FCS_ENABLE_MEMD
/*    IF_PARAM_THEN_FUNC3_GOTO_NEXT("", fcs_memd_set_box_size,                  PARSE_VAL(fcs_float), PARSE_VAL(fcs_float), PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_time_step,                 PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_total_number_of_particles, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_local_number_of_particles, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_init_flag,                 PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_mesh_size_1D,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_speed_of_light,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_permittivity,              PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_temperature,               PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_bjerrum_length,            PARSE_VAL(fcs_float));*/
#endif
#ifdef FCS_ENABLE_MMM1D
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_far_switch_radius", mmm1d_set_far_switch_radius, PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_bessel_cutoff",     mmm1d_set_bessel_cutoff,     PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_maxPWerror",        mmm1d_set_maxPWerror,        PARSE_VAL(fcs_float));
#endif
#ifdef FCS_ENABLE_MMM2D
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_maxPWerror",           mmm2d_set_maxPWerror,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_far_cutoff",           mmm2d_set_far_cutoff,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("mmm2d_dielectric_contrasts", mmm2d_set_dielectric_contrasts, PARSE_VAL(fcs_float), PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_layers_per_node",      mmm2d_set_layers_per_node,      PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_skin",                 mmm2d_set_skin,                 PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("",                           mmm2d_require_total_energy,     PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_P2NFFT
    /* P2NFFT specific parameters */
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_r_cut",                p2nfft_set_r_cut,                     PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_epsI",                 p2nfft_set_epsI,                      PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_epsB",                 p2nfft_set_epsB,                      PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_c",                    p2nfft_set_c,                         PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_alpha",                p2nfft_set_alpha,                     PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_intpol_order",         p2nfft_set_interpolation_order,       PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_near",             p2nfft_set_reg_near,                  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_near_name",        p2nfft_set_reg_near_by_name,          PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_far",              p2nfft_set_reg_far,                   PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_far_name",         p2nfft_set_reg_far_by_name,           PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_p",                    p2nfft_set_p,                         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_require_virial",       p2nfft_require_virial,                PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_ignore_tolerance",     p2nfft_set_ignore_tolerance,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("p2nfft_grid",                 p2nfft_set_grid,                      PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("p2nfft_oversampled_grid",     p2nfft_set_oversampled_grid,          PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_cao",                  p2nfft_set_cao,                       PARSE_VAL(fcs_int));

    /* PNFFT specific parameters */
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("pnfft_N",                     p2nfft_set_pnfft_N,                   PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("pnfft_n",                     p2nfft_set_pnfft_n,                   PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_window_name",           p2nfft_set_pnfft_window_by_name,      PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_window",                p2nfft_set_pnfft_window,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_m",                     p2nfft_set_pnfft_m,                   PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_intpol_order",          p2nfft_set_pnfft_interpolation_order, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_phi_hat",           p2nfft_set_pnfft_pre_phi_hat,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_fg_psi",                p2nfft_set_pnfft_fg_psi,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_fft_in_place",          p2nfft_set_pnfft_fft_in_place,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_sort_nodes",            p2nfft_set_pnfft_sort_nodes,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_interlaced",            p2nfft_set_pnfft_interlaced,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_grad_ik",               p2nfft_set_pnfft_grad_ik,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_psi",               p2nfft_set_pnfft_pre_psi,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_fg_psi",            p2nfft_set_pnfft_pre_fg_psi,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_full_psi",          p2nfft_set_pnfft_pre_full_psi,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_full_fg_psi",       p2nfft_set_pnfft_pre_full_fg_psi,     PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_real_f",                p2nfft_set_pnfft_real_f,              PARSE_VAL(fcs_int));

    /* PFFT specific parameters */
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_patience",               p2nfft_set_pfft_patience,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_patience_name",          p2nfft_set_pfft_patience_by_name,      PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_preserve_input",         p2nfft_set_pfft_preserve_input,       PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_tune",                   p2nfft_set_pfft_tune,                 PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_P3M
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_r_cut",           p3m_set_r_cut,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_alpha",           p3m_set_alpha,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_grid",            p3m_set_grid,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_cao",             p3m_set_cao,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_require_total_energy",\
                                                         p3m_require_total_energy, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_PEPC
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_epsilon",           pepc_set_epsilon,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_theta",             pepc_set_theta,             PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_num_walk_threads",  pepc_set_num_walk_threads,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_dipole_correction", pepc_set_dipole_correction, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_load_balancing",    pepc_set_load_balancing,    PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_npm",               pepc_set_npm,               PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_debug_level",       pepc_set_debug_level,       PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_PP3MG
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_cells_x",        pp3mg_set_cells_x,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_cells_y",        pp3mg_set_cells_y,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_cells_z",        pp3mg_set_cells_z,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_ghosts",         pp3mg_set_ghosts,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_degree",         pp3mg_set_degree,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_max_particles",  pp3mg_set_max_particles,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_max_iterations", pp3mg_set_max_iterations, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_tol",            pp3mg_set_tol,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_distribution",   pp3mg_set_distribution,   PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_discretization", pp3mg_set_discretization, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_VMG
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_max_level",            vmg_set_max_level,            PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_max_iterations",       vmg_set_max_iterations,       PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_smoothing_steps",      vmg_set_smoothing_steps,      PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_cycle_type",           vmg_set_cycle_type,           PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_precision",            vmg_set_precision,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_near_field_cells",     vmg_set_near_field_cells,     PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_interpolation_order",  vmg_set_interpolation_order,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_discretization_order", vmg_set_discretization_order, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_WOLF
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("wolf_cutoff", wolf_set_cutoff, PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("wolf_alpha", wolf_set_alpha, PARSE_VAL(fcs_float));
#endif

    DUMMY_REFERENCE_TO_STATIC_FUNCTIONS(param);

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
