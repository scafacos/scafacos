/*
  Copyright (C) 2011-2012 Rene Halver, Olaf Lenz
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

#include <stdio.h>

#include "fcs_p3m.h"
#include "FCSCommon.h"
#include "p3m/scafacos.h"


#define P3M_CHECK_RETURN_RESULT(_h_, _f_)  do { \
  CHECK_HANDLE_RETURN_RESULT(_h_, _f_); \
  CHECK_METHOD_RETURN_RESULT(_h_, _f_, FCS_METHOD_P3M, "p3m"); \
  } while (0)

#define P3M_CHECK_RETURN_VAL(_h_, _f_, _v_)  do { \
  CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_); \
  CHECK_METHOD_RETURN_VAL(_h_, _f_, FCS_METHOD_P3M, "p3m", _v_); \
  } while (0)


/* initialization function for basic p3m parameters */
FCSResult fcs_p3m_init(FCS handle)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  handle->shift_positions = 1;

  handle->destroy = fcs_p3m_destroy;
  handle->set_r_cut = fcs_p3m_set_r_cut;
  handle->unset_r_cut = fcs_p3m_set_r_cut_tune;
  handle->get_r_cut = fcs_p3m_get_r_cut;
  handle->set_tolerance = fcs_p3m_set_tolerance;
  handle->get_tolerance = fcs_p3m_get_tolerance;
  handle->set_parameter = fcs_p3m_set_parameter;
  handle->print_parameters = fcs_p3m_print_parameters;
  handle->tune = fcs_p3m_tune;
  handle->run = fcs_p3m_run;

  ifcs_p3m_init(&handle->method_context, handle->communicator);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

/* internal p3m-specific tuning function */
FCSResult fcs_p3m_tune(FCS handle, 
		       fcs_int local_particles,
		       fcs_float *positions, fcs_float *charges)
{
  FCSResult result;

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  /* Handle periodicity */
  const fcs_int *periodicity = fcs_get_periodicity(handle);
  if (! (periodicity[0] && periodicity[1] && periodicity[2]))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, 
      "p3m requires periodic boundary conditions.");
    
  /* Handle box size */
  const fcs_float *a = fcs_get_box_a(handle);
  const fcs_float *b = fcs_get_box_b(handle);
  const fcs_float *c = fcs_get_box_c(handle);
  if (!fcs_is_orthogonal(a, b, c)) {
    if (ifcs_p3m_check_triclinic_box(a[1],a[2],b[2])) {
      result = ifcs_p3m_set_triclinic_flag(handle->method_context);
      CHECK_RESULT_RETURN(result);
    }
    else
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
          "p3m triclinic requires the box to be as follows: \n \
          the first box vector is parallel to the x axis\n \
          the second box vector is in the yz plane.");
  } else {
    if (!fcs_uses_principal_axes(a, b, c))
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
          "p3m requires the box vectors to be parallel to the principal axes.");
  }

  ifcs_p3m_set_box_a(handle->method_context, a[0]);
  ifcs_p3m_set_box_b(handle->method_context, b[1]);
  ifcs_p3m_set_box_c(handle->method_context, c[2]);

  ifcs_p3m_set_box_geometry(handle->method_context, a, b, c);

  ifcs_p3m_set_near_field_flag(handle->method_context, 
				 fcs_get_near_field_flag(handle));

  fcs_int max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  /* Effectively, tune initializes the algorithm. */
  result = ifcs_p3m_tune(handle->method_context, local_particles, max_local_particles, positions, charges);

  FCS_DEBUG_FUNC_OUTRO(__func__, result);

  return result;
}

/* internal p3m-specific run function */
FCSResult fcs_p3m_run(FCS handle, 
		      fcs_int local_particles,
		      fcs_float *positions, fcs_float *charges,
		      fcs_float *fields, fcs_float *potentials)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  fcs_p3m_tune(handle, local_particles, positions, charges);

  fcs_int max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  ifcs_p3m_run(handle->method_context, local_particles, max_local_particles, positions, charges, fields, potentials);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

/* clean-up function for p3m */
FCSResult fcs_p3m_destroy(FCS handle)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_destroy(handle->method_context);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

/******************************************************************************************************
 *
 *						Setter and Getter functions for p3m parameters
 *
 ******************************************************************************************************/


FCSResult fcs_p3m_set_r_cut(FCS handle, fcs_float r_cut) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_r_cut(handle->method_context, r_cut);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_set_r_cut_tune(FCS handle) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_r_cut_tune(handle->method_context);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_r_cut(FCS handle, fcs_float *r_cut) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_get_r_cut(handle->method_context, r_cut);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p3m_set_alpha(FCS handle, fcs_float alpha) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_alpha(handle->method_context, alpha);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p3m_set_alpha_tune(FCS handle) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_alpha_tune(handle->method_context);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_alpha(FCS handle, fcs_float *alpha) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_get_alpha(handle->method_context, alpha);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p3m_set_grid(FCS handle, fcs_int grid) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_grid(handle->method_context, grid);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_set_grid_tune(FCS handle) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_grid_tune(handle->method_context);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_grid(FCS handle, fcs_int *grid) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_get_grid(handle->method_context, grid);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p3m_set_cao(FCS handle, fcs_int cao) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_cao(handle->method_context, cao);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_set_cao_tune(FCS handle) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_cao_tune(handle->method_context);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_cao(FCS handle, fcs_int *cao) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_get_cao(handle->method_context, cao);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_require_total_energy(FCS handle, fcs_int flag) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_require_total_energy(handle->method_context, flag);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_total_energy(FCS handle, fcs_float *total_energy) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  FCSResult result = ifcs_p3m_get_total_energy(handle->method_context, total_energy);

  FCS_DEBUG_FUNC_OUTRO(__func__, result);

  return result;
}


FCSResult fcs_p3m_set_tolerance_field(FCS handle, fcs_float tolerance_field) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_tolerance_field(handle->method_context, tolerance_field);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

void fcs_p3m_set_tolerance_field_f(void *handle, fcs_float tolerance_field, fcs_int *return_value)
{
  FCSResult result = fcs_p3m_set_tolerance_field((FCS)handle, tolerance_field);

  FCS_DEBUG_MOP(printf("p3m got as accuracy: %10.4" FCS_LMOD_FLOAT "f \n", tolerance_field););

  if (result == FCS_RESULT_SUCCESS)
    *return_value = 0;
  else
    *return_value = fcs_result_get_return_code(result);
}

FCSResult fcs_p3m_set_tolerance_field_tune(FCS handle) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_tolerance_field_tune(handle->method_context);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_tolerance_field(FCS handle, fcs_float *tolerance_field) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_get_tolerance_field(handle->method_context, tolerance_field);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_near_parameters(FCS handle,
				      fcs_p3m_near_parameters_t *near_params) {

  ifcs_p3m_get_near_params(handle->method_context, &near_params->alpha, &near_params->potentialOffset);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  FCSResult result;

  if (tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
  {
    result = fcs_p3m_set_tolerance_field(handle, tolerance);

  } else
  {
    result = fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__, "Unsupported tolerance type. P3M only supports FCS_TOLERANCE_TYPE_FIELD.");
  }

  FCS_DEBUG_FUNC_OUTRO(__func__, result);

  return result;
}

FCSResult fcs_p3m_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
  FCSResult result = fcs_p3m_get_tolerance_field(handle, tolerance);

  FCS_DEBUG_FUNC_OUTRO(__func__, result);

  return result;
}

FCSResult fcs_p3m_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_r_cut",                p3m_set_r_cut,            FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_alpha",                p3m_set_alpha,            FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_grid",                 p3m_set_grid,             FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_cao",                  p3m_set_cao,              FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_require_total_energy", p3m_require_total_energy, FCS_PARSE_VAL(fcs_int));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_print_parameters(FCS handle)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  fcs_float tolerance;
  fcs_p3m_get_tolerance_field(handle, &tolerance);
  printf("p3m absolute field tolerance: %" FCS_LMOD_FLOAT "e\n", tolerance);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_set_potential_shift(FCS handle, fcs_int flag) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_p3m_set_potential_shift(handle->method_context, flag);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}

fcs_int fcs_p3m_get_potential_shift(FCS handle) {

  FCS_DEBUG_FUNC_INTRO(__func__);

  P3M_CHECK_RETURN_VAL(handle, __func__, -1);

  int r = ifcs_p3m_get_potential_shift(handle->method_context);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return r;
}
