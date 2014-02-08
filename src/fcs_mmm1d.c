/*
  Copyright (C) 2011-2012 Pedro Sanchez

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

#include "fcs_mmm1d.h"
#include "mmm1d/mmm1d.h"
#include "mmm1d/types.h"

/* forward declarations */
static FCSResult fcs_mmm1d_check(FCS handle, const char* fnc_name);

/* initialization function for basic p3m parameters */
FCSResult fcs_mmm1d_init(FCS handle)
{
  char* fnc_name = "fcs_mmm1d_init";
  FCSResult result;
  
  result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  
  handle->tune = fcs_mmm1d_tune;
  handle->run = fcs_mmm1d_run;

  mmm1d_init(&handle->method_context, handle->communicator);
  return NULL;
}

FCSResult fcs_mmm1d_tune(FCS handle, 
		       fcs_int local_particles,
		       fcs_float *positions, fcs_float *charges)
{
  char* fnc_name = "fcs_mmm1d_tune";
  FCSResult result;
  
  result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  
  /* Check box periodicity */
  const fcs_int *periodicity = fcs_get_periodicity(handle);
  if (periodicity[0] != 0 ||
      periodicity[1] != 0 ||
      periodicity[2] != 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "mmm1d requires z-axis periodic boundary.");
  
  /* Check box shape */
  const fcs_float *a = fcs_get_box_a(handle);
  const fcs_float *b = fcs_get_box_b(handle);
  const fcs_float *c = fcs_get_box_c(handle);
  if (!fcs_is_cubic(a, b, c))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "mmm1d requires a cubic box.");
  
  mmm1d_set_box_a(handle->method_context, a[0]);
  mmm1d_set_box_b(handle->method_context, b[1]);
  mmm1d_set_box_c(handle->method_context, c[2]);
  
  /* Effectively, tune initializes the algorithm */
  result = mmm1d_tune(handle->method_context, 
        local_particles,
        positions, charges);
  
  return result;
}

/* internal mmm1d-specific run function */
FCSResult fcs_mmm1d_run(FCS handle, 
		      fcs_int local_particles,
		      fcs_float *positions, fcs_float *charges,
		      fcs_float *fields, fcs_float *potentials)
{
//   char* fnc_name = "fcs_mmm1d_run";
//   FCSResult result;

  fcs_mmm1d_tune(handle, local_particles, positions, charges);
  
  fcs_int max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  mmm1d_run(handle->method_context,
	  local_particles, max_local_particles, positions, charges, fields, potentials);
  return NULL;
}

/* clean-up function for mmm1d */
FCSResult fcs_mmm1d_destroy(FCS handle)
{
  mmm1d_destroy(handle->method_context);
  return NULL;
}

/* method to check if mmm1d parameters are entered into FCS */
FCSResult fcs_mmm1d_check(FCS handle, const char* fnc_name) {
  if (handle == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_METHOD_MMM1D)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"mmm1d\".");
  
  return NULL;
}

/******************************************************************************************************
 *
 *            Setter and Getter functions for mmm1d parameters
 *
 ******************************************************************************************************/
FCSResult fcs_mmm1d_set_far_switch_radius(FCS handle, fcs_float rad2) {
  const char *fnc_name = "fcs_mmm1d_set_far_switch_radius";
  
  FCSResult result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  mmm1d_set_far_switch_radius_2(handle->method_context, rad2*rad2);
  return NULL;
}

FCSResult fcs_mmm1d_get_far_switch_radius(FCS handle, fcs_float *rad2) {
  const char *fnc_name = "fcs_mmm1d_get_far_switch_radius";
  
  FCSResult result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  mmm1d_get_far_switch_radius_2(handle->method_context, rad2);
  *rad2=sqrt(*rad2);
  return NULL;
}

FCSResult fcs_mmm1d_set_maxPWerror(FCS handle, fcs_float maxPWerror) {
  const char *fnc_name = "fcs_mmm1d_set_maxPWerror";
  
  FCSResult result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  mmm1d_set_maxPWerror(handle->method_context, maxPWerror);
  return NULL;
}

FCSResult fcs_mmm1d_get_maxPWerror(FCS handle, fcs_float *maxPWerror) {
  const char *fnc_name = "fcs_mmm1d_get_maxPWerror";
  
  FCSResult result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  mmm1d_get_maxPWerror(handle->method_context, maxPWerror);
  return NULL;
}

FCSResult fcs_mmm1d_set_bessel_cutoff(FCS handle, fcs_int cutoff) {
  const char *fnc_name = "fcs_mmm1d_set_bessel_cutoff";
  
  FCSResult result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  mmm1d_set_bessel_cutoff(handle->method_context, cutoff);
  return NULL;
}

FCSResult fcs_mmm1d_get_bessel_cutoff(FCS handle, fcs_int *cutoff) {
  const char *fnc_name = "fcs_mmm1d_get_bessel_cutoff";
  
  FCSResult result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  mmm1d_get_bessel_cutoff(handle->method_context, cutoff);
  return NULL;
}

FCSResult fcs_mmm1d_require_virial(FCS handle, fcs_int flag) { return NULL; }

FCSResult fcs_mmm1d_get_virial(FCS handle, fcs_float *virial) {
  const char* fnc_name =  "fcs_mmm1d_get_virial";
  FCSResult result = fcs_mmm1d_check(handle, fnc_name);
  if (result != NULL) return result;
  if (virial == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"null pointer supplied for virial"); 

  fcs_int i;
  for (i=0; i < 9; i++)
    virial[i] = 0.0;
  return NULL;
}

FCSResult fcs_mmm1d_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  const char *fnc_name = "fcs_mmm1d_set_parameter";

  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_far_switch_radius", mmm1d_set_far_switch_radius, FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_bessel_cutoff",     mmm1d_set_bessel_cutoff,     FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_maxPWerror",        mmm1d_set_maxPWerror,        FCS_PARSE_VAL(fcs_float));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_mmm1d_print_parameters(FCS handle)
{
  fcs_float radius;
  fcs_float PWerror;
  fcs_int cutoff;
  fcs_mmm1d_get_far_switch_radius(handle, &radius);
  fcs_mmm1d_get_bessel_cutoff(handle, &cutoff);
  fcs_mmm1d_get_maxPWerror(handle, &PWerror);
  printf("mmm1d bessel cutoff: %" FCS_LMOD_INT "d\n", cutoff);
  printf("mmm1d far switch radius: %e\n", radius);
  printf("mmm1d maximum PWerror: %e\n", PWerror);

  return FCS_RESULT_SUCCESS;
}
