/*
  Copyright (C) 2011-2012 Rene Halver, Olaf Lenz

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

/* method to check if p3m parameters are entered into FCS */
static FCSResult fcs_p3m_check(FCS handle, const char* fnc_name) {
  if (handle == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_METHOD_P3M)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"p3m\".");

  return NULL;
}

/* initialization function for basic p3m parameters */
FCSResult fcs_p3m_init(FCS handle)
{
  char* fnc_name = "fcs_p3m_init";
  FCSResult result;
  
  result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;

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
  handle->set_compute_virial = fcs_p3m_require_virial;
  handle->get_virial = fcs_p3m_get_virial;

  ifcs_p3m_init(&handle->method_context, handle->communicator);
  
  return NULL;
}

/* internal p3m-specific tuning function */
FCSResult fcs_p3m_tune(FCS handle, 
		       fcs_int local_particles,
		       fcs_float *positions, fcs_float *charges)
{
  char* fnc_name = "fcs_p3m_tune";
  FCSResult result;
  result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;

  /* Handle periodicity */
  const fcs_int *periodicity = fcs_get_periodicity(handle);
  if (! (periodicity[0] && periodicity[1] && periodicity[2]))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, 
      "p3m requires periodic boundary conditions.");
    
  /* Handle box size */
  const fcs_float *a = fcs_get_box_a(handle);
  const fcs_float *b = fcs_get_box_b(handle);
  const fcs_float *c = fcs_get_box_c(handle);
  if (!fcs_is_orthogonal(a, b, c)){
        if (ifcs_p3m_check_triclinic_box(a[1],a[2],b[2])){
            
            if(ifcs_p3m_set_triclinic_flag(handle->method_context)!=NULL)
           return ifcs_p3m_set_triclinic_flag(handle->method_context);           
        }
        else
            return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
                "p3m triclinic requires the box to be as follows: \n \
                the first box vector is parallel to the x axis\n \
                the second box vector is in the yz plane.");
    } else {
        if (!fcs_uses_principal_axes(a, b, c))
            return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
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
  result = ifcs_p3m_tune(handle->method_context, 
                         local_particles, max_local_particles,
			 positions, charges);
  return result;
}

/* internal p3m-specific run function */
FCSResult fcs_p3m_run(FCS handle, 
		      fcs_int local_particles,
		      fcs_float *positions, fcs_float *charges,
		      fcs_float *fields, fcs_float *potentials)
{
//   char* fnc_name = "fcs_p3m_run";
//   FCSResult result;

  fcs_p3m_tune(handle, local_particles, positions, charges);

  fcs_int max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  ifcs_p3m_run(handle->method_context,
		 local_particles, max_local_particles, positions, charges, fields, potentials);

  return NULL;
}

/* clean-up function for p3m */
FCSResult fcs_p3m_destroy(FCS handle)
{
  ifcs_p3m_destroy(handle->method_context);
  return NULL;
}

/******************************************************************************************************
 *
 *						Setter and Getter functions for p3m parameters
 *
 ******************************************************************************************************/


FCSResult fcs_p3m_set_r_cut(FCS handle, fcs_float r_cut) {
  const char *fnc_name = "fcs_p3m_set_r_cut";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_r_cut(handle->method_context, r_cut);
  return NULL;
}

FCSResult fcs_p3m_set_r_cut_tune(FCS handle) {
  const char *fnc_name = "fcs_p3m_set_r_cut_tune";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_r_cut_tune(handle->method_context);
  return NULL;
}

FCSResult fcs_p3m_get_r_cut(FCS handle, fcs_float *r_cut) {
  const char *fnc_name = "fcs_p3m_get_r_cut";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_get_r_cut(handle->method_context, r_cut);
  return NULL;
}


FCSResult fcs_p3m_set_alpha(FCS handle, fcs_float alpha) {
  const char *fnc_name = "fcs_p3m_set_alpha";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_alpha(handle->method_context, alpha);
  return NULL;
}

FCSResult fcs_p3m_set_alpha_tune(FCS handle) {
  const char *fnc_name = "fcs_p3m_set_alpha_tune";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_alpha_tune(handle->method_context);
  return NULL;
}

FCSResult fcs_p3m_get_alpha(FCS handle, fcs_float *alpha) {
  const char *fnc_name = "fcs_p3m_get_alpha";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_get_alpha(handle->method_context, alpha);
  return NULL;
}


FCSResult fcs_p3m_set_grid(FCS handle, fcs_int grid) {
  const char *fnc_name = "fcs_p3m_set_grid";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_grid(handle->method_context, grid);
  return NULL;
}

FCSResult fcs_p3m_set_grid_tune(FCS handle) {
  const char *fnc_name = "fcs_p3m_set_grid_tune";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_grid_tune(handle->method_context);
  return NULL;
}

FCSResult fcs_p3m_get_grid(FCS handle, fcs_int *grid) {
  const char *fnc_name = "fcs_p3m_get_grid";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_get_grid(handle->method_context, grid);
  return NULL;
}


FCSResult fcs_p3m_set_cao(FCS handle, fcs_int cao) {
  const char *fnc_name = "fcs_p3m_set_cao";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_cao(handle->method_context, cao);
  return NULL;
}

FCSResult fcs_p3m_set_cao_tune(FCS handle) {
  const char *fnc_name = "fcs_p3m_set_cao_tune";


  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_cao_tune(handle->method_context);
  return NULL;
}

FCSResult fcs_p3m_get_cao(FCS handle, fcs_int *cao) {
  const char *fnc_name = "fcs_p3m_get_cao";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_get_cao(handle->method_context, cao);
  return NULL;
}

FCSResult fcs_p3m_require_total_energy(FCS handle, fcs_int flag) {
  const char *fnc_name = "fcs_p3m_require_total_energy";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_require_total_energy(handle->method_context, flag);
  return NULL;
}

FCSResult fcs_p3m_get_total_energy(FCS handle, fcs_float *total_energy) {
  const char *fnc_name = "fcs_p3m_require_total_energy";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  return ifcs_p3m_get_total_energy(handle->method_context, total_energy);
}


FCSResult fcs_p3m_set_tolerance_field(FCS handle, fcs_float tolerance_field) {
  const char *fnc_name = "fcs_p3m_set_tolerance_field";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_tolerance_field(handle->method_context, tolerance_field);
  return NULL;
}

void fcs_p3m_set_tolerance_field_f(void *handle, fcs_float tolerance_field, fcs_int *return_value)
{
  FCSResult result = fcs_p3m_set_tolerance_field((FCS)handle, tolerance_field);
  printf("p3m got as accuracy: %10.4" FCS_LMOD_FLOAT "f \n", tolerance_field);
  if (NULL == result)
    *return_value = 0;
  else
    *return_value = fcs_result_get_return_code(result);
}

FCSResult fcs_p3m_set_tolerance_field_tune(FCS handle) {
  const char *fnc_name = "fcs_p3m_set_tolerance_field_tune";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_set_tolerance_field_tune(handle->method_context);
  return NULL;
}

FCSResult fcs_p3m_get_tolerance_field(FCS handle, fcs_float *tolerance_field) {
  const char *fnc_name = "fcs_p3m_get_tolerance_field";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  ifcs_p3m_get_tolerance_field(handle->method_context, tolerance_field);
  return NULL;
}

FCSResult fcs_p3m_get_near_parameters(FCS handle,
				      fcs_p3m_near_parameters_t *near_params) {
/*
  ifcs_p3m_get_alpha(handle->method_context, near_params);
*/
    fcs_float *alpha=(fcs_float*)malloc(sizeof(fcs_float));
    fcs_float *offset=(fcs_float*)malloc(sizeof(fcs_float));
    ifcs_p3m_get_near_params(handle->method_context, alpha, offset);
    near_params->alpha=*alpha;
    near_params->potentialOffset=*offset;
  return NULL;
}

FCSResult fcs_p3m_require_virial(FCS handle, fcs_int flag) { return NULL; }

FCSResult fcs_p3m_get_virial(FCS handle, fcs_float *virial) {
  const char* fnc_name =  "fcs_p3m_get_virial";
  FCSResult result = fcs_p3m_check(handle, fnc_name);
  if (result != NULL) return result;
  if (virial == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"null pointer supplied for virial"); 

  fcs_int i;
  for (i=0; i < 9; i++)
    virial[i] = 0.0;
  return NULL;
}

FCSResult fcs_p3m_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance)
{
  const char *fnc_name = "fcs_p3m_set_tolerance";

  if (tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
  {
    fcs_p3m_set_tolerance_field(handle, tolerance);
    return FCS_RESULT_SUCCESS;

  }
  
  return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "Unsupported tolerance type. P3M only supports FCS_TOLERANCE_TYPE_FIELD.");
}

FCSResult fcs_p3m_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance)
{
  *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
  fcs_p3m_get_tolerance_field(handle, tolerance);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  const char *fnc_name = "fcs_p3m_set_parameter";

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
  fcs_float tolerance;
  fcs_p3m_get_tolerance_field(handle, &tolerance);
  printf("p3m absolute field tolerance: %" FCS_LMOD_FLOAT "e\n", tolerance);

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_set_potential_shift(FCS handle, fcs_int flag){
    ifcs_p3m_set_potential_shift(handle->method_context, flag);
    return FCS_RESULT_SUCCESS;
}

FCSResult fcs_p3m_get_potential_shift(FCS handle, fcs_int *flag){
    ifcs_p3m_get_potential_shift(handle->method_context, flag);
    return FCS_RESULT_SUCCESS;
}