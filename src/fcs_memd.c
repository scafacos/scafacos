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

#include <math.h>
#include "mpi.h"
#include <stdio.h>

#include "FCSDefinitions.h"
#include "FCSResult_p.h"

#include "fcs_memd.h"
#include "memd/memd.h"
#include "memd/data_types.h"

static FCSResult fcs_memd_check(FCS handle, const char* fnc_name) {
    if (handle == NULL)
        return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");
    
    if (fcs_get_method(handle) != FCS_METHOD_MEMD)
        return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"memd\".");
    
    return NULL;
}

FCSResult fcs_memd_init(FCS handle)
{
    char* fnc_name = "fcs_memd_init";
    FCSResult result;
    result = fcs_memd_check(handle, fnc_name);
    if (result != NULL) return result;

    handle->tune = fcs_memd_tune;
    handle->run = fcs_memd_run;

    ifcs_memd_init(&handle->method_context, handle->communicator);
    
    return NULL;
}


FCSResult fcs_memd_tune(FCS handle, fcs_int local_particles, fcs_float *positions,  fcs_float *charges)
{
    char* fnc_name = "fcs_memd_tune";
    FCSResult result;
    
    /* Handle periodicity */
    const fcs_int *periodicity = fcs_get_periodicity(handle);
    if (! (periodicity[0] && periodicity[1] && periodicity[2]))
        return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
                                "memd requires periodic boundary conditions.");
    
    /* Handle box size */
    const fcs_float *a = fcs_get_box_a(handle);
    const fcs_float *b = fcs_get_box_b(handle);
    const fcs_float *c = fcs_get_box_c(handle);
    if (!fcs_is_orthogonal(a, b, c))
        return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
                                "memd requires the box to be orthorhombic.");
    
    if (!fcs_uses_principal_axes(a, b, c))
        return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
                                "memd requires the box vectors to be parallel to the principal axes.");
    
    fcs_memd_set_box_size(handle->method_context, a[0], b[1], c[2]);
    
    /* tune mesh and f_mass, calculate initial fields */
    return fcs_memd_tune_method(handle->method_context, local_particles, positions, charges);
}


FCSResult fcs_memd_run(FCS handle, fcs_int local_particles, fcs_float *positions,  fcs_float *charges, fcs_float *fields, fcs_float *potentials)
{
    /* retune if needed */
/*    if (fcs_memd_needs_retuning(handle->method_context, local_particles, positions, charges))
        fcs_memd_tune(handle, local_particles, local_max_particles, positions, charges);
*/
    fcs_int max_local_particles = fcs_get_max_local_particles(handle);
    if (local_particles > max_local_particles) max_local_particles = local_particles;

    ifcs_memd_run(handle->method_context,
            local_particles, max_local_particles, positions, charges, fields, potentials);
        
    return NULL;
}


FCSResult fcs_memd_destroy(FCS handle)
{
  FCSResult result = fcs_memd_exit(handle->method_context);

  free(handle->memd_param);
  handle->memd_param = NULL;

  return result;
}

FCSResult fcs_memd_require_virial(FCS handle, fcs_int flagvalue)
{
    return NULL;
}

FCSResult fcs_memd_get_virial(FCS handle, fcs_float* virial)
{
    return NULL;
}

FCSResult fcs_memd_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
/*  const char *fnc_name = "fcs_memd_set_parameter";*/

  char *param = *current;
  char *cur = *next;

  *matched = 0;

/*  FCS_PARSE_IF_PARAM_THEN_FUNC3_GOTO_NEXT("", fcs_memd_set_box_size,                  FCS_PARSE_VAL(fcs_float), FCS_PARSE_VAL(fcs_float), FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_time_step,                 FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_total_number_of_particles, FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_local_number_of_particles, FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_init_flag,                 FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_mesh_size_1D,              FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_speed_of_light,            FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_permittivity,              FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_temperature,               FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_bjerrum_length,            FCS_PARSE_VAL(fcs_float));*/

  return FCS_RESULT_SUCCESS;

/*next_param:*/
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_memd_print_parameters(FCS handle)
{
  printf("not implemented\n");
  return FCS_RESULT_SUCCESS;
}
