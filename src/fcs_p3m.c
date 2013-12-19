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
#include "p3m/p3m.hpp"
#include "p3m/types.hpp"

/* method to check if p3m parameters are entered into FCS */
static FCSResult fcs_p3m_check(FCS handle, const char* fnc_name) {
    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

    if (fcs_get_method(handle) != FCS_P3M)
        return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"p3m\".");

    return NULL;
}

/* initialization function for basic p3m parameters */
FCSResult fcs_p3m_init(FCS handle, MPI_Comm communicator) {
    char* fnc_name = "fcs_p3m_init";
    FCSResult result;

    result = fcs_p3m_check(handle, fnc_name);
    if (result != NULL) return result;

    ifcs_p3m_init(&handle->method_context, communicator);

    return NULL;
}

/* internal p3m-specific tuning function */
FCSResult fcs_p3m_tune(FCS handle,
        fcs_int local_particles, fcs_int local_max_particles,
        fcs_float *positions, fcs_float *charges) {
    char* fnc_name = "fcs_p3m_tune";
    FCSResult result;
    result = fcs_p3m_check(handle, fnc_name);
    if (result != NULL) return result;

    /* Handle periodicity */
    fcs_int *periodicity = fcs_get_periodicity(handle);
    if (!(periodicity[0] && periodicity[1] && periodicity[2]))
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
            "p3m requires periodic boundary conditions.");

    /* Handle box size */
    fcs_float *a = fcs_get_box_a(handle);
    fcs_float *b = fcs_get_box_b(handle);
    fcs_float *c = fcs_get_box_c(handle);
    if (!fcs_is_orthogonal(a, b, c)){
/*
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
            "p3m requires the box to be orthorhombic.");
*/
/*
        #define P3M_TRICLINC
*/
        ifcs_p3m_set_triclinic_flag(handle->method_context);
}
/*
    if (!fcs_uses_principal_axes(a, b, c))
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
            "p3m requires the box vectors to be parallel to the principal axes.");
*/
  
    ifcs_p3m_set_box_a(handle->method_context, a[0]);
    ifcs_p3m_set_box_b(handle->method_context, b[1]);
    ifcs_p3m_set_box_c(handle->method_context, c[2]);

    ifcs_p3m_set_box_geometry(handle->method_context, a, b, c);
    
    
    ifcs_p3m_set_near_field_flag(handle->method_context,
            fcs_get_near_field_flag(handle));

    /* Effectively, tune initializes the algorithm. */
    result = ifcs_p3m_tune(handle->method_context,
            local_particles, local_max_particles,
            positions, charges);
    return result;
}

/* internal p3m-specific run function */
FCSResult fcs_p3m_run(FCS handle,
        fcs_int local_particles, fcs_int local_max_particles,
        fcs_float *positions, fcs_float *charges,
        fcs_float *fields, fcs_float *potentials) {
    //   char* fnc_name = "fcs_p3m_run";
    //   FCSResult result;

    fcs_p3m_tune(handle, local_particles, local_max_particles, positions, charges);

    ifcs_p3m_run(handle->method_context,
            local_particles, local_max_particles,
            positions, charges,
            fields, potentials);

    return NULL;
}

/* clean-up function for p3m */
FCSResult fcs_p3m_destroy(FCS handle) {
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

void fcs_p3m_set_tolerance_field_f(void *handle, fcs_float tolerance_field, fcs_int *return_value) {
    FCSResult result = fcs_p3m_set_tolerance_field((FCS) handle, tolerance_field);
    printf("p3m got as accuracy: %10.4" FCS_LMOD_FLOAT "f \n", tolerance_field);
    if (NULL == result)
        *return_value = 0;
    else
        *return_value = fcsResult_getReturnCode(result);
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
    ifcs_p3m_get_alpha(handle->method_context, near_params);
    return NULL;
}

FCSResult fcs_p3m_require_virial(FCS handle, fcs_int flag) {
    return NULL;
}

FCSResult fcs_p3m_get_virial(FCS handle, fcs_float *virial) {
    const char* fnc_name = "fcs_p3m_get_virial";
    FCSResult result = fcs_p3m_check(handle, fnc_name);
    if (result != NULL) return result;
    if (virial == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied for virial");

    fcs_int i;
    for (i = 0; i < 9; i++)
        virial[i] = 0.0;
    return NULL;
}

