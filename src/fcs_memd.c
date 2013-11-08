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
        return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");
    
    if (fcs_get_method(handle) != FCS_MEMD)
        return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"memd\".");
    
    return NULL;
}

extern FCSResult fcs_memd_init(FCS handle, MPI_Comm communicator)
{
    char* fnc_name = "fcs_memd_init";
    FCSResult result;
    result = fcs_memd_check(handle, fnc_name);
    if (result != NULL) return result;

    ifcs_memd_init(&handle->method_context, communicator);
    
    return NULL;
}


extern FCSResult fcs_memd_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions,  fcs_float *charges)
{
    char* fnc_name = "fcs_memd_tune";
    FCSResult result;
    
    /* Handle periodicity */
    fcs_int *periodicity = fcs_get_periodicity(handle);
    if (! (periodicity[0] && periodicity[1] && periodicity[2]))
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
                                "memd requires periodic boundary conditions.");
    
    /* Handle box size */
    fcs_float *a = fcs_get_box_a(handle);
    fcs_float *b = fcs_get_box_b(handle);
    fcs_float *c = fcs_get_box_c(handle);
    if (!fcs_is_orthogonal(a, b, c))
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
                                "memd requires the box to be orthorhombic.");
    
    if (!fcs_uses_principal_axes(a, b, c))
        return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
                                "memd requires the box vectors to be parallel to the principal axes.");
    
    fcs_memd_set_box_size(handle->method_context, a[0], b[1], c[2]);
    
    /* tune mesh and f_mass, calculate initial fields */
    return fcs_memd_tune_method(handle->method_context, local_particles, positions, charges);
}


extern FCSResult fcs_memd_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions,  fcs_float *charges, fcs_float *fields, fcs_float *potentials)
{
    /* retune if needed */
/*    if (fcs_memd_needs_retuning(handle->method_context, local_particles, positions, charges))
        fcs_memd_tune(handle, local_particles, local_max_particles, positions, charges);
*/
    ifcs_memd_run(handle->method_context,
            local_particles, local_max_particles, positions, charges, fields, potentials);
        
    return NULL;
}


extern FCSResult fcs_memd_destroy(FCS handle)
{
    return fcs_memd_exit(handle->method_context);
}

FCSResult fcs_memd_require_virial(FCS handle, fcs_int flagvalue)
{
    return NULL;
}

FCSResult fcs_memd_get_virial(FCS handle, fcs_float* virial)
{
    return NULL;
}
