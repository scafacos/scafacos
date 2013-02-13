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

    maggs_init(&handle->method_context, communicator);

    return NULL;
}


extern FCSResult fcs_memd_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions,  fcs_float *charges)
{
    /* tune mesh and f_mass, calculate initial fields */
    return memd_tune_method(handle->method_context, local_particles, positions, charges);
}


extern FCSResult fcs_memd_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions,  fcs_float *charges, fcs_float *fields, fcs_float *potentials)
{
    /* retune if needed */
/*    if (memd_needs_retuning(handle->method_context, local_particles, positions, charges))
        fcs_memd_tune(handle, local_particles, local_max_particles, positions, charges);
*/    
    memd_run(handle->method_context,
            local_particles, local_max_particles, positions, charges, fields, potentials);
        
    return NULL;
}


extern FCSResult fcs_memd_destroy(FCS handle)
{
    return maggs_exit(handle->method_context);
}

FCSResult fcs_memd_require_virial(FCS handle, fcs_int flagvalue)
{
    return NULL;
}

FCSResult fcs_memd_get_virial(FCS handle, fcs_float* virial)
{
    return NULL;
}
