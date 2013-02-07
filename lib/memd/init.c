/*
 Copyright (C) 2010/2011/2012 Florian Fahrenberger
 
 This file is part of ScaFaCoS.
 
 ScaFaCoS is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 ScaFaCoS is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "init.h"
#include "data_types.h"
#include "communication.h"
#include "helper_functions.h"
#include "initial_solution.h"

/***************************/
/****** init and exit ******/
/***************************/

/** Initialization function.
 Sets maggs structure variables.
 Calls calculation of initial D-field. */
FCSResult maggs_init(void** rawdata, MPI_Comm communicator)
{
    memd_struct* memd;
    if (*rawdata == NULL) {
        memd = calloc(1,sizeof(memd_struct));
        *rawdata = memd;
        printf("New handle created!\n"); fflush(stdout);
    } else {
        memd = (memd_struct*) rawdata;
    }


    FCSResult result;
    if ( (memd->parameters.mesh>0) && (memd->parameters.box_length[0]>0.0) ) {
        memd->parameters.inva  = (fcs_float) memd->parameters.mesh/memd->parameters.box_length[0];
        memd->parameters.a     = 1.0/memd->parameters.inva;
        maggs_setup_local_lattice(memd);
    } else {
        int k;
//        fprintf(stdout, "box_l: %f\n", memd->parameters.box_length[0]); fflush(stdout);
        memd->parameters.mesh=32;
        if (memd->parameters.box_length[0]<ROUND_ERROR_PREC) {
            FOR3D(k) memd->parameters.box_length[k] = 10.0;
        }
        memd->parameters.inva  = (fcs_float) memd->parameters.mesh/memd->parameters.box_length[0];
        memd->parameters.a     = 1.0/memd->parameters.inva;
        maggs_setup_local_lattice(memd);        
    }
    
    maggs_setup_communicator(memd, communicator);
    
    result = maggs_sanity_checks(memd);

    //if(!this_node) fprintf(stderr, "%d: Electric field is initialized\n", this_node);
    return result;
}

/** Frees the dynamically allocated memory
 Currently not called from anywhere. */
FCSResult maggs_exit(void* rawdata)
{
    memd_struct* memd = (memd_struct*) rawdata;
    free(memd->lattice);
    free(memd->neighbor);
    free(memd->local_cells.cell);
    free(memd->ghost_cells.cell);
    free(memd->Dfield);
    free(memd->Bfield);
    return NULL;
}
