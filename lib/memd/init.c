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
#include <string.h>

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
FCSResult ifcs_memd_init(void** rawdata, MPI_Comm communicator)
{
    memd_struct* memd;
    if (*rawdata == NULL) {
        memd = calloc(1,sizeof(memd_struct));
        memset(memd, 0, sizeof(memd_struct));
        *rawdata = memd;
        printf("New handle created!\n"); fflush(stdout);
    } else {
        memd = (memd_struct*) rawdata;
    }


    fcs_memd_setup_communicator(memd, communicator);
    
    //if(!this_node) fprintf(stderr, "%d: Electric field is initialized\n", this_node);
    return ifcs_memd_sanity_checks(memd);
}

/** Frees the dynamically allocated memory
 Currently not called from anywhere. */
FCSResult fcs_memd_exit(void* rawdata)
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
