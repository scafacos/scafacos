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

#ifndef _MEMD_COMMUNICATION_H
#define _MEMD_COMMUNICATION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>
#include "data_types.h"

/** check and set up communicator */
void fcs_memd_setup_communicator(memd_struct* memd, MPI_Comm communicator);
/** set up lattice structure and all parameters */
void fcs_memd_setup_local_lattice(memd_struct* memd);
/** communicate surface patches */
void fcs_memd_exchange_surface_patch(memd_struct* memd, fcs_float *field, fcs_int dim, fcs_int e_equil);

#endif