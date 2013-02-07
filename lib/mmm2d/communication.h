/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group

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
#ifndef _MMM2D_COMMUNICATION_H
#define _MMM2D_COMMUNICATION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <mpi.h>

typedef struct {
  /* The MPI communicator to use */
  MPI_Comm mpicomm;
  /* The original MPI communicator to use (possibly not cartesian) */
  MPI_Comm mpicomm_orig;
  /* The size of the communicator */
  fcs_int size;
  /* The rank within the communicator */
  fcs_int rank;

  /** The number of nodes in each spatial dimension. */
  fcs_int node_grid[3];
  /** position of this node in the node grid */
  fcs_int node_pos[3];
  /** the six nearest neighbors of a node in the node grid. */
  fcs_int node_neighbors[6];

  /** Size of the local box. */
  fcs_float local_box_l[3];
  /** Left (bottom, front) corner of this nodes local box. */ 
  fcs_float my_left[3];
  /** Right (top, back) corner of this nodes local box. */ 
  fcs_float my_right[3];
} mmm2d_comm_struct;

void mmm2d_comm_init(mmm2d_comm_struct *comm, MPI_Comm communicator);
void mmm2d_comm_destroy(mmm2d_comm_struct *comm);
#endif
