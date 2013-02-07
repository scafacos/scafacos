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
#include "communication.h"

/****************************************************/
/* IMPLEMENTATION */
/****************************************************/
void
mmm1d_comm_init(mmm1d_comm_struct *comm, MPI_Comm communicator) {
  /* store the original communicator */
  comm->mpicomm_orig = communicator;
  MPI_Comm_size(communicator, &comm->size);

  /* Test whether the communicator is cartesian and correct dimensionality */
  fcs_int comm_is_cart = 0;
  fcs_int status;
  MPI_Topo_test(communicator, &status);
  if (status == MPI_CART) {
    /* Communicator is cartesian, so test dimensionality */
    fcs_int ndims;
    MPI_Cartdim_get(communicator, &ndims);
    if (ndims == 3) {
      /* Correct dimensionality, so get grid and test periodicity */
      fcs_int periodicity[3];
      MPI_Cart_get(communicator, 3, comm->node_grid, periodicity, comm->node_pos);
      if (periodicity[0] && periodicity[1] && periodicity[2]) {
        /* If periodicity is correct, we can just use this communicator */
        comm->mpicomm = communicator;
        /* get the rank */
        MPI_Comm_rank(communicator, &comm->rank);
        comm_is_cart = 1;
      }
    }
  }

  /* otherwise, we have to set up the cartesian communicator */
  if (!comm_is_cart) {

    comm->node_grid[0] = 0.0;
    comm->node_grid[1] = 0.0;
    comm->node_grid[2] = 0.0;

    /* compute node grid */
    MPI_Dims_create(comm->size, 3, comm->node_grid);

    /* create communicator */
    fcs_int periodicity[3] = {0, 0, 1};
    MPI_Cart_create(comm->mpicomm_orig, 3, comm->node_grid, periodicity, 1, &comm->mpicomm);

    /* get the rank */
    MPI_Comm_rank(communicator, &comm->rank);
    /* get node pos */
    MPI_Cart_coords(comm->mpicomm, comm->rank, 3, comm->node_pos);
  }
}

void mmm1d_comm_prepare(mmm1d_comm_struct *comm, fcs_float *box_l) {
  /* compute box limits */
  /*
  for(fcs_int i = 0; i < 3; i++) {
    comm->local_box_l[i] = box_l[i]/(fcs_float)comm->node_grid[i];
    comm->my_left[i]   = comm->node_pos[i]    *comm->local_box_l[i];
    comm->my_right[i]  = (comm->node_pos[i]+1)*comm->local_box_l[i];
  }
  */
}

void mmm1d_comm_destroy(mmm1d_comm_struct *comm) {
  /**/
}
