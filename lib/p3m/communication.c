/*
  Copyright (C) 2011,2012 Olaf Lenz
  
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
#include "utils.h"

/****************************************************/
/* IMPLEMENTATION */
/****************************************************/
void
ifcs_p3m_comm_init(ifcs_p3m_comm_struct *comm, MPI_Comm communicator) {
  P3M_DEBUG(printf( "  ifcs_p3m_comm_init() started...\n"));

  /* store the original communicator */
  comm->mpicomm_orig = communicator;
  MPI_Comm_size(communicator, &comm->size);

  P3M_DEBUG(printf( "  ifcs_p3m_comm_init() finished.\n"));
}

void 
ifcs_p3m_comm_prepare(ifcs_p3m_comm_struct *comm, fcs_float *box_l) {
  P3M_DEBUG(printf( "  ifcs_p3m_comm_prepare() started...\n"));

  /* Test whether the communicator is cartesian and correct dimensionality */
  int comm_is_cart = 0;
  int status;

  MPI_Topo_test(comm->mpicomm_orig, &status);
  if (status == MPI_CART) {
    /* Communicator is cartesian, so test dimensionality */
    int ndims;
    MPI_Cartdim_get(comm->mpicomm_orig, &ndims);
    if (ndims == 3) {
      /* Correct dimensionality, so get grid and test periodicity */
      int periodicity[3];
      MPI_Cart_get(comm->mpicomm_orig, 3, comm->node_grid, 
		   periodicity, comm->node_pos);
      if (periodicity[0] && periodicity[1] && periodicity[2]) {
	/* If periodicity is correct, we can just use this communicator */
	comm->mpicomm = comm->mpicomm_orig;
	/* get the rank */
	MPI_Comm_rank(comm->mpicomm, &comm->rank);
	comm_is_cart = 1;
      }
    }
  }

  /* otherwise, we have to set up the cartesian communicator */
  if (!comm_is_cart) {
    P3M_DEBUG(printf( "    Setting up cartesian communicator...\n"));

    comm->node_grid[0] = 0;
    comm->node_grid[1] = 0;
    comm->node_grid[2] = 0;

    /* compute node grid */
    MPI_Dims_create(comm->size, 3, comm->node_grid);

    P3M_INFO(printf("    node_grid=%dx%dx%d\n",				\
		    comm->node_grid[0],					\
		    comm->node_grid[1],					\
		    comm->node_grid[2]));

    /* create communicator */
    int periodicity[3] = {1, 1, 1};
    MPI_Cart_create(comm->mpicomm_orig, 3, comm->node_grid, 
		    periodicity, 1, &comm->mpicomm);

    /* get the rank */
    MPI_Comm_rank(comm->mpicomm, &comm->rank);
    /* get node pos */
    MPI_Cart_coords(comm->mpicomm, comm->rank, 3, comm->node_pos);
  }
    
  /* fetch neighborhood info */
  for (int dir = 0; dir < 3; dir++) {
    MPI_Cart_shift(comm->mpicomm, dir, 1, 
		   &comm->node_neighbors[2*dir],
		   &comm->node_neighbors[2*dir+1]); 
    P3M_DEBUG_LOCAL(printf( "    %d: dir=%d: n1=%d n2=%d\n", comm->rank, dir, \
			    comm->node_neighbors[2*dir],		\
			    comm->node_neighbors[2*dir+1]));
  }

  /* init local points */
  for (int i=0; i< 3; i++) {
    comm->local_box_l[i] = 0.0;
    comm->my_left[i] = 0.0;
    comm->my_right[i] = 0.0;
  }

  /* compute box limits */
  for(fcs_int i = 0; i < 3; i++) {
    comm->local_box_l[i] = box_l[i]/(fcs_float)comm->node_grid[i]; 
    comm->my_left[i]   = comm->node_pos[i]    *comm->local_box_l[i];
    comm->my_right[i]  = (comm->node_pos[i]+1)*comm->local_box_l[i];    
  }
  P3M_DEBUG(printf("    local_box_l=(%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f)\n"                  \
                   "    my_left=(%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f)\n"                      \
                   "    my_right=(%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f)\n",                    \
                   comm->local_box_l[0],                                \
                   comm->local_box_l[1],                                \
                   comm->local_box_l[2],                                \
                   comm->my_left[0], comm->my_left[1], comm->my_left[2], \
                   comm->my_right[0], comm->my_right[1], comm->my_right[2] \
                   ));
  


  P3M_DEBUG(printf("  ifcs_p3m_comm_prepare() finished.\n"));
}

void ifcs_p3m_comm_destroy(ifcs_p3m_comm_struct *comm) {
  /* if (comm->mpicomm != comm->mpicomm_orig) */
  /*   MPI_Comm_free(comm->mpicomm); */
}
