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
#include "Communication.hpp"
#include "utils.hpp"

namespace P3M {
  /****************************************************/
  /* IMPLEMENTATION */
  /****************************************************/
  Communication::Communication(MPI_Comm mpicomm) {
	  /* store the original communicator */
	  this->mpicomm = mpicomm;
	  mpicomm_orig = mpicomm;
	  MPI_Comm_size(mpicomm, &size);
	  MPI_Comm_rank(mpicomm, &rank);
  }

  Communication::~Communication() {
	  if (mpicomm != mpicomm_orig)
		  MPI_Comm_free(&mpicomm);
  }

  void Communication::prepare(p3m_float box_l[3]) {
    P3M_DEBUG(printf( "  P3M::Communication::prepare() started...\n"));

    /* Test whether the communicator is cartesian and correct dimensionality */
    bool comm_is_cart = false;
    int status;

    MPI_Topo_test(mpicomm_orig, &status);
    if (status == MPI_CART) {
      /* Communicator is cartesian, so test dimensionality */
      int ndims;
      MPI_Cartdim_get(mpicomm_orig, &ndims);
      if (ndims == 3) {
        /* Correct dimensionality, so get grid and test periodicity */
        int periodicity[3];
        MPI_Cart_get(mpicomm_orig, 3, node_grid,
                     periodicity, node_pos);
        if (periodicity[0] && periodicity[1] && periodicity[2]) {
          /* If periodicity is correct, we can just use this communicator */
          mpicomm = mpicomm_orig;
          /* get the rank */
          MPI_Comm_rank(mpicomm, &rank);
          comm_is_cart = true;
        }
      }
    }

    /* otherwise, we have to set up the cartesian communicator */
    if (!comm_is_cart) {
      P3M_DEBUG(printf( "    Setting up cartesian communicator...\n"));

      node_grid[0] = 0;
      node_grid[1] = 0;
      node_grid[2] = 0;

      /* compute node grid */
      MPI_Dims_create(size, 3, node_grid);

#ifdef P3M_ENABLE_INFO
      if (onMaster())
    	  printf("    node_grid=%dx%dx%d\n", node_grid[0], node_grid[1], node_grid[2]);
#endif

      /* create communicator */
      int periodicity[3] = {1, 1, 1};
      MPI_Cart_create(mpicomm_orig, 3, node_grid,
                      periodicity, 1, &mpicomm);

      /* get the rank */
      MPI_Comm_rank(mpicomm, &rank);
      /* get node pos */
      MPI_Cart_coords(mpicomm, rank, 3, node_pos);
    }
    
    /* fetch neighborhood info */
    for (int dir = 0; dir < 3; dir++) {
      MPI_Cart_shift(mpicomm, dir, 1,
                     &node_neighbors[2*dir],
                     &node_neighbors[2*dir+1]);
      P3M_DEBUG_LOCAL(printf( "    %d: dir=%d: n1=%d n2=%d\n", rank, dir, \
                              node_neighbors[2*dir],		\
                              node_neighbors[2*dir+1]));
    }

    /* init local points */
    for (int i=0; i< 3; i++) {
      local_box_l[i] = 0.0;
      my_left[i] = 0.0;
      my_right[i] = 0.0;
    }

    /* compute box limits */
    for(p3m_int i = 0; i < 3; i++) {
      local_box_l[i] = box_l[i]/(p3m_float)node_grid[i];
      my_left[i]   = node_pos[i]    *local_box_l[i];
      my_right[i]  = (node_pos[i]+1)*local_box_l[i];
    }
    P3M_DEBUG(printf("    local_box_l=" F3FLOAT "\n"                      \
                     "    my_left=" F3FLOAT "\n"                          \
                     "    my_right=" F3FLOAT "\n",                        \
                     local_box_l[0],                                \
                     local_box_l[1],                                \
                     local_box_l[2],                                \
                     my_left[0], my_left[1], my_left[2], \
                     my_right[0], my_right[1], my_right[2] \
                     ));

    P3M_DEBUG(printf("  P3M::Communication::prepare() finished.\n"));
  }
}
