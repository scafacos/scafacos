/*
  Copyright (C) 2011
  
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

#include "init.h"
#include "types.h"
#include <stdlib.h>
// #include <stdio.h>

static MPI_Comm cart_comm_1d(MPI_Comm communicator) {
  /* test whether the communicator is cartesian and correct dimensionality */
  fcs_int status;
  MPI_Topo_test(communicator, &status);
  if (status == MPI_CART) {
    /* Communicator is cartesian, so test dimensionality */
    fcs_int ndims;
    MPI_Cartdim_get(communicator, &ndims);
    if (ndims == 3) {
      /* Correct dimensionality, so get grid and test periodicity */
      fcs_int node_grid[3], periodicity[3], node_pos[3];
      MPI_Cart_get(communicator, 3, node_grid, periodicity, node_pos);
      if (!periodicity[0] && !periodicity[1] && periodicity[2]) {
        return communicator;
      }
    }
  }
  /* Otherwise, we have to set up the cartesian communicator */
  fcs_int comm_size;
  MPI_Comm_size(communicator, &comm_size);
  fcs_int node_grid[3] = {0, 0, 0};
  MPI_Dims_create(comm_size, 3, node_grid);
  fcs_int periodicity[3] = {1, 1, 0};

  MPI_Comm new_communicator;
  MPI_Cart_create(communicator, 3, node_grid, periodicity, 1, &new_communicator);

  return new_communicator;
}

void mmm1d_init(void **rd, MPI_Comm communicator) {
  
  mmm1d_data_struct *d;
  
  if (*rd == NULL) {
    /* allocate the memory for the mmm1d data structure */
    d = malloc(sizeof(mmm1d_data_struct));

    /* store the new pointer in rd */
    *rd = d;
  } else {
    d = (mmm1d_data_struct*)rd;
  }
  
  /* Init the communication stuff */
  d->comm = cart_comm_1d(communicator);

  /* Init the default MMM1D parameters */
  d->far_switch_radius_2=-1;
  d->bessel_cutoff=MMM1D_DEFAULT_MAXIMAL_B_CUT;
  d->maxPWerror=MMM1D_DEFAULT_REQUIRED_ACCURACY;
  
  d->box_l[0]=0.; d->box_l[1]=0.; d->box_l[2]=0.;
  
  d->polTaylor=NULL;
  d->polTaylor = malloc(sizeof(mmm_data_struct));
  (d->polTaylor)->modPsi=NULL;
  (d->polTaylor)->n_modPsi=0;
  d->bessel_radii = NULL;
}

void mmm1d_destroy(void *rd) {
  if (rd != NULL) {
    mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
    free(d);
  }
}
