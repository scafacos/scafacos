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
#include <stdio.h>

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
  mmm1d_comm_init(&d->comm, communicator);
  
  /* Init the default MMM1D parameters */
  d->far_switch_radius_2=-1;
  d->bessel_cutoff=MMM1D_DEFAULT_MAXIMAL_B_CUT;
  d->maxPWerror=MMM1D_DEFAULT_REQUIRED_ACCURACY;
  d->total_charge=0.;
  
  d->box_l[0]=0.; d->box_l[1]=0.; d->box_l[2]=0.;
  
  d->local_charges= NULL;
  d->local_positions= NULL;
  d->n_localpart=0;
  
  d->polTaylor=NULL;
  d->polTaylor = malloc(sizeof(mmm_data_struct));
  (d->polTaylor)->modPsi=NULL;
  (d->polTaylor)->n_modPsi=0;
  d->bessel_radii = NULL;
}

/* safe free */
static void sfree(void* ptr) {
  if (ptr != NULL)
    free(ptr);
}

void mmm1d_destroy(void *rd) {
  if (rd != NULL) {
    mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
    sfree(d);
  }
}
