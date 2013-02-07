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

#include "init.h"
#include "types.h"
#include <stdlib.h>
#include <stdio.h>

void mmm2d_init(void **rd, MPI_Comm communicator) {
   mmm2d_data_struct *d;
  
  if (*rd == NULL) {
    /* allocate the memory for the mmm1d data structure */
    d = malloc(sizeof(mmm2d_data_struct));

    /* store the new pointer in rd */
    *rd = d;
  } else {
    d = (mmm2d_data_struct*)rd;
  }
  
  /* Init the communication stuff */
  mmm2d_comm_init(&d->comm, communicator);
  
  d->maxPWerror=MMM2D_DEFAULT_REQUIRED_ACCURACY;
  d->far_cut=MMM2D_DEFAULT_FAR_CUTOFF;
  d->far_cut2=d->far_cut*d->far_cut;
  d->calculate_far=MMM2D_DEFAULT_FAR_CUTOFF_CALCULATED;
  d->delta_mid_top=MMM2D_DEFAULT_DCONTRAST_TOP;
  d->delta_mid_bot=MMM2D_DEFAULT_DCONTRAST_BOT;
  d->layers_per_node=MMM2D_DEFAULT_N_LAYERS;
  d->n_total_layers=MMM2D_DEFAULT_N_LAYERS*d->comm.size+2;
  d->skin=MMM2D_DEFAULT_SKIN;
  d->my_bottom=0;
  
  d->require_total_energy=1;
  
  d->box_l[0]=0.; d->box_l[1]=0.; d->box_l[2]=0.;
  
  d->besselCutoff.e = NULL;
  d->besselCutoff.n = 0;
  d->besselCutoff.max = 0;
  
  d->bon.e = NULL;
  d->bon.n = 0;
  d->bon.max = 0;
  
  d->polTaylor=NULL;
  d->polTaylor = malloc(sizeof(mmm_data_struct));
  (d->polTaylor)->modPsi=NULL;
  (d->polTaylor)->n_modPsi=0;
  
  d->local_charges= NULL;
  d->local_positions= NULL;
  d->zslices_nparticles= NULL;
  d->n_localpart=0;
  
  d->partblk = NULL;
  d->lclcblk = NULL;
  d->gblcblk = NULL;

  d->scxcache = NULL;
  d->n_scxcache = 0;
  d->scycache = NULL;
  d->n_scycache = 0;
  
  d->needs_tuning=1;
  
  if(d->delta_mid_top!=0. || d->delta_mid_bot!=0.){
    d->dielectric_contrast_on=1;
    d->delta_mult = d->delta_mid_top*d->delta_mid_bot;
  } else {
    d->dielectric_contrast_on=0;
    d->delta_mult = 0;
  }
}

/* safe free */
static void sfree(void* ptr) {
  if (ptr != NULL)
    free(ptr);
}

void mmm2d_destroy(void *rd) {
  if (rd != NULL) {
    mmm2d_data_struct *d = (mmm2d_data_struct*)rd;
    sfree(d);
  }
}
