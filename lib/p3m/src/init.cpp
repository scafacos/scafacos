/*
  Copyright (C) 2013,2014 Olaf Lenz
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
#include "utils.hpp"
#include "p3m.hpp"

namespace P3M {

  void init(data_struct *d, MPI_Comm communicator) {
    P3M_DEBUG(printf( "P3M::init() started...\n"));

    /* Init the communication stuff */
    comm_init(&d->comm, communicator);

    /* Init the P3M parameters */
    d->box_l[0] = 1.0;
    d->box_l[1] = 1.0;
    d->box_l[2] = 1.0;
    d->skin = 0.0;
    d->tolerance_field = P3M_DEFAULT_TOLERANCE_FIELD;
    d->n_interpol = P3M_DEFAULT_N_INTERPOL;

    /* Tunable parameters */
    d->r_cut = 0.0;
    d->alpha = 0.0;
    d->grid[0] = 0;
    d->grid[1] = 0;
    d->grid[2] = 0;
    d->cao = 0;

    /* Everything needs to be retuned at the beginning */
    d->needs_retune = 1;
    d->tune_r_cut = 1;
    d->tune_alpha = 1;
    d->tune_grid = 1;
    d->tune_cao = 1;

    /* Which components to compute? */
    d->require_total_energy = 0;
    d->total_energy = 0.0;

#ifdef P3M_PRINT_TIMINGS
    d->require_timings = FULL;
#else
    d->require_timings = NONE;
#endif
    for (int i=0; i < NUM_TIMINGS; i++)
      d->timings[i] = 0.0;

    /* Init the derived params */
    d->grid_off[0] = P3M_DEFAULT_GRIDOFF;
    d->grid_off[1] = P3M_DEFAULT_GRIDOFF;
    d->grid_off[2] = P3M_DEFAULT_GRIDOFF;
    d->cao_cut[0] = 0.0;
    d->cao_cut[1] = 0.0;
    d->cao_cut[2] = 0.0;
    d->a[0] = 0.0;
    d->a[1] = 0.0;
    d->a[2] = 0.0;
    d->ai[0] = 0.0;
    d->ai[1] = 0.0;
    d->ai[2] = 0.0;
    d->additional_grid[0] = 0.0;
    d->additional_grid[1] = 0.0;
    d->additional_grid[2] = 0.0;

    /* init the P3M data */
    d->sum_qpart = 0;
    d->sum_q2 = 0.0;
    d->square_sum_q = 0.0;

    d->caf = NULL;
    d->cafx = NULL;
    d->cafy = NULL;
    d->cafz = NULL;
    d->caf_d = NULL;
    d->cafx_d = NULL;
    d->cafy_d = NULL;
    d->cafz_d = NULL;

    d->pos_shift = 0.0;
    d->meshift_x = NULL;
    d->meshift_y = NULL;
    d->meshift_z = NULL;

    d->d_op[0] = NULL;
    d->d_op[1] = NULL;
    d->d_op[2] = NULL;
    d->g_force = NULL;
    d->g_energy = NULL;

    d->ks_pnum = 0;

    d->send_grid = NULL;
    d->recv_grid = NULL;

    /* init the fft */
    d->fft = new Parallel3DFFT(&d->comm);

    P3M_DEBUG(printf( "P3M::init() finished.\n"));
  }

  void destroy(data_struct *d) {
    d->fft->free_data(d->rs_grid);
    d->fft->free_data(d->ks_grid);
    delete d->fft;
    comm_destroy(&d->comm);
    sfree(d->send_grid);
    sfree(d->recv_grid);
    sdelete(d->caf);
    sdelete(d->cafx);
    sdelete(d->cafy);
    sdelete(d->cafz);
    sdelete(d->caf_d);
    sdelete(d->cafx_d);
    sdelete(d->cafy_d);
    sdelete(d->cafz_d);
    sfree(d->g_energy);
    sfree(d->g_force);
    for (p3m_int i=0; i<3; i++) {
      sfree(d->d_op[i]);
    }
    
    sfree(d);
  }

}
