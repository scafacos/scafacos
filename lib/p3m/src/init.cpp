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

	data_struct::data_struct(MPI_Comm mpicomm) :
		comm(mpicomm), fft(comm) {
		P3M_DEBUG(printf( "P3M::P3M() started...\n"));

		/* Init the P3M parameters */
		box_l[0] = 1.0;
		box_l[1] = 1.0;
		box_l[2] = 1.0;
		skin = 0.0;
		tolerance_field = P3M_DEFAULT_TOLERANCE_FIELD;
		n_interpol = P3M_DEFAULT_N_INTERPOL;

		/* Tunable parameters */
		r_cut = 0.0;
		alpha = 0.0;
		grid[0] = 0;
		grid[1] = 0;
		grid[2] = 0;
		cao = 0;

		/* Everything needs to be retuned at the beginning */
		needs_retune = 1;
		tune_r_cut = 1;
		tune_alpha = 1;
		tune_grid = 1;
		tune_cao = 1;

		/* Which components to compute? */
		require_total_energy = 0;
		total_energy = 0.0;

#ifdef P3M_PRINT_TIMINGS
		require_timings = FULL;
#else
		require_timings = NONE;
#endif
		for (int i=0; i < NUM_TIMINGS; i++)
			timings[i] = 0.0;

		/* Init the derived params */
		grid_off[0] = P3M_DEFAULT_GRIDOFF;
		grid_off[1] = P3M_DEFAULT_GRIDOFF;
		grid_off[2] = P3M_DEFAULT_GRIDOFF;
		cao_cut[0] = 0.0;
		cao_cut[1] = 0.0;
		cao_cut[2] = 0.0;
		a[0] = 0.0;
		a[1] = 0.0;
		a[2] = 0.0;
		ai[0] = 0.0;
		ai[1] = 0.0;
		ai[2] = 0.0;
		additional_grid[0] = 0.0;
		additional_grid[1] = 0.0;
		additional_grid[2] = 0.0;

		/* init the P3M data */
		sum_qpart = 0;
		sum_q2 = 0.0;
		square_sum_q = 0.0;

		caf = NULL;
		cafx = NULL;
		cafy = NULL;
		cafz = NULL;
		caf_d = NULL;
		cafx_d = NULL;
		cafy_d = NULL;
		cafz_d = NULL;

		pos_shift = 0.0;
		meshift_x = NULL;
		meshift_y = NULL;
		meshift_z = NULL;

		d_op[0] = NULL;
		d_op[1] = NULL;
		d_op[2] = NULL;
		g_force = NULL;
		g_energy = NULL;

		ks_pnum = 0;

		send_grid = NULL;
		recv_grid = NULL;


		P3M_DEBUG(printf( "P3M::P3M() finished.\n"));
	}

	data_struct::~data_struct() {
		fft.free_data(rs_grid);
		fft.free_data(ks_grid);
		sfree(send_grid);
		sfree(recv_grid);
		sdelete(caf);
		sdelete(cafx);
		sdelete(cafy);
		sdelete(cafz);
		sdelete(caf_d);
		sdelete(cafx_d);
		sdelete(cafy_d);
		sdelete(cafz_d);
		sfree(g_energy);
		sfree(g_force);
		for (p3m_int i=0; i<3; i++)
			sfree(d_op[i]);
	}

}
