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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "error_estimate.h"
#include "tune.h"
#include "utils.h"
#include "prepare.h"
#include <stdlib.h>
#include <stdio.h>

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
static FCSResult 
ifcs_p3m_tuneit(ifcs_p3m_data_struct *d, 
		fcs_int num_particles,
		fcs_float *positions, 
		fcs_float *charges);
static void  
ifcs_p3m_count_charges(ifcs_p3m_data_struct *d, 
		       fcs_int num_particles, 
		       fcs_float *charges);

static fcs_float 
ifcs_p3m_get_error(ifcs_p3m_data_struct *d,
		   fcs_int grid[3], 
		   fcs_int cao,
		   fcs_float r_cut,		   
		   fcs_float *_rs_err, 
		   fcs_float *_ks_err);

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
FCSResult 
ifcs_p3m_tune(void* rd,
	      fcs_int num_particles,
	      fcs_int max_particles,
	      fcs_float *positions, 
	      fcs_float *charges) {
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  const char* fnc_name = "ifcs_p3m_tune";

  /* Distinguish between two types of parameters:
     Input params:
     * box length
     * #charges
     * cutoff
     * skin (only if comm is cartesian)
     * error (default=1e-4)
     * component flags

     Performance related params:
     * caf interpolation

     Tuned params:
     * alpha
     * cao
     * grid
     */
  P3M_INFO(printf( "ifcs_p3m_tune() started...\n"));

  /* Prepare the communicator before tuning */
  ifcs_p3m_comm_prepare(&d->comm, d->box_l);

  /* Count the charges */
  fcs_float sum_q2_before = d->sum_q2;
  ifcs_p3m_count_charges(d, num_particles, charges);

  /* Retune if the number of charges has changed */
  if (!d->needs_retune && !(fcs_float_is_equal(sum_q2_before, d->sum_q2))) {
    P3M_INFO(printf( "  Number of charges changed.\n"));
    d->needs_retune = 1;
  }

  /* Do not retune if there are no charges */
  if (fcs_float_is_zero(d->sum_q2)) {
    P3M_INFO(printf( "  No charges in the system.\n"));
    P3M_INFO(printf( "  Retuning is not required.\n"));
    P3M_INFO(printf( "ifcs_p3m_tune() finished.\n"));
    return NULL;
  }

  /* Exit if retuning is unnecessary */
  if (!d->needs_retune) {
    P3M_INFO(printf( "  Retuning is not required.\n"));
    P3M_INFO(printf( "ifcs_p3m_tune() finished.\n"));
    return NULL;
  }

  P3M_INFO(printf( "  Retuning is required.\n"));

  /* check whether the input parameters are sane */
  if (!d->tune_r_cut) {
    if (d->r_cut < 0.0) {
      return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "r_cut is negative!");
    }
    
    if (fcs_float_is_zero(d->r_cut)) {
      return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "r_cut is too small!");
    }
    
    /* * check whether cutoff is larger than half a box length */
    if ((d->r_cut > 0.5*d->box_l[0]) ||
	(d->r_cut > 0.5*d->box_l[1]) ||
	(d->r_cut > 0.5*d->box_l[2]))
      return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "r_cut is larger than half a system box length.");

    /* * check whether cutoff is larger than domain size */
  }

  FCSResult result = ifcs_p3m_tuneit(d, num_particles, positions, charges);
  if (result != NULL) return result;

  P3M_INFO(printf( "  Tuning was successful, preparing the method.\n"));

  /* when the parameters were retuned, prepare the method */
  ifcs_p3m_prepare(d, max_particles);

  /* mark that method was retuned */
  d->needs_retune = 0;

  P3M_INFO(printf( "ifcs_p3m_tune() finished.\n"));
  return NULL;
}

static FCSResult
ifcs_p3m_tuneit(ifcs_p3m_data_struct *d,
		fcs_int num_particles,
		fcs_float *positions,
		fcs_float *charges) {
  const char* fnc_name="ifcs_p3m_tuneit";
  fcs_int grid1d;
  const fcs_int grid1d_min = 4;
  const fcs_int grid1d_max = 1024;
  #ifdef P3M_AD
  const fcs_int cao_min = 2;
  #else
  const fcs_int cao_min = 1;
  #endif
  const fcs_int cao_max = 7;

  P3M_DEBUG(printf( "  ifcs_p3m_tuneit() started...\n"));

  P3M_INFO(printf("  Tuning P3M to p3m_tolerance_field=%" FCS_LMOD_FLOAT "g.\n", \
		  d->tolerance_field));
  if (d->tune_r_cut) {
    /* compute the average distance between two charges  */
    fcs_float avg_dist = pow((d->box_l[0]*d->box_l[1]*d->box_l[2]) / d->sum_qpart, 0.33333);
    /* set r_cut to 3 times the average distance between charges */
    d->r_cut = 3.0 * avg_dist;
    
    /* tune r_cut to half the box length */
    if (0.5*d->box_l[1]-d->skin < d->r_cut)
      d->r_cut = 0.5*d->box_l[0] - d->skin;
    if (0.5*d->box_l[1]-d->skin < d->r_cut)
      d->r_cut = 0.5*d->box_l[1] - d->skin;
    if (0.5*d->box_l[2]-d->skin < d->r_cut)
      d->r_cut = 0.5*d->box_l[2] - d->skin;

    P3M_INFO(printf( "    tuned r_cut to %" FCS_LMOD_FLOAT "f\n", d->r_cut));
  }

  /* On master */
  if (d->comm.rank == 0) {
    /* @todo Non-cubic case */
    /* @todo Tune only grid or only cao */
    if (!d->tune_grid && !d->tune_cao) {
      /* TUNE ONLY ALPHA */
      d->error = ifcs_p3m_get_error(d, d->grid, d->cao, d->r_cut, 
                                    &d->rs_error, &d->ks_error);
    } else if (d->tune_cao && d->tune_grid) {
      /* TUNE CAO AND GRID */
      P3M_INFO(printf( "    Tuning grid and cao.\n"));
      d->cao = cao_max;
      /* find smallest possible grid */
      P3M_DEBUG(printf( "    Finding minimal grid at maximal cao...\n"));
      for (grid1d = grid1d_min; 
	   grid1d <= grid1d_max;
	   grid1d *= 2) {
	d->grid[0] = grid1d;
	d->grid[1] = grid1d;
	d->grid[2] = grid1d;
	P3M_DEBUG(printf( "      Testing grid=%d...\n", grid1d));
	d->error = ifcs_p3m_get_error(d, d->grid, d->cao, d->r_cut, 
                                      &d->rs_error, &d->ks_error);
	if (d->error < d->tolerance_field) break;
      }

      if (d->error < d->tolerance_field) {
	P3M_DEBUG(printf( "    => minimal grid=(%d, %d, %d)\n", \
			  d->grid[0], d->grid[1], d->grid[2]));
      
	/* find smallest possible value of cao */
	P3M_DEBUG(printf( "    Finding minimal cao...\n"));
	for (d->cao = cao_min; d->cao <= cao_max; d->cao++) {
	  P3M_DEBUG(printf( "      Testing cao=%d\n", d->cao));
	  d->error = ifcs_p3m_get_error(d, d->grid, d->cao, d->r_cut, 
                                        &d->rs_error, &d->ks_error);
	  if (d->error < d->tolerance_field) break;
	}
	if (d->error < d->tolerance_field)
	  P3M_DEBUG(printf( "    => minimal cao=%d\n", d->cao));
      }
    }

    // send the final parameter set and end slave loop
    ifcs_p3m_param_broadcast(d->r_cut, d->grid, d->alpha, d->cao,
				   d->error, d->rs_error, d->ks_error,
				   d->comm.mpicomm);

  } else {
    // start the slave loop and wait for parallel requests
    ifcs_p3m_param_broadcast_slave(d->sum_qpart, d->sum_q2, 
					 d->box_l,
					 &d->r_cut, d->grid, 
					 &d->alpha, &d->cao,
					 &d->error, 
					 &d->rs_error, &d->ks_error,
					 d->comm.mpicomm);
  }

  if (d->error > d->tolerance_field) {
    char msg[255];
    sprintf(msg, 
	    "Cannot achieve required accuracy (p3m_tolerance_field=%" FCS_LMOD_FLOAT "e) for given parameters.", 
	    d->tolerance_field);
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, msg);
  }

  P3M_INFO(printf( "    r_cut=%" FCS_LMOD_FLOAT "g cao=%d grid=(%d, %d, %d) alpha=%" FCS_LMOD_FLOAT "g\n", \
		   d->r_cut, d->cao, d->grid[0], d->grid[1], d->grid[2], d->alpha));
  P3M_INFO(printf( "    error=%" FCS_LMOD_FLOAT "e rs_error=%" FCS_LMOD_FLOAT "e ks_error=%" FCS_LMOD_FLOAT "e\n", \
		   d->error, d->rs_error, d->ks_error));
  P3M_DEBUG(printf( "  ifcs_p3m_tuneit() finished.\n"));
  return NULL;
}

/** Calculate number of charged particles, the sum of the squared
    charges and the squared sum of the charges. */
static void 
ifcs_p3m_count_charges(ifcs_p3m_data_struct *d, 
		       fcs_int num_particles, fcs_float *charges) {  
  int i;
  fcs_float node_sums[3], tot_sums[3];

  for(i=0;i<3;i++) { 
    node_sums[i]=0.0; 
    tot_sums[i]=0.0;
  }

  for (i = 0; i < num_particles; i++) {
    if (!fcs_float_is_zero(charges[i])) {
      node_sums[0] += 1.0;
      node_sums[1] += SQR(charges[i]);
      node_sums[2] += charges[i];
    }
  }
  
  MPI_Allreduce(node_sums, tot_sums, 3, FCS_MPI_FLOAT, MPI_SUM, d->comm.mpicomm);
  d->sum_qpart    = (fcs_int)(tot_sums[0]+0.1);
  d->sum_q2       = tot_sums[1];
  d->square_sum_q = SQR(tot_sums[2]);

  P3M_DEBUG(printf("  ifcs_p3m_count_charges(): "			\
		   "num_charged_particles=%d sum_squared_charges=%" FCS_LMOD_FLOAT "f sum_charges=%lf\n", \
		   d->sum_qpart, d->sum_q2, sqrt(d->square_sum_q)));
}

/** Get the error for this combination of parameters and tune alpha if
    required. In fact, the real space error is tuned such that it
    contributes half of the total error, and then the Fourier space
    error is calculated. Returns the achieved error and the optimal
    alpha.
    */
static fcs_float 
ifcs_p3m_get_error(ifcs_p3m_data_struct *d,
		   fcs_int grid[3], fcs_int cao, fcs_float r_cut,
		   fcs_float *_rs_err, fcs_float *_ks_err) {
  fcs_float err;
  if (d->tune_alpha)
    ifcs_p3m_determine_good_alpha(d->sum_qpart, d->sum_q2, d->box_l, 
					r_cut, grid, cao, 
					d->tolerance_field, &d->alpha,
					&err, _rs_err, _ks_err, d->comm.mpicomm);
  else
    ifcs_p3m_compute_error_estimate(d->sum_qpart, d->sum_q2, d->box_l, 
					  r_cut, grid, d->alpha, cao, 
					  &err, _rs_err, _ks_err, d->comm.mpicomm);

  return err;
}




