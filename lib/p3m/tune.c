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

/** Good mesh sizes for fftw. 
    -1 denotes end of list.
**/
static const int good_gridsize[] = 
  {4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 0,
   36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 0,
   66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128, 0,
   130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192, 0,
   196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256, 0,
   260, 264, 270, 280, 288, 294, 0,
   300, 360, 400, 480, 512, 0,
   576, 648, 729, 768, 864, 972, 1024, 0,
   1152, 1296, 1458, 1536, 0,
   1728, 1944, 2048, 0,
   2187, 2304, 2592, 0,
   2916, 3072, 0,
   3456, 0,
   3888, 0,
   -1};

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
      ifcs_p3m_compute_error_and_tune_alpha(d);
    } else if (d->tune_cao && d->tune_grid) {
      /* TUNE CAO AND GRID */
      P3M_INFO(printf( "    Tuning grid and cao.\n"));
      d->cao = cao_max;
      /* find smallest possible grid */
      P3M_DEBUG(printf( "    Finding the minimal grid at maximal cao...\n"));
      fcs_int lower_ix = 0;
      fcs_int upper_ix = 0;
      while (1) {
        /* advance to next good gridsize step */
        while (good_gridsize[upper_ix] > 0) upper_ix++;
        fcs_int grid1d = good_gridsize[upper_ix-1];
        if (grid1d == -1) break;

        /* test the gridsize */
	d->grid[0] = grid1d;
	d->grid[1] = grid1d;
	d->grid[2] = grid1d;
	P3M_INFO(printf( "      Testing grid=%d...\n", grid1d));
	ifcs_p3m_compute_error_and_tune_alpha(d);
	if (d->error < d->tolerance_field) {
          upper_ix--;
          break;
        }

        lower_ix = upper_ix+1;
        upper_ix++;
      }

      if (d->error < d->tolerance_field) {
        /* now the right gridsize is between lower_ix and upper_ix */
        /* bisect to find the right size */
        while (lower_ix+1 < upper_ix) {
          fcs_int test_ix = (lower_ix+upper_ix)/2;
          fcs_int grid1d = good_gridsize[test_ix];

          /* test the gridsize */
          d->grid[0] = grid1d;
          d->grid[1] = grid1d;
          d->grid[2] = grid1d;
          P3M_INFO(printf( "      Testing grid=%d...\n",grid1d));
          fcs_float error_before = d->error;
          fcs_float rs_error_before = d->rs_error;
          fcs_float ks_error_before = d->ks_error;
          ifcs_p3m_compute_error_and_tune_alpha(d);
          if (d->error < d->tolerance_field) {
            // parameters achieve error
            upper_ix = test_ix;
          } else {
            // parameters do not achieve error
            lower_ix = test_ix;
            // return old errors
            d->error = error_before;
            d->ks_error = ks_error_before;
            d->rs_error = rs_error_before;
          }
        }
        /* now the right size is at upper_ix */
        fcs_int grid1d = good_gridsize[upper_ix];
        d->grid[0] = grid1d;
        d->grid[1] = grid1d;
        d->grid[2] = grid1d;
	P3M_DEBUG(printf( "    => minimal grid=(%d, %d, %d)\n", \
			  d->grid[0], d->grid[1], d->grid[2]));
      
	/* find smallest possible value of cao */
	P3M_DEBUG(printf( "    Finding minimal cao...\n"));
	for (d->cao = cao_min; d->cao <= cao_max; d->cao++) {
	  P3M_INFO(printf( "      Testing cao=%d\n", d->cao));
	  ifcs_p3m_compute_error_and_tune_alpha(d);
	  if (d->error < d->tolerance_field) break;
	}
	if (d->error < d->tolerance_field)
	  P3M_DEBUG(printf( "    => minimal cao=%d\n", d->cao));
      }
    }

    // send the final parameter set and end slave loop
    ifcs_p3m_param_broadcast(d);

  } else {
    // start the slave loop and wait for parallel requests
    ifcs_p3m_param_broadcast_slave(d);
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





