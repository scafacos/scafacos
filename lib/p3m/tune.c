/*
  Copyright (C) 2011,2012,2013 Olaf Lenz
  
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

#include "tune.h"

#include "FCSCommon.h"
#include "tune_broadcast.h"
#include "error_estimate.h"
#include "timing.h"
#include "run.h"
#include "utils.h"
#include "prepare.h"
#include <stdlib.h>
#include <stdio.h>

/***************************************************/
/* TYPES AND CONSTANTS */
/***************************************************/
/** Good mesh sizes for fftw 
 */
static const fcs_int good_gridsize[] = 
  {0, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32,
   36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64,
   66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128,
   130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192,
   196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256,
   260, 264, 270, 280, 288, 294,
   300, 360, 400, 480, 512,
   576, 648, 729, 768, 864, 972, 1024,
   1152, 1296, 1458, 1536,
   1728, 1944, 2048,
   2187, 2304, 2592,
   2916, 3072,
   3456,
   3888};

/** Steps to do when trying to determine smallest possible grid size. */
static const fcs_int step_good_gridsize[] =
  { 0, 16, 27, 45, 59, 73, 79, 84, 91, 95, 99, 102, 104, 105 } ;
static const fcs_int num_steps_good_gridsize = 14;

struct tune_params_t;
typedef struct tune_params_t {
  /** cutoff radius */
  fcs_float r_cut;
  /** charge assignment order ([0,P3M_MAX_CAO]). */
  fcs_int cao;
  /** number of grid points per coordinate direction (>0). */
  fcs_int grid[3];
  /** Ewald splitting parameter */
  fcs_float alpha;

  /** Errors */
  fcs_float rs_error, ks_error, error;
  /** Timings */
  fcs_float timing, timing_near, timing_far;

  /** pointer to next param set */
  struct tune_params_t *next_params;
} tune_params;

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
static FCSResult
ifcs_p3m_tune_broadcast_master(ifcs_p3m_data_struct *d, 
                               fcs_int num_particles, fcs_int max_num_particles,
                               fcs_float *positions, fcs_float *charges);
static FCSResult
ifcs_p3m_tune_r_cut_cao_grid(ifcs_p3m_data_struct *d, 
                             fcs_int num_particles, fcs_int max_num_particles,
                             fcs_float *positions, fcs_float *charges,
                             tune_params **params);
static FCSResult
ifcs_p3m_tune_cao_grid(ifcs_p3m_data_struct *d, 
                       fcs_int num_particles, fcs_int max_num_particles,
                       fcs_float *positions, fcs_float *charges,
                       tune_params **params);
static FCSResult
ifcs_p3m_tune_grid(ifcs_p3m_data_struct *d, 
                   fcs_int num_particles, fcs_int max_num_particles,
                   fcs_float *positions, fcs_float *charges,
                   tune_params **params);
static FCSResult
ifcs_p3m_time_params(ifcs_p3m_data_struct *d,
                     fcs_int num_particles, fcs_int max_particles,
                     fcs_float *positions, fcs_float *charges,
                     tune_params **params_to_try);

static void  
ifcs_p3m_count_charges
(ifcs_p3m_data_struct *d, fcs_int num_particles, fcs_float *charges);

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
    P3M_INFO(printf( "  Number of charges changed, retuning is needed.\n"));
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
    
    /* check whether cutoff is larger than half a box length */
    if ((d->r_cut > 0.5*d->box_l[0]) ||
	(d->r_cut > 0.5*d->box_l[1]) ||
	(d->r_cut > 0.5*d->box_l[2]))
      return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, 
                              "r_cut is larger than half a system box length.");

    /* check whether cutoff is larger than domain size */
  }

  FCSResult result = NULL;
  if (d->comm.rank == 0)
    result = ifcs_p3m_tune_broadcast_master(d, num_particles, max_particles, positions, charges);
  else 
    result = ifcs_p3m_tune_broadcast_slave(d, num_particles, max_particles, positions, charges);
  if (result != NULL) {
    P3M_INFO(printf( "  Tuning failed.\n"));
  } else {
    P3M_INFO(printf( "  Tuning was successful.\n"));

    /* At the end of retuning, prepare the method */
    ifcs_p3m_prepare(d, max_particles);

    /* mark that the method was retuned */
    d->needs_retune = 0;
  }

  P3M_INFO(printf( "ifcs_p3m_tune() finished.\n"));
  return result;
}


static FCSResult
ifcs_p3m_tune_broadcast_master(ifcs_p3m_data_struct *d, 
                               fcs_int num_particles, fcs_int max_num_particles,
                               fcs_float *positions, fcs_float *charges) {
  P3M_INFO(printf("  Tuning P3M to p3m_tolerance_field=" FFLOATE "\n",  \
                  d->tolerance_field));

  tune_params *params = NULL;
  FCSResult result = NULL;
  
  result = ifcs_p3m_tune_r_cut_cao_grid(d, num_particles, max_num_particles, 
                                        positions, charges, &params);
  
  /* Set and broadcast the final parameters. */
  d->r_cut = params->r_cut;
  d->alpha = params->alpha;
  d->grid[0] = params->grid[0];
  d->grid[1] = params->grid[1];
  d->grid[2] = params->grid[2];
  d->cao = params->cao;
  ifcs_p3m_tune_broadcast_command(d, CMD_FINISHED);

  /* free the final param set */
  free(params);

  return result;
}

static FCSResult
ifcs_p3m_tune_r_cut_cao_grid(ifcs_p3m_data_struct *d, 
                             fcs_int num_particles, fcs_int max_num_particles,
                             fcs_float *positions, fcs_float *charges,
                             tune_params **params) {
  *params = malloc(sizeof(tune_params));
  (*params)->next_params = NULL;
  
  if (d->tune_r_cut) {
    /* compute the average distance between two charges  */
    fcs_float avg_dist = 
      pow((d->box_l[0]*d->box_l[1]*d->box_l[2]) 
          / d->sum_qpart, 0.33333);

    /* FIRST ESTIMATE */
    /* set the initial r_cut to 3 times the average distance between
       charges */
    (*params)->r_cut = 3.0 * avg_dist;
    
    /* tune r_cut to half the box length */
    if (0.5*d->box_l[1]-d->skin < (*params)->r_cut)
      (*params)->r_cut = 0.5*d->box_l[0] - d->skin;
    if (0.5*d->box_l[1]-d->skin < (*params)->r_cut)
      (*params)->r_cut = 0.5*d->box_l[1] - d->skin;
    if (0.5*d->box_l[2]-d->skin < (*params)->r_cut)
      (*params)->r_cut = 0.5*d->box_l[2] - d->skin;

    P3M_INFO(printf( "    r_cut=" FFLOAT " (first estimate)\n", (*params)->r_cut));
    FCSResult result = 
      ifcs_p3m_tune_cao_grid(d, num_particles, max_num_particles, 
                             positions, charges, params);
    if (result) return result;

    /* SECOND ESTIMATE */
    /* use the fact that we know that timing_near scales like r_cut**3
       and timing_far like r_cut**(-3) to get the second estimate */

    /* @ToDo: check whether it converges */
    fcs_float rel_timing_diff = 
      fabs((*params)->timing_near - (*params)->timing_far) / 
      ((*params)->timing_near + (*params)->timing_far);
    P3M_INFO(printf( "    rel_timing_diff=" FFLOAT, rel_timing_diff));

    /* /\* @ToDo: Replace with constant *\/ */
    /* if (rel_timing_diff < 0.1) { */
    /*   P3M_INFO(printf( " => found r_cut\n")); */
    /*   return NULL; */
    /* } else P3M_INFO(printf( " => tune r_cut\n")); */

    fcs_float rcut3 = (*params)->r_cut * (*params)->r_cut * (*params)->r_cut;
    fcs_float c_near = (*params)->timing_near/rcut3;
    fcs_float c_far = (*params)->timing_far*rcut3;
    fcs_float rcut_new = pow(c_far/c_near, 1./6.);
    
    (*params)->r_cut = rcut_new;
    P3M_INFO(printf( "    r_cut=" FFLOAT " (second estimate)\n", (*params)->r_cut));
    result = 
      ifcs_p3m_tune_cao_grid(d, num_particles, max_num_particles, 
                             positions, charges, params);
    if (result) return result;

    rel_timing_diff = 
      fabs((*params)->timing_near - (*params)->timing_far) / 
      ((*params)->timing_near + (*params)->timing_far);
    P3M_INFO(printf("    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));
    P3M_INFO(printf("    Finished tuning.\n"));

    return NULL;
  } else {
    (*params)->r_cut = d->r_cut;
    P3M_INFO(printf( "    r_cut=" FFLOAT " (fixed)\n", (*params)->r_cut));
    return ifcs_p3m_tune_cao_grid(d, num_particles, max_num_particles, 
                                  positions, charges, params);
  }
   
}

static FCSResult
ifcs_p3m_tune_cao_grid(ifcs_p3m_data_struct *d, 
                       fcs_int num_particles, fcs_int max_num_particles,
                       fcs_float *positions, fcs_float *charges,
                       tune_params **params) {
#ifdef P3M_AD
  const fcs_int cao_min = 2;
#else
  const fcs_int cao_min = 1;
#endif

  if (d->tune_cao) {
    for (tune_params *p = *params; p != NULL; p = p->next_params) {
      p->cao = P3M_MAX_CAO;
      P3M_INFO(printf("    r_cut=" FFLOAT ", cao={ ", p->r_cut));
      for (fcs_int cao = P3M_MAX_CAO-1; cao >= cao_min; cao--) {
        // Insert new param set
        tune_params *pnew = malloc(sizeof(tune_params));
        pnew->r_cut = p->r_cut;
        pnew->cao = cao;
        pnew->next_params = p->next_params;
        p->next_params = pnew;
        p = pnew;
        P3M_INFO(printf(FINT " ", cao));
      }
      P3M_INFO(printf("}\n"));
    }
  } else {
    // Set cao in param set
    for (tune_params *p = *params; p != NULL; p = p->next_params)
      p->cao = d->cao;
    P3M_INFO(printf( "    cao=" FINT " (fixed)\n", (*params)->cao));
  }

  return ifcs_p3m_tune_grid(d, num_particles, max_num_particles, 
                            positions, charges,
                            params);
}

/* params with identical r_cut should be adjacent and have decreasing cao */
static FCSResult
ifcs_p3m_tune_grid(ifcs_p3m_data_struct *d, 
                   fcs_int num_particles, fcs_int max_num_particles,
                   fcs_float *positions, fcs_float *charges,
                   tune_params **params) {
  if (d->tune_grid) {
    fcs_int step_ix = 0;
    // store the minimal grid size seen so far
    fcs_int min_grid1d = good_gridsize[step_good_gridsize[num_steps_good_gridsize-1]];
    // store the last param set, so that we can compare it to the current one
    tune_params **last_p = NULL;
    tune_params **p = params;
    while (*p != NULL) {

      // Find smallest possible grid for this set of parameters
      // Test whether accuracy can be achieved with given parameters
      d->r_cut = (*p)->r_cut;
      d->cao = (*p)->cao;

      // reset the step only if r_cut has changed
      if (last_p == NULL || (*last_p)->r_cut != (*p)->r_cut)
        step_ix = 0;
      
      fcs_int upper_ix;
      // test this step
      P3M_INFO(printf("    Trying to find grid for r_cut=" FFLOAT ", "  \
                      "cao=" FINT "\n",                                 \
                      d->r_cut, d->cao));
      do {
        step_ix++;
        if (step_ix >= num_steps_good_gridsize) break;
        upper_ix = step_good_gridsize[step_ix];
        fcs_int grid1d = good_gridsize[upper_ix];
        d->grid[0] = grid1d;
        d->grid[1] = grid1d;
        d->grid[2] = grid1d;
        P3M_DEBUG(printf("      rough grid=" FINT "\n", grid1d));
        ifcs_p3m_compute_error_and_tune_alpha(d);
        (*p)->error = d->error;
        (*p)->rs_error = d->rs_error;
        (*p)->ks_error = d->ks_error;
        P3M_DEBUG(printf("        => alpha=" FFLOAT ", "                \
                         "error=" FFLOATE "\n", d->alpha, d->error));
      } while (d->error > d->tolerance_field);

      // reached largest possible grid, remove this parameter set
      if (step_ix >= num_steps_good_gridsize) {
        P3M_INFO(printf("    Too large grid size, skipping parameter set.\n"));
        tune_params *ptmp = *p;
        *p = (*p)->next_params;
        free(ptmp);
        // also remove all params with the same r_cut
        while (*p != NULL && (*p)->r_cut == d->r_cut) {
          ptmp = *p;
          *p = (*p)->next_params;
          free(ptmp);
        }
        continue;
      }

      // error would be small enough at upper_ix, but not at lower_ix,
      // so bisect to find optimal gridsize
      fcs_int lower_ix = step_good_gridsize[step_ix-1];
      while (lower_ix+1 < upper_ix) {
        fcs_int test_ix = (lower_ix+upper_ix)/2;
        fcs_int grid1d = good_gridsize[test_ix];
        d->grid[0] = grid1d;
        d->grid[1] = grid1d;
        d->grid[2] = grid1d;
        P3M_DEBUG(printf("      fine grid=" FINT "\n", grid1d));
        ifcs_p3m_compute_error_and_tune_alpha(d);
        P3M_DEBUG(printf("          => alpha=" FFLOAT ", "              \
                         "error=" FFLOATE "\n", d->alpha, d->error));
        if (d->error < d->tolerance_field) {
          // parameters achieve error
          upper_ix = test_ix;
          (*p)->error = d->error;
          (*p)->rs_error = d->rs_error;
          (*p)->ks_error = d->ks_error;
        } else {
          // parameters do not achieve error
          lower_ix = test_ix;
        }
      }

      // now the right size is at upper_ix
      fcs_int grid1d = good_gridsize[upper_ix];

      // store the new grid size and alpha
      if (min_grid1d > grid1d) min_grid1d = grid1d;

      (*p)->alpha = d->alpha;
      (*p)->grid[0] = grid1d;
      (*p)->grid[1] = grid1d;
      (*p)->grid[2] = grid1d;
      P3M_INFO(printf( "      => grid=" F3INT ", "                      \
                       "alpha=" FFLOAT ", "                             \
                       "error=" FFLOATE "\n",                           \
                       (*p)->grid[0], (*p)->grid[1], (*p)->grid[2],     \
                       (*p)->alpha, (*p)->error));
      
      // decrease step_ix so that the same step_ix is tested for the
      // next param set
      step_ix--;

      // compare grid size to previous data set
      // if it is larger than any previous size + P3M_MAX_GRID_DIFF, remove it
      if (min_grid1d + P3M_MAX_GRID_DIFF < grid1d) {
        P3M_INFO(printf("      grid too large => removing data set\n"));
        tune_params *tmp = *p;
        *p = (*p)->next_params;
        free(tmp);
        continue;
      } 
      // compare grid size and cao to previous data set
      // if previous param set has identical r_cut but larger or equal cao and/or grid, 
      // remove it
      else if (last_p != NULL && (*last_p)->r_cut == (*p)->r_cut 
               && (*last_p)->cao >= (*p)->cao && (*last_p)->grid[0] >= (*p)->grid[0]) {
        P3M_INFO(printf("      better than previous => removing previous data set\n"));
        tune_params *tmp = *last_p;
        *last_p = *p;
        free(tmp);
      } else {
        last_p = p;
      }

      // advance to next parameter set
      p = &((*p)->next_params);
    }

  } else {

    // fixed grid
    tune_params **p = params;

    while (*p != NULL) {
      // test whether accuracy can be achieved with given parameters
      // d->grid already contains the wanted grid
      d->cao = (*p)->cao;
      d->r_cut = (*p)->r_cut;
      P3M_INFO(printf("    r_cut=" FFLOAT ", "                   \
                      "cao=" FINT ", "                           \
                      "grid=" F3INT " (fixed)\n",                \
                      d->r_cut, d->cao,                          \
                      d->grid[0], d->grid[1], d->grid[2]));
      ifcs_p3m_compute_error_and_tune_alpha(d);
      if (d->error < d->tolerance_field) {
        // error is small enough for this parameter set, so keep it
        P3M_INFO(printf("    error is small enough, keeping params\n"));
        (*p)->grid[0] = d->grid[0];
        (*p)->grid[1] = d->grid[1];
        (*p)->grid[2] = d->grid[2];
        (*p)->alpha = d->alpha;
        // advance to next param set
        p = &((*p)->next_params);
      } else {
        // otherwise remove this parameter set
        P3M_INFO(printf("    error too large, removing params\n"));
        tune_params *ptmp = *p;
        // change pointer to next param set
        *p = (*p)->next_params;
        free(ptmp);
      }
    }
  }

  return ifcs_p3m_time_params(d, num_particles, max_num_particles, 
                              positions, charges,
                              params);
}

static FCSResult
ifcs_p3m_time_params(ifcs_p3m_data_struct *d,
                     fcs_int num_particles, fcs_int max_particles,
                     fcs_float *positions, fcs_float *charges,
                     tune_params **params) {
  const char* fnc_name = "ifcs_p3m_time_params";

  /* Now time the different parameter sets */
  tune_params *best_params = NULL;
  double best_timing = 1.e100;

#ifdef FCS_ENABLE_INFO
  fcs_int num_params=0;
  for (tune_params *current_params = *params;
       current_params != NULL;  
       current_params = current_params->next_params)
    num_params++;
  printf("Timing %d param sets...\n", num_params);
#endif

  for (tune_params *current_params = *params;
       current_params != NULL;  
       current_params = current_params->next_params) {
    /* use the parameters */
    d->r_cut = current_params->r_cut;
    d->alpha = current_params->alpha;
    d->grid[0] = current_params->grid[0];
    d->grid[1] = current_params->grid[1];
    d->grid[2] = current_params->grid[2];
    d->cao = current_params->cao;
    
    ifcs_p3m_timing(d, num_particles, max_particles, positions, charges);
    current_params->timing = d->timings[TIMING];
    current_params->timing_near = d->timings[TIMING_NEAR];
    current_params->timing_far = d->timings[TIMING_FAR];
    P3M_INFO(printf( "  Timing r_cut=" FFLOAT ", "                      \
                     "alpha=" FFLOAT ", "                               \
                     "grid=" F3INT ", "                                 \
                     "cao=" FINT " "                                    \
                     "=> timing=" FFLOAT " "                            \
                     "(" FFLOAT " near, " FFLOAT " far)\n",              \
                     d->r_cut, d->alpha,                                \
                     d->grid[0], d->grid[1], d->grid[2],                \
                     d->cao, current_params->timing,                    \
                     current_params->timing_near, current_params->timing_far));

    if (current_params->timing < best_timing) {
      best_timing = current_params->timing;
      best_params = current_params;
    }
  }

  if (best_params == NULL)
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name, 
                            "Internal error: No best timing.");

  /* Free the list, keep only the best */
  tune_params *current_params = *params;
  while (current_params != NULL) {
    tune_params *next_params = current_params->next_params;
    if (current_params != best_params)
      free(current_params);
    current_params = next_params;
  }

  best_params->next_params = NULL;
  *params = best_params;

  return NULL;
}

/** Calculate number of charged particles, the sum of the squared
    charges and the squared sum of the charges. Called in parallel at
    the beginning of tuning. */
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
		   "num_charged_particles=" FINT ", "                   \
                   "sum_squared_charges=" FFLOAT ", "                   \
                   "sum_charges=" FFLOAT "\n",                          \
		   d->sum_qpart, d->sum_q2, sqrt(d->square_sum_q)));
}
