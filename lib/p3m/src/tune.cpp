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

#include "tune_broadcast.hpp"
#include "error_estimate.hpp"
#include "timing.hpp"
#include "utils.hpp"
#include "prepare.hpp"
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <list>

namespace ScaFaCoS {
  namespace P3M {
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

    struct tune_params {
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
    };
    typedef std::list<tune_params*> tune_params_l;

    /***************************************************/
    /* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
    /***************************************************/
    tune_params* 
    tune_far(data_struct *d,
             fcs_int num_particles, fcs_float *positions, fcs_float *charges,
             fcs_float r_cut);
    
    void
    tune_broadcast_master(data_struct *d, fcs_int num_particles,
                          fcs_float *positions, fcs_float *charges);
    
    tune_params*
    tune_alpha_cao_grid(data_struct *d, 
                        fcs_int num_particles, fcs_float *positions, fcs_float *charges, 
                        fcs_float r_cut);
      
    tune_params* 
    tune_cao_grid(data_struct *d, 
                  fcs_int num_particles, fcs_float *positions, fcs_float *charges,
                  fcs_float r_cut, fcs_float alpha);

    tune_params*
    tune_grid(data_struct *d, fcs_int num_particles,
              fcs_float *positions, fcs_float *charges,
              tune_params_l &params_to_try);

    tune_params*
    time_params(data_struct *d, 
                fcs_int num_particles, fcs_float *positions, fcs_float *charges,
                tune_params_l &params_to_try);
    
    void  
    count_charges(data_struct *d, 
                  fcs_int num_particles, fcs_float *charges);

    /***************************************************/
    /* IMPLEMENTATION */
    /***************************************************/
    void
    tune(data_struct *d, fcs_int num_particles,
         fcs_float *positions, fcs_float *charges) {

      /* Prepare the communicator before tuning */
      comm_prepare(&d->comm, d->box_l);

      /* Count the charges */
      fcs_float sum_q2_before = d->sum_q2;
      count_charges(d, num_particles, charges);
      
      /* Retune if the number of charges has changed */
      if (!d->needs_retune && !(float_is_equal(sum_q2_before, d->sum_q2))) {
        P3M_INFO(printf( "  Number of charges changed, retuning is needed.\n"));
        d->needs_retune = 1;
      }
      
      /* Do not retune if there are no charges */
      if (float_is_zero(d->sum_q2)) {
        P3M_INFO(printf( "  No charges in the system.\n"));
        P3M_INFO(printf( "  Retuning is not required.\n"));
        P3M_INFO(printf( "tune_far() finished.\n"));
        return;
      }
      
      /* Exit if retuning is unnecessary */
      if (!d->needs_retune) {
        P3M_INFO(printf( "  Retuning is not required.\n"));
        P3M_INFO(printf( "tune_far() finished.\n"));
        return;
      }
      
      P3M_INFO(printf( "  Retuning is required.\n"));

      if (d->tune_r_cut) {

        /* compute the average distance between two charges  */
        fcs_float avg_dist = 
          pow((d->box_l[0]*d->box_l[1]*d->box_l[2]) 
              / d->sum_qpart, 0.33333);
        
        /* FIRST ESTIMATE */
        /* set the initial r_cut to 3 times the average distance between
           charges */
        fcs_float r_cut = 3.0 * avg_dist;
        
        /* tune r_cut to half the box length */
        if (0.5*d->box_l[1]-d->skin < r_cut)
          r_cut = 0.5*d->box_l[0] - d->skin;
        if (0.5*d->box_l[1]-d->skin < r_cut)
          r_cut = 0.5*d->box_l[1] - d->skin;
        if (0.5*d->box_l[2]-d->skin < r_cut)
          r_cut = 0.5*d->box_l[2] - d->skin;

        P3M_INFO(printf( "    r_cut=" FFLOAT " (first estimate)\n", r_cut));

        // @todo get near timing
        tune_params *p = 
          tune_far(d, num_particles, positions, charges, r_cut);
        
        /* SECOND ESTIMATE */
        /* use the fact that we know that timing_near scales like r_cut**3
           and timing_far like r_cut**(-3) to get the second estimate */

        fcs_float rel_timing_diff = 
          fabs(p->timing_near - p->timing_far) / 
          (p->timing_near + p->timing_far);
        P3M_INFO(printf( "    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));

        fcs_float rcut3 = pow(r_cut, 3);
        fcs_float c_near = p->timing_near/rcut3;
        fcs_float c_far = p->timing_far*rcut3;
        fcs_float rcut_new = pow(c_far/c_near, 1./6.);
    
        r_cut = rcut_new;
        P3M_INFO(printf( "    r_cut=" FFLOAT " (second estimate)\n", r_cut));

        // @todo get near timing
        p = tune_far(d, num_particles, positions, charges, r_cut);
        
        rel_timing_diff = 
          fabs(p->timing_near - p->timing_far) / 
          (p->timing_near + p->timing_far);
        P3M_INFO(printf("    rel_timing_diff=" FFLOAT "\n", rel_timing_diff));
        P3M_INFO(printf("    Finished tuning.\n"));
      } else {
        P3M_INFO(printf( "    r_cut=" FFLOAT " (fixed)\n", d->r_cut));
        tune_far(d, num_particles, positions, charges, d->r_cut);
      }

      /* mark that the method was retuned */
      d->needs_retune = 0;
    }

    tune_params* 
    tune_far(data_struct *d,
             fcs_int num_particles, fcs_float *positions, fcs_float *charges,
             fcs_float r_cut) {
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
      P3M_INFO(printf( "tune_far() started...\n"));

      /* check whether the input parameters are sane */
      if (r_cut < 0.0)
        throw std::domain_error("r_cut is negative!");
    
      if (float_is_zero(r_cut)) 
        throw std::domain_error("r_cut is too small!");
    
      /* check whether cutoff is larger than half a box length */
      if ((r_cut > 0.5*d->box_l[0]) ||
          (r_cut > 0.5*d->box_l[1]) ||
          (r_cut > 0.5*d->box_l[2]))
        throw std::domain_error("r_cut is larger than half a system box length.");

      /* check whether cutoff is larger than domain size */

      tune_params *final_params = NULL;
      try {
        d->r_cut = r_cut;
        if (d->comm.rank == 0) {
          final_params = 
            tune_alpha_cao_grid(d, num_particles, positions, charges, r_cut);
          
          // @todo: move to tune function
          /* Set and broadcast the final parameters. */
          d->r_cut = r_cut;
          d->alpha = final_params->alpha;
          d->grid[0] = final_params->grid[0];
          d->grid[1] = final_params->grid[1];
          d->grid[2] = final_params->grid[2];
          d->cao = final_params->cao;
          tune_broadcast_command(d, CMD_FINISHED);
        } else {
          tune_broadcast_slave(d, num_particles, positions, charges);
        }
      } catch (...) {
        P3M_INFO(printf( "  Tuning failed.\n"));
        P3M_INFO(printf( "tune_far() finished.\n"));
        throw;
      }

      P3M_INFO(printf( "  Tuning was successful.\n"));

      /* At the end of retuning, prepare the method */
      prepare(d);
      
      
      P3M_INFO(printf( "tune_far() finished.\n"));

      return final_params;
    }


    /* Tune alpha */
    tune_params*
    tune_alpha_cao_grid(data_struct *d, 
                        fcs_int num_particles, fcs_float *positions, fcs_float *charges, 
                        fcs_float r_cut) {
      if (d->tune_alpha) {
        determine_good_alpha(d);
        P3M_INFO(printf("    => alpha=" FFLOAT "\n", d->alpha));
      } else {
        P3M_INFO(printf("    alpha=" FFLOAT " (fixed)\n", d->alpha));
      }

      return tune_cao_grid(d, num_particles, positions, charges, r_cut, d->alpha);
    }

    /* Tune cao */
    tune_params* 
    tune_cao_grid(data_struct *d, 
                  fcs_int num_particles, fcs_float *positions, fcs_float *charges,
                  fcs_float r_cut, fcs_float alpha) {
#ifdef P3M_AD
      const fcs_int min_cao = 2;
#else
      const fcs_int min_cao = 1;
#endif

      tune_params_l params_to_try;

      if (d->tune_cao) {
        P3M_INFO(printf("    Testing cao={ "));
        for (fcs_int cao = CAF::max_cao; cao >= min_cao; --cao) {
          tune_params *p = new tune_params;
          p->cao = cao;
          p->alpha = alpha;
          params_to_try.push_back(p);
          P3M_INFO(printf(FINT " ", cao));
        }
        P3M_INFO(printf("}\n"));
      } else {
        tune_params *p = new tune_params;
        p->cao = d->cao;
        p->alpha = alpha;
        P3M_INFO(printf( "    cao=" FINT " (fixed)\n", d->cao));
      }

      return tune_grid(d, num_particles, positions, charges, params_to_try);
    }

    /* params should have decreasing cao */
    tune_params*
    tune_grid(data_struct *d, fcs_int num_particles,
              fcs_float *positions, fcs_float *charges,
              tune_params_l &params_to_try) {
      if (d->tune_grid) {
        fcs_int step_ix = 0;
        // store the minimal grid size seen so far
        fcs_int min_grid1d = good_gridsize[step_good_gridsize[num_steps_good_gridsize-1]];

        for (tune_params_l::iterator p = params_to_try.begin(); 
             p != params_to_try.end(); ++p) {
          // Find smallest possible grid for this set of parameters
          // Test whether accuracy can be achieved with given parameters
          d->cao = (*p)->cao;
      
          fcs_int upper_ix;
          // test this step
          P3M_INFO(printf("    Trying to find grid for r_cut=" FFLOAT ", " \
                          "alpha=" FFLOAT ", "                          \
                          "cao=" FINT "\n",                             \
                          d->r_cut, d->alpha, d->cao));
          do {
            step_ix++;
            if (step_ix >= num_steps_good_gridsize) break;
            upper_ix = step_good_gridsize[step_ix];
            fcs_int grid1d = good_gridsize[upper_ix];
            d->grid[0] = grid1d;
            d->grid[1] = grid1d;
            d->grid[2] = grid1d;
            P3M_DEBUG(printf("      rough grid=" FINT "\n", grid1d));
            compute_error_estimate(d);
            (*p)->error = d->error;
            (*p)->rs_error = d->rs_error;
            (*p)->ks_error = d->ks_error;
            P3M_DEBUG(printf("        => error=" FFLOATE "\n", d->error));
          } while (d->error > d->tolerance_field);

          // reached largest possible grid, remove the rest of the parameter sets
          if (step_ix >= num_steps_good_gridsize) {
            P3M_INFO(printf("    Too large grid size, skipping rest of parameter sets.\n"));
            for (tune_params_l::iterator pit = params_to_try.end(); 
                 pit != p; --pit) {
              delete[] *pit;
              params_to_try.erase(pit);
            }
            break;
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
            compute_error_estimate(d);
            P3M_DEBUG(printf("          => error=" FFLOATE "\n", d->error));
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

          (*p)->grid[0] = grid1d;
          (*p)->grid[1] = grid1d;
          (*p)->grid[2] = grid1d;
          P3M_INFO(printf( "      => grid=" F3INT ", "                      \
                           "error=" FFLOATE "\n",                           \
                           (*p)->grid[0], (*p)->grid[1], (*p)->grid[2],     \
                           (*p)->error));
      
          // decrease step_ix so that the same step_ix is tested for the
          // next param set
          step_ix--;

          // compare grid size to previous data set
          // if it is larger than any previous size + P3M_MAX_GRID_DIFF, skip it
          if (min_grid1d + P3M_MAX_GRID_DIFF < grid1d) {
            P3M_INFO(printf("      grid too large => removing data set\n"));
            // remove the rest of the params
            for (tune_params_l::iterator pit = p;
                 pit != params_to_try.end(); ++pit) {
              delete[] *pit;
              *pit = NULL;
            }
            break;
          }

          if (p != params_to_try.begin()) {
            tune_params_l::iterator prev = p;
            prev--;
            if ((*prev)->grid[0] >= (*p)->grid[0]) {
              P3M_INFO(printf("      better than previous => removing previous data set.\n"));
              delete[] *prev;
              params_to_try.erase(prev);
            }
          }
        }

        // remove empty param sets
        params_to_try.remove(NULL);


      } else {

        for (tune_params_l::iterator p = params_to_try.begin(); 
             p != params_to_try.end(); ++p) {
          // test whether accuracy can be achieved with given parameters
          // d->grid already contains the wanted grid
          d->cao = (*p)->cao;
          P3M_INFO(printf("    r_cut=" FFLOAT ", "                   \
                          "cao=" FINT ", "                           \
                          "grid=" F3INT " (fixed)\n",                \
                          d->r_cut, d->cao,                          \
                          d->grid[0], d->grid[1], d->grid[2]));
          compute_error_estimate(d);

          if (d->error < d->tolerance_field) {
            // error is small enough for this parameter set, so keep it
            P3M_INFO(printf("    error (%le) < tolerance (%le), keeping params\n", \
                            d->error, d->tolerance_field));
            (*p)->grid[0] = d->grid[0];
            (*p)->grid[1] = d->grid[1];
            (*p)->grid[2] = d->grid[2];
            (*p)->alpha = d->alpha;
          } else {
            // otherwise remove this parameter set
            P3M_INFO(printf("    error (%le) > tolerance (%le), removing params\n", \
                            d->error, d->tolerance_field));
            delete[] *p;
            *p = NULL;
          }
        }

        // really remove the removed param sets
        params_to_try.remove(NULL);
      }

      return time_params(d, num_particles, positions, charges, params_to_try);
    }

    tune_params*
    time_params(data_struct *d, 
                fcs_int num_particles, fcs_float *positions, fcs_float *charges,
                tune_params_l &params_to_try) {
      /* Now time the different parameter sets */
      tune_params *best_params = NULL;
      double best_timing = 1.e100;
      
#ifdef FCS_ENABLE_INFO
      printf("Timing %ld param sets...\n", params_to_try.size());
#endif

      for (tune_params_l::iterator it = params_to_try.begin();
           it != params_to_try.end();  
           it++) {
        tune_params &p = **it;
        /* use the parameters */
        d->alpha = p.alpha;
        d->grid[0] = p.grid[0];
        d->grid[1] = p.grid[1];
        d->grid[2] = p.grid[2];
        d->cao = p.cao;
    
        timing(d, num_particles, positions, charges);
        p.timing = d->timings[TIMING];
        p.timing_near = d->timings[TIMING_NEAR];
        p.timing_far = d->timings[TIMING_FAR];
        P3M_INFO(printf( "  Timing r_cut=" FFLOAT ", "                  \
                         "alpha=" FFLOAT ", "                           \
                         "grid=" F3INT ", "                             \
                         "cao=" FINT " "                                \
                         "=> timing=" FFLOAT " "                        \
                         "(" FFLOAT " near, " FFLOAT " far)\n",         \
                         d->r_cut, d->alpha,                            \
                         d->grid[0], d->grid[1], d->grid[2],            \
                         d->cao, p.timing,                \
                         p.timing_near, p.timing_far));
        
        if (p.timing < best_timing) {
          best_timing = p.timing;
          best_params = &p;
        }
      }

      if (best_params == NULL)
        throw std::runtime_error("Internal error: Tuning could not get a best timing.");

      for (tune_params_l::iterator p = params_to_try.begin();
           p != params_to_try.end(); ++p) {
        if (*p != best_params)
          delete[] *p;
      }

      return best_params;
    }

    /** Calculate number of charged particles, the sum of the squared
        charges and the squared sum of the charges. Called in parallel at
        the beginning of tuning. */
    void 
    count_charges(data_struct *d, fcs_int num_particles, fcs_float *charges) {  
      int i;
      fcs_float node_sums[3], tot_sums[3];

      for(i=0;i<3;i++) { 
        node_sums[i]=0.0; 
        tot_sums[i]=0.0;
      }

      for (i = 0; i < num_particles; i++) {
        if (!float_is_zero(charges[i])) {
          node_sums[0] += 1.0;
          node_sums[1] += SQR(charges[i]);
          node_sums[2] += charges[i];
        }
      }
  
      MPI_Allreduce(node_sums, tot_sums, 3, FCS_MPI_FLOAT, MPI_SUM, d->comm.mpicomm);
      d->sum_qpart    = (fcs_int)(tot_sums[0]+0.1);
      d->sum_q2       = tot_sums[1];
      d->square_sum_q = SQR(tot_sums[2]);

      P3M_DEBUG(printf("  count_charges(): "			\
                       "num_charged_particles=" FINT ", "                   \
                       "sum_squared_charges=" FFLOAT ", "                   \
                       "sum_charges=" FFLOAT "\n",                          \
                       d->sum_qpart, d->sum_q2, sqrt(d->square_sum_q)));
    }
  }
}
