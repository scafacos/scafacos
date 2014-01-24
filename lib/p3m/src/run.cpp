/*
  Copyright (C) 2013 Olaf Lenz, Florian Weik
  Copyright (C) 2011,2012 Olaf Lenz
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "utils.hpp"
#include "types.hpp"
#include "p3m.hpp"
#include "caf.hpp"
#include <sys/time.h>
#include <sys/resource.h>
#include <cstdlib>
#include <cstdio>

#include "fcs_p3m_p.h"
#include "FCSCommon.h"
#include "common/near/near.h"
#include "common/gridsort/gridsort.h"

using namespace ScaFaCoS::P3M;

namespace ScaFaCoS {
  namespace P3M {
    /***************************************************/
    /* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
    /***************************************************/
    /* domain decomposition */
    static void
    domain_decompose(data_struct *d, fcs_gridsort_t *gridsort,
                     fcs_int _num_particles, 
                     fcs_float *_positions, fcs_float *_charges,
                     fcs_int *num_real_particles,
                     fcs_float **positions, fcs_float **charges,
                     fcs_gridsort_index_t **indices, 
                     fcs_int *num_ghost_particles,
                     fcs_float **ghost_positions, fcs_float **ghost_charges,
                     fcs_gridsort_index_t **ghost_indices
                     );

    /* callback function for near field computations */
    static inline void 
    compute_near(const void *param, fcs_float dist, fcs_float *f, fcs_float *p)
    {
      fcs_float alpha = *((fcs_float *) param);
      fcs_p3m_compute_near(alpha, dist, p, f);
    }
    
    /* callback function for performing a whole loop of near field computations (using compute_near) */
    FCS_NEAR_LOOP_FP(compute_near_loop, compute_near);
    
    /* charge assignment */
    static void assign_charges(data_struct* d,
                               fcs_float *data,
                               fcs_int num_particles,
                               fcs_float *positions, 
                               fcs_float *charges,
                               fcs_int shifted);

    /* collect grid from neighbor processes */
    static void gather_grid(data_struct* d, fcs_float* rs_grid);
    static void add_block(fcs_float *in, fcs_float *out, int start[3], int size[3], int dim[3]);
    /* spread grid to neighbor processors */
    static void spread_grid(data_struct* d, fcs_float* rs_grid);

    /* apply energy optimized influence function */
    static void apply_energy_influence_function(data_struct* d);
    /* apply force optimized influence function */
    static void apply_force_influence_function(data_struct* d);
    /* differentiate kspace in direction dim */
    static void ik_diff(data_struct* d, int dim);

    /* compute the total energy (in k-space, so no backtransform) */
    static fcs_float compute_total_energy(data_struct* d);
    /* assign the potentials to the positions */
    static void 
    assign_potentials(data_struct* d, fcs_float *data, 
                      fcs_int num_particles, 
                      fcs_float* positions, fcs_float* charges,
                      fcs_int shifted,
                      fcs_float* potentials);
    
    static void collect_print_timings(data_struct *d);

#ifdef P3M_IK
    /* assign the fields to the positions in dimension dim [IK]*/
    static void 
    assign_fields_ik(data_struct* d, 
                     fcs_float *data,
                     fcs_int dim,
                     fcs_int num_particles, 
                     fcs_float* positions,
                     fcs_int shifted,
                     fcs_float* fields);
#endif
#ifdef P3M_AD
    /* Backinterpolate the forces obtained from k-space to the positions [AD]*/
    static void 
    assign_fields_ad(data_struct* d,
                     fcs_float *data,
                     fcs_int num_real_particles, 
                     fcs_float* positions,
                     fcs_int shifted,
                     fcs_float* fields);
#endif


    /***************************************************/
    /* IMPLEMENTATION */
    /***************************************************/
#define START(ID)                                               \
    if (d->require_timings!=NONE) d->timings[ID] += -MPI_Wtime();
#define STOP(ID)                                                \
    if (d->require_timings!=NONE) d->timings[ID] += MPI_Wtime();
#define STOPSTART(ID1, ID2)                     \
    if (d->require_timings!=NONE) {                   \
      d->timings[ID1] += MPI_Wtime();           \
      d->timings[ID2] += -MPI_Wtime();          \
    }                                           \


#if defined(P3M_INTERLACE) && defined(P3M_AD)
    void run(data_struct* d,
             fcs_int _num_particles,
             fcs_float *_positions, 
             fcs_float *_charges,
             fcs_float *_fields,
             fcs_float *_potentials) {
      /* Here we assume, that the method is already tuned and that all
         parameters are valid */
      P3M_INFO(printf( "P3M::run() started...\n"));
      P3M_DEBUG(printf("  type of timing: %d\n", d->require_timings));
      /* reset all timers */
      if (d->require_timings!=NONE) {
        for (int i = 0; i < NUM_TIMINGS; i++)
          d->timings[i] = 0.0;
      }

      P3M_INFO(printf("    system parameters: box_l=" F3FLOAT "\n", \
                      d->box_l[0], d->box_l[1], d->box_l[2]));
      P3M_INFO(printf(                                                      \
                      "    p3m params: "                                    \
                      "r_cut=" FFLOAT ", grid=" F3INT ", cao=" FINT ", "    \
                      "alpha=" FFLOAT ", grid_off=" F3FLOAT "\n",           \
                      d->r_cut, d->grid[0], d->grid[1], d->grid[2], d->cao, \
                      d->alpha, d->grid_off[0], d->grid_off[1], d->grid_off[2]));

      P3M_DEBUG_LOCAL(MPI_Barrier(d->comm.mpicomm));
      P3M_DEBUG_LOCAL(printf("    %d: num_particles=%d\n",	\
                             d->comm.rank, _num_particles));


      /* decompose system */
      fcs_int num_real_particles;
      fcs_int num_ghost_particles;
      fcs_float *positions, *ghost_positions;
      fcs_float *charges, *ghost_charges;
      fcs_gridsort_index_t *indices, *ghost_indices;
      fcs_gridsort_t gridsort;

      START(TIMING_DECOMP)
        domain_decompose(d, &gridsort, 
                         _num_particles, _positions, _charges,
                         &num_real_particles,
                         &positions, &charges, &indices,
                         &num_ghost_particles,
                         &ghost_positions, &ghost_charges, &ghost_indices);

      /* allocate local fields and potentials */
      fcs_float *fields = NULL; 
      fcs_float *potentials = NULL; 
      if (_fields != NULL)
        fields = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*3*num_real_particles));
      if (_potentials != NULL || d->require_total_energy)
        potentials = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*num_real_particles));
  
      STOPSTART(TIMING_DECOMP, TIMING_CA);

      /* charge assignment */
      assign_charges(d, d->fft.data_buf, num_real_particles, 
                     positions, charges, 0);
      STOPSTART(TIMING_CA, TIMING_GATHER);
      /* gather the ca grid */
      gather_grid(d, d->fft.data_buf);

      // Complexify
      for (fcs_int i = d->local_grid.size-1; i >= 0; i--)
        d->rs_grid[2*i] = d->fft.data_buf[i];

      STOPSTART(TIMING_GATHER, TIMING_CA);

      // Second (shifted) run
      /* charge assignment */
      assign_charges(d, d->fft.data_buf, num_real_particles, 
                     positions, charges, 1);

      STOPSTART(TIMING_CA, TIMING_GATHER);

      /* gather the ca grid */
      gather_grid(d, d->fft.data_buf);
      /* now d->rs_grid should contain the local ca grid */

      // Complexify
      for (fcs_int i = d->local_grid.size-1; i >= 0; i--)
        d->rs_grid[2*i+1] = d->fft.data_buf[i];
      STOP(TIMING_GATHER);

      /* forward transform */
    P3M_DEBUG(printf("  calling fft_perform_forw()...\n"));
    START(TIMING_FORWARD);
    fft_perform_forw(&d->fft, &d->comm, d->rs_grid);
    STOP(TIMING_FORWARD);
    P3M_DEBUG(printf("  returned from fft_perform_forw().\n"));
    
      /********************************************/
      /* POTENTIAL COMPUTATION */
      /********************************************/
      if (d->require_total_energy || _potentials != NULL) {
        /* apply energy optimized influence function */
        START(TIMING_INFLUENCE)
          apply_energy_influence_function(d);
        /* result is in d->ks_grid */
        STOP(TIMING_INFLUENCE)

          /* compute total energy, but not potentials */
          if (d->require_total_energy && potentials == NULL) {
            START(TIMING_POTENTIALS)
              d->total_energy = compute_total_energy(d);
            STOP(TIMING_POTENTIALS)
              }

        if (_potentials != NULL) {
          /* backtransform the grid */
          
          if( d->require_timings != ESTIMATE_ALL
           && d->require_timings != ESTIMATE_FFT ){
              P3M_DEBUG(printf( "  calling fft_perform_back (potentials)...\n"));
              START(TIMING_BACK)
              fft_perform_back(&d->fft, &d->comm, d->ks_grid);
              STOP(TIMING_BACK)
              P3M_DEBUG(printf( "  returned from fft_perform_back.\n"));
          }
          /** First (unshifted) run */
          START(TIMING_SPREAD)
            for (fcs_int i=0; i<d->local_grid.size; i++) {
              d->fft.data_buf[i] = d->ks_grid[2*i];
            } 

          spread_grid(d, d->fft.data_buf);

          STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS)

            assign_potentials(d, d->fft.data_buf,
                              num_real_particles, positions, 
                              charges, 0, potentials);

          STOPSTART(TIMING_POTENTIALS, TIMING_SPREAD)

            /** Second (shifted) run */
            for (fcs_int i=0; i<d->local_grid.size; i++) {
              d->fft.data_buf[i] = d->ks_grid[2*i+1];
            }
          spread_grid(d, d->fft.data_buf);

          STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS)

            assign_potentials(d, d->fft.data_buf,
                              num_real_particles, positions,
                              charges, 1, potentials);

          STOP(TIMING_POTENTIALS)
            }
      }

      /********************************************/
      /* FIELD COMPUTATION */
      /********************************************/
      if (_fields != NULL) {
        /* apply force optimized influence function */
        START(TIMING_INFLUENCE);
        apply_force_influence_function(d);
        STOP(TIMING_INFLUENCE);
    
        /* backtransform the grid */
        if(d->require_timings != ESTIMATE_ALL
        && d->require_timings != ESTIMATE_FFT){
            P3M_DEBUG(printf("  calling fft_perform_back...\n"));
            START(TIMING_BACK);
            fft_perform_back(&d->fft, &d->comm, d->ks_grid);
            STOP(TIMING_BACK);
            P3M_DEBUG(printf("  returned from fft_perform_back.\n"));            
        }
        
        START(TIMING_SPREAD);
          /* First (unshifted) run */
          P3M_INFO(printf("  computing unshifted grid\n"));
        for (fcs_int i=0; i<d->local_grid.size; i++) {
          d->fft.data_buf[i] = d->ks_grid[2*i];
        } 
    
        spread_grid(d, d->fft.data_buf);

        STOPSTART(TIMING_SPREAD, TIMING_FIELDS)

          assign_fields_ad(d, d->fft.data_buf, num_real_particles, 
                           positions, 0, fields);

        STOPSTART(TIMING_FIELDS, TIMING_SPREAD)
    
          /* Second (shifted) run */
          P3M_INFO(printf("  computing shifted grid\n"));
        for (fcs_int i=0; i<d->local_grid.size; i++) {
          d->fft.data_buf[i] = d->ks_grid[2*i+1];
        }
    
        spread_grid(d, d->fft.data_buf);
    
        STOPSTART(TIMING_SPREAD, TIMING_FIELDS)

          assign_fields_ad(d, d->fft.data_buf, num_real_particles, 
                           positions, 1, fields);

            STOP(TIMING_FIELDS)
        }

    if (d->near_field_flag) {
        /* start near timer */
        START(TIMING_NEAR)

          /* compute near field */
          fcs_near_t near;
        fcs_float alpha = d->alpha;
  
        fcs_near_create(&near);
        /*  fcs_near_set_field_potential(&near, compute_near);*/
        fcs_near_set_loop(&near, compute_near_loop);

        fcs_float box_base[3] = {0.0, 0.0, 0.0 };
        fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
        fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
        fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
        fcs_near_set_system(&near, box_base, box_a, box_b, box_c, NULL);

        fcs_near_set_particles(&near, num_real_particles, num_real_particles,
                               positions, charges, indices,
                               (_fields != NULL) ? fields : NULL, 
                               (_potentials != NULL) ? potentials : NULL);

        fcs_near_set_ghosts(&near, num_ghost_particles,
                            ghost_positions, ghost_charges, ghost_indices);

        P3M_DEBUG(printf( "  calling fcs_near_compute()...\n"));
        fcs_near_compute(&near, d->r_cut, &alpha, d->comm.mpicomm);
        P3M_DEBUG(printf( "  returning from fcs_near_compute().\n"));
 
        fcs_near_destroy(&near);

        STOP(TIMING_NEAR)
          }

      START(TIMING_COMP)
        /* sort particles back */
        P3M_DEBUG(printf( "  calling fcs_gridsort_sort_backward()...\n"));
      fcs_gridsort_sort_backward(&gridsort,
                                 fields, potentials,
                                 _fields, _potentials, 1,
                                 d->comm.mpicomm);
      P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_backward().\n"));
  
      fcs_gridsort_free(&gridsort);
      fcs_gridsort_destroy(&gridsort);
      
      if (fields != NULL) free(fields);
      if (potentials != NULL) free(potentials);

      STOP(TIMING_COMP)
              
      /* estimate FFT back timing*/
      if(d->require_timings == ESTIMATE_ALL
      || d->require_timings == ESTIMATE_FFT) {
          if (_potentials != NULL)
              d->timings[TIMING_BACK] = 2 * d->timings[TIMING_FORWARD];
          else
              d->timings[TIMING_BACK] = d->timings[TIMING_FORWARD];
      }  
      
        /* collect timings from the different nodes */
        if (d->require_timings!=NONE) collect_print_timings(d);

      P3M_INFO(printf( "P3M::run() finished.\n"));
    }

#elif !defined(P3M_INTERLACE) && defined(P3M_IK)
    void run(data_struct* rd,
             fcs_int _num_particles,
             fcs_float *_positions, 
             fcs_float *_charges,
             fcs_float *_fields,
             fcs_float *_potentials) {
      /* Here we assume, that the method is already tuned and that all
         parameters are valid */
      P3M_INFO(printf( "ifcs_p3m_run() started...\n"));
      data_struct *d = (data_struct*)rd;
      
      P3M_DEBUG(printf("type of timing: %d\n",d->require_timings));

      /* reset all timers */
      if (d->require_timings!=NONE) {
        for (int i = 0; i < NUM_TIMINGS; i++)
          d->timings[i] = 0.0;
      }

      P3M_INFO(printf("    system parameters: box_l=" F3FLOAT "\n", \
                      d->box_l[0], d->box_l[1], d->box_l[2]));
      P3M_INFO(printf(                                                      \
                      "    p3m params: "                                    \
                      "r_cut=" FFLOAT ", grid=" F3INT ", cao=" FINT ", "    \
                      "alpha=" FFLOAT ", grid_off=" F3FLOAT "\n",           \
                      d->r_cut, d->grid[0], d->grid[1], d->grid[2], d->cao, \
                      d->alpha, d->grid_off[0], d->grid_off[1], d->grid_off[2]));

      P3M_DEBUG_LOCAL(MPI_Barrier(d->comm.mpicomm));
      P3M_DEBUG_LOCAL(printf("    %d: num_particles=%d\n",	\
                             d->comm.rank, _num_particles));


      /* decompose system */
      fcs_int num_real_particles;
      fcs_int num_ghost_particles;
      fcs_float *positions, *ghost_positions;
      fcs_float *charges, *ghost_charges;
      fcs_gridsort_index_t *indices, *ghost_indices;
      fcs_gridsort_t gridsort;

      START(TIMING_DECOMP)
        domain_decompose(d, &gridsort, 
                         _num_particles, _positions, _charges,
                         &num_real_particles,
                         &positions, &charges, &indices,
                         &num_ghost_particles,
                         &ghost_positions, &ghost_charges, &ghost_indices);

      /* allocate local fields and potentials */
      fcs_float *fields = NULL; 
      fcs_float *potentials = NULL; 
      if (_fields != NULL)
        fields = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*3*num_real_particles));
      if (_potentials != NULL || d->require_total_energy)
        potentials = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*num_real_particles));
  
      STOPSTART(TIMING_DECOMP, TIMING_CA);

      /* charge assignment */
      assign_charges(d, d->rs_grid, num_real_particles, 
                     positions, charges, 0);
      STOPSTART(TIMING_CA, TIMING_GATHER);
      /* gather the ca grid */
      gather_grid(d, d->rs_grid);
      /* now d->rs_grid should contain the local ca grid */
      STOP(TIMING_GATHER);

      /* forward transform */
      
      P3M_DEBUG(printf( "  calling fft_perform_forw()...\n"));
      START(TIMING_FORWARD);
      fft_perform_forw(&d->fft, &d->comm, d->rs_grid);
      STOP(TIMING_FORWARD);
      P3M_DEBUG(printf("  returned from fft_perform_forw().\n"));
      
      /********************************************/
      /* POTENTIAL COMPUTATION */
      /********************************************/
      if (d->require_total_energy || _potentials != NULL) {
        /* apply energy optimized influence function */
        START(TIMING_INFLUENCE)
          apply_energy_influence_function(d);
        /* result is in d->ks_grid */
        STOP(TIMING_INFLUENCE)

          /* compute total energy, but not potentials */
          if (d->require_total_energy && potentials == NULL) {
            START(TIMING_POTENTIALS)
              d->total_energy = compute_total_energy(d);
            STOP(TIMING_POTENTIALS)
              }

        if (_potentials != NULL) {
          /* backtransform the grid */
          
          if(d->require_timings != ESTIMATE_ALL
          && d->require_timings != ESTIMATE_FFT){
              P3M_DEBUG(printf( "  calling fft_perform_back (potentials)...\n"));
              START(TIMING_BACK)
              fft_perform_back(&d->fft, &d->comm, d->ks_grid);
              STOP(TIMING_BACK)
              P3M_DEBUG(printf( "  returned from fft_perform_back.\n"));
          }
          /* redistribute energy grid */
          START(TIMING_SPREAD);
          spread_grid(d, d->ks_grid);

          STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS);

          /* compute potentials */
          assign_potentials(d, d->ks_grid,
                            num_real_particles, positions,
                            charges, 0,
                            potentials);

          STOP(TIMING_POTENTIALS);
        }
      }

      /********************************************/
      /* FIELD COMPUTATION */
      /********************************************/
      if (_fields != NULL) {
        /* apply force optimized influence function */
        START(TIMING_INFLUENCE);
        apply_force_influence_function(d);
        STOP(TIMING_INFLUENCE);
    
        /* result is in d->ks_grid */
        for (int dim = 0; dim < 3; dim++) {
          /* differentiate in direction dim */
          /* result is stored in d->rs_grid */
          START(TIMING_FIELDS);
          ik_diff(d, dim);     
          STOP(TIMING_FIELDS);
          
        if(d->require_timings != ESTIMATE_ALL
        && d->require_timings != ESTIMATE_FFT){
            
           /* backtransform the grid */
            P3M_DEBUG(printf( "  calling fft_perform_back (field dim=%d)...\n", dim));
            START(TIMING_BACK);
            fft_perform_back(&d->fft, &d->comm, d->rs_grid);
            STOP(TIMING_BACK);
            P3M_DEBUG(printf("  returned from fft_perform_back.\n"));            
        }
          
          /* redistribute force grid */
        START(TIMING_SPREAD);          
          spread_grid(d, d->rs_grid);
          STOPSTART(TIMING_SPREAD, TIMING_FIELDS);
          assign_fields_ik(d, d->rs_grid, dim, num_real_particles, 
                           positions, 0, fields);
          STOP(TIMING_FIELDS);
          }
        }

      if (d->near_field_flag) {
        /* start near timer */
        START(TIMING_NEAR)

          /* compute near field */
          fcs_near_t near;
        fcs_float alpha = d->alpha;
  
        fcs_near_create(&near);
        /*  fcs_near_set_field_potential(&near, compute_near);*/
        fcs_near_set_loop(&near, compute_near_loop);

        fcs_float box_base[3] = {0.0, 0.0, 0.0 };
        fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
        fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
        fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
        fcs_near_set_system(&near, box_base, box_a, box_b, box_c, NULL);

        fcs_near_set_particles(&near, num_real_particles, num_real_particles,
                               positions, charges, indices,
                               (_fields != NULL) ? fields : NULL, 
                               (_potentials != NULL) ? potentials : NULL);

        fcs_near_set_ghosts(&near, num_ghost_particles,
                            ghost_positions, ghost_charges, ghost_indices);

        P3M_DEBUG(printf( "  calling fcs_near_compute()...\n"));
        fcs_near_compute(&near, d->r_cut, &alpha, d->comm.mpicomm);
        P3M_DEBUG(printf( "  returning from fcs_near_compute().\n"));
 
        fcs_near_destroy(&near);

        STOP(TIMING_NEAR)
          }

      START(TIMING_COMP)
        /* sort particles back */
        P3M_DEBUG(printf( "  calling fcs_gridsort_sort_backward()...\n"));
      fcs_gridsort_sort_backward(&gridsort,
                                 fields, potentials,
                                 _fields, _potentials, 1,
                                 d->comm.mpicomm);
      P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_backward().\n"));
  
      fcs_gridsort_free(&gridsort);
      fcs_gridsort_destroy(&gridsort);

      if (fields != NULL) free(fields);
      if (potentials != NULL) free(potentials);

      STOP(TIMING_COMP)
              
              /* estimate FFT back timing*/
if (d->require_timings == ESTIMATE_ALL || d->require_timings == ESTIMATE_FFT) {
            if (_potentials != NULL)
                d->timings[TIMING_BACK] = 4 * d->timings[TIMING_FORWARD];
            else
                d->timings[TIMING_BACK] = 3 * d->timings[TIMING_FORWARD];
        
            /* collect timings from the different nodes */
        if (d->require_timings!=NONE) collect_print_timings(d);       

      P3M_INFO(printf( "ifcs_p3m_run() finished.\n"));
            }
        }
#endif

    /***************************************************/
    /* RUN COMPONENTS */
    static void
    domain_decompose(data_struct *d, fcs_gridsort_t *gridsort,
                     fcs_int _num_particles,
                     fcs_float *_positions, fcs_float *_charges,
                     fcs_int *num_real_particles,
                     fcs_float **positions, fcs_float **charges,
                     fcs_gridsort_index_t **indices, 
                     fcs_int *num_ghost_particles,
                     fcs_float **ghost_positions, fcs_float **ghost_charges,
                     fcs_gridsort_index_t **ghost_indices
                     ) {
      fcs_float box_base[3] = {0.0, 0.0, 0.0 };
      fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
      fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
      fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
      fcs_int num_particles;

      fcs_gridsort_create(gridsort);
  
      fcs_gridsort_set_system(gridsort, box_base, box_a, box_b, box_c, NULL);
      fcs_gridsort_set_particles(gridsort, _num_particles, _num_particles, 
                                 _positions, _charges);

      P3M_DEBUG(printf( "  calling fcs_gridsort_sort_forward()...\n"));
      /* @todo: Set skin to r_cut only, when near field is wanted! */
      fcs_gridsort_sort_forward(gridsort,
                                (d->near_field_flag ? d->r_cut : 0.0),
                                d->comm.mpicomm);
      P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_forward().\n"));
      fcs_gridsort_separate_ghosts(gridsort, 
                                   num_real_particles, 
                                   num_ghost_particles);

      fcs_gridsort_get_sorted_particles(gridsort, 
                                        &num_particles, NULL, NULL, NULL, NULL);

      fcs_gridsort_get_real_particles(gridsort, 
                                      num_real_particles, 
                                      positions, charges, 
                                      indices);
  
      fcs_gridsort_get_ghost_particles(gridsort, num_ghost_particles, 
                                       ghost_positions, ghost_charges, 
                                       ghost_indices);

      P3M_DEBUG_LOCAL(MPI_Barrier(d->comm.mpicomm));
      P3M_DEBUG_LOCAL(printf(                                               \
                             "    %d: num_particles=%d"                     \
                             " num_real_particles=%d"                       \
                             " num_ghost_particles=%d\n",                   \
                             d->comm.rank, num_particles,                   \
                             *num_real_particles, *num_ghost_particles));
    }


    /***************************************************/
    /* CHARGE ASSIGNMENT */
    /** Compute the data of the charge assignment grid points.

        The function returns the linear index of the top left grid point
        in the charge assignment grid that corresponds to real_pos. When
        "shifted" is set, it uses the shifted position for interlacing.
        After the call, caf_cache contains a cache of the values of the
        charge assignment fraction (caf) for x,y,z.
    */
    static fcs_int 
    get_ca_points(data_struct *d, 
                  fcs_float real_pos[3], 
                  fcs_int shifted) {
      /* linear index of the grid point */
      fcs_int linind = 0;

      for (fcs_int dim=0; dim<3; dim++) {
        /* position in normalized coordinates in [0,1] */
        fcs_float pos = (real_pos[dim] - d->local_grid.ld_pos[dim]) * d->ai[dim];
        /* shift position to the corner of the charge assignment area */
        pos -= d->pos_shift;
        /* if using the interlaced grid, shift it more */
        if (shifted) pos -= 0.5;
        /* nearest grid point in the ca grid */
        fcs_int grid_ind  = (fcs_int)floor(pos);
#ifdef ADDITIONAL_CHECKS
        if (grid_ind < 0) {
          printf("grid_ind[%d]=%d < 0\n", dim, grid_ind);
          printf("pos=%lf, pos_shift=%lf, real_pos=%lf, ld_pos=%lf, ai=%lf, shifted=%d", 
                 pos, d->pos_shift, real_pos[dim], d->local_grid.ld_pos[dim],  d->ai[dim], 
                 shifted);
        } else if (grid_ind > d->local_grid.dim[dim]) {
          printf("grid_ind[%d]=%d > %d\n", dim, grid_ind, d->local_grid.dim[dim]);
          printf("pos=%lf, pos_shift=%lf, real_pos=%lf, ld_pos=%lf, ai=%lf, shifted=%d", 
                 pos, d->pos_shift, real_pos[dim], d->local_grid.ld_pos[dim],  d->ai[dim], 
                 shifted);
        }
#endif
        /* linear index of grid point */
        linind = grid_ind + linind*d->local_grid.dim[dim];
        /* normalized distance to grid point */
        fcs_float dist = (pos-grid_ind)-0.5;

        switch (dim) {
        case 0: d->cafx->update(dist); break;
        case 1: d->cafy->update(dist); break;
        case 2: d->cafz->update(dist); break;
        }
      
#ifdef P3M_AD
        switch (dim) {
        case 0: d->cafx_d->update(dist); break;
        case 1: d->cafy_d->update(dist); break;
        case 2: d->cafz_d->update(dist); break;
        }
#endif

#ifdef ADDITIONAL_CHECKS
        if (real_pos[dim] < d->comm.my_left[dim] 
            || real_pos[dim] > d->comm.my_right[dim]) {
          printf("%d: dim %d: position not in domain! " F3FLOAT "\n", 
                 d->comm.rank, dim, real_pos[0], real_pos[1], real_pos[2]);
        }
#endif
      }

#ifdef ADDITIONAL_CHECKS
      if (linind < 0) {
        printf("ERROR: %d: linind %d < 0\n", d->comm.rank, linind);
      } else if (linind >= d->local_grid.size) {
        printf("ERROR: %d: linind %d > %d\n", d->comm.rank, linind, d->local_grid.size);
      }
#endif

      return linind;
    }

    /** Assign the charges to the grid */
    static void 
    assign_charges(data_struct* d,
                   fcs_float *data,
                   fcs_int num_real_particles,
                   fcs_float *positions, 
                   fcs_float *charges,
                   fcs_int shifted) {
      P3M_DEBUG(printf( "  assign_charges() started...\n"));

      const fcs_int q2off = d->local_grid.q_2_off;
      const fcs_int q21off = d->local_grid.q_21_off;

      /* init local charge grid */
      for (fcs_int i=0; i<d->local_grid.size; i++) data[i] = 0.0;

      /* now assign the charges */
      for (fcs_int pid=0; pid < num_real_particles; pid++) {
        const fcs_float q = charges[pid];
        fcs_int linind_grid = get_ca_points(d, &positions[pid*3], shifted);

        /* Loop over all ca grid points nearby and compute charge assignment fraction */
        for (fcs_float *caf_x = d->cafx->begin(); caf_x < d->cafx->end(); caf_x++) {
          for (fcs_float *caf_y = d->cafy->begin(); caf_y < d->cafy->end(); caf_y++) {
            fcs_float caf_xy = *caf_x * *caf_y;
            for (fcs_float *caf_z = d->cafz->begin(); caf_z < d->cafz->end(); caf_z++) {
              /* add it to the grid */
              data[linind_grid] += q * caf_xy * *caf_z;
              linind_grid++;
            }
            linind_grid += q2off;
          }
          linind_grid += q21off;
        }
      }

      P3M_DEBUG(printf( "  assign_charges() finished...\n"));
    }

    /* Gather information for FFT grid inside the nodes domain (inner local grid) */
    static void gather_grid(data_struct* d, fcs_float* rs_grid) {
      MPI_Status status;
      fcs_float *tmp_ptr;

      P3M_DEBUG(printf( "  gather_grid() started...\n"));

      /* direction loop */
      for(fcs_int s_dir=0; s_dir<6; s_dir++) {
        fcs_int r_dir;
        if(s_dir%2==0) r_dir = s_dir+1;
        else           r_dir = s_dir-1;
        /* pack send block */
        if (d->sm.s_size[s_dir]>0)
          fft_pack_block(rs_grid, d->send_grid, d->sm.s_ld[s_dir], 
                         d->sm.s_dim[s_dir], d->local_grid.dim, 1);
      
        /* communication */
        /** @todo Replace with MPI_Sendrecv */
        if (d->comm.node_neighbors[s_dir] != d->comm.rank) {
          for (fcs_int evenodd=0; evenodd<2; evenodd++) {
            if ((d->comm.node_pos[s_dir/2]+evenodd)%2 == 0) {
              if (d->sm.s_size[s_dir] > 0) {
                P3M_DEBUG_LOCAL(printf("    %d: sending %d floats to %d (s_dir=%d)\n", \
                                       d->comm.rank, d->sm.s_size[s_dir],	\
                                       d->comm.node_neighbors[s_dir], s_dir));
                MPI_Send(d->send_grid, d->sm.s_size[s_dir], FCS_MPI_FLOAT,
                         d->comm.node_neighbors[s_dir], REQ_P3M_GATHER, d->comm.mpicomm);
              }
            } else {
              if (d->sm.r_size[r_dir] > 0) {
                P3M_DEBUG_LOCAL(printf( "    %d: receiving %d floats from %d (r_dir=%d)\n", \
                                        d->comm.rank, d->sm.r_size[r_dir],	\
                                        d->comm.node_neighbors[r_dir], r_dir));
                MPI_Recv(d->recv_grid, d->sm.r_size[r_dir], FCS_MPI_FLOAT,
                         d->comm.node_neighbors[r_dir], REQ_P3M_GATHER, d->comm.mpicomm, &status);
              }
            }
          }
        } else {
          tmp_ptr = d->recv_grid;
          d->recv_grid = d->send_grid;
          d->send_grid = tmp_ptr;
        }
        /* add recv block */
        if(d->sm.r_size[r_dir]>0) {
          add_block(d->recv_grid, rs_grid, d->sm.r_ld[r_dir], 
                    d->sm.r_dim[r_dir], d->local_grid.dim);
        }
      }

      P3M_DEBUG(printf( "  gather_grid() finished.\n"));
    }

    static void add_block(fcs_float *in, fcs_float *out, int start[3], int size[3], int dim[3]) {
      /* fast,mid and slow changing indices */
      int f,m,s;
      /* linear index of in grid, linear index of out grid */
      int li_in=0,li_out=0;
      /* offsets for indizes in output grid */
      int m_out_offset,s_out_offset;

      li_out = start[2] + ( dim[2]*( start[1] + (dim[1]*start[0]) ) );
      m_out_offset  = dim[2] - size[2];
      s_out_offset  = (dim[2] * (dim[1] - size[1]));

      for(s=0 ;s<size[0]; s++) {
        for(m=0; m<size[1]; m++) {
          for(f=0; f<size[2]; f++) {
            out[li_out++] += in[li_in++];
          }
          li_out += m_out_offset;
        }
        li_out += s_out_offset;
      }
    }


    /* apply the influence function */
    static void apply_energy_influence_function(data_struct* d) {
      P3M_DEBUG(printf( "  apply_energy_influence_function() started...\n"));
      const fcs_int size = d->fft.plan[3].new_size;
      for (fcs_int i=0; i < size; i++) {
        d->ks_grid[2*i] = d->g_energy[i] * d->rs_grid[2*i]; 
        d->ks_grid[2*i+1] = d->g_energy[i] * d->rs_grid[2*i+1]; 
      }
      P3M_DEBUG(printf( "  apply_energy_influence_function() finished.\n"));
    }

    /* apply the influence function */
    static void apply_force_influence_function(data_struct* d) {
      P3M_DEBUG(printf( "  apply_force_influence_function() started...\n"));
      const fcs_int size = d->fft.plan[3].new_size;
      for (fcs_int i=0; i < size; i++) {
        d->ks_grid[2*i] = d->g_force[i] * d->rs_grid[2*i]; 
        d->ks_grid[2*i+1] = d->g_force[i] * d->rs_grid[2*i+1]; 
      }
      P3M_DEBUG(printf( "  apply_force_influence_function() finished.\n"));
    }
    
    /* Add up the measured timings and collect information from all nodes.
     * Print the timings if requested so.
     */
        static void collect_print_timings(data_struct *d) {
            d->timings[TIMING_FAR] += d->timings[TIMING_CA];
            d->timings[TIMING_FAR] += d->timings[TIMING_GATHER];
            d->timings[TIMING_FAR] += d->timings[TIMING_FORWARD];
            d->timings[TIMING_FAR] += d->timings[TIMING_BACK];
            d->timings[TIMING_FAR] += d->timings[TIMING_INFLUENCE];
            d->timings[TIMING_FAR] += d->timings[TIMING_SPREAD];
            d->timings[TIMING_FAR] += d->timings[TIMING_POTENTIALS];
            d->timings[TIMING_FAR] += d->timings[TIMING_FIELDS];

            d->timings[TIMING] += d->timings[TIMING_DECOMP];
            d->timings[TIMING] += d->timings[TIMING_FAR];
            d->timings[TIMING] += d->timings[TIMING_NEAR];
            d->timings[TIMING] += d->timings[TIMING_COMP];

            if (on_root())
                MPI_Reduce(MPI_IN_PLACE, d->timings,
                    NUM_TIMINGS, MPI_DOUBLE, MPI_MAX,
                    0, d->comm.mpicomm);
            else
                MPI_Reduce(d->timings, 0,
                    NUM_TIMINGS, MPI_DOUBLE, MPI_MAX,
                    0, d->comm.mpicomm);

#ifdef P3M_PRINT_TIMINGS
          printf("  P3M TIMINGS:\n");
          printf("    total=%le (%lf)\n", d->timings[TIMING], 1.0);
          printf("      far=%le (%lf)\n", d->timings[TIMING_FAR], 
                 d->timings[TIMING_FAR]/d->timings[TIMING]);
          printf("     near=%le (%lf)\n", d->timings[TIMING_NEAR], 
                 d->timings[TIMING_NEAR]/d->timings[TIMING]);
          printf("       ca=%le (%lf)\n", d->timings[TIMING_CA], 
                 d->timings[TIMING_CA]/d->timings[TIMING]);
          printf("      pot=%le (%lf)", d->timings[TIMING_POTENTIALS], 
                 d->timings[TIMING_POTENTIALS]/d->timings[TIMING]);
          //if(d->require_timings == ESTIMATE_ALL //not yet implemented
          //|| d->require_timings == ESTIMATE_ASSIGNMENT)
          //    printf(" (empirical estimate)");
          printf("\n");
          printf("   fields=%le (%lf)", d->timings[TIMING_FIELDS], 
                 d->timings[TIMING_FIELDS]/d->timings[TIMING]);
          //if(d->require_timings == ESTIMATE_ALL //not yet implemented
          //|| d->require_timings == ESTIMATE_ASSIGNMENT)
          //    printf(" (empirical estimate)");
          printf("\n");
          printf("   gather=%le (%lf)\n", d->timings[TIMING_GATHER], 
                 d->timings[TIMING_GATHER]/d->timings[TIMING]);
          printf("   spread=%le (%lf)\n", d->timings[TIMING_SPREAD], 
                 d->timings[TIMING_SPREAD]/d->timings[TIMING]);
          printf("  forward=%le (%lf)\n", d->timings[TIMING_FORWARD], 
                 d->timings[TIMING_FORWARD]/d->timings[TIMING]);
          printf("     back=%le (%lf)", d->timings[TIMING_BACK], 
                 d->timings[TIMING_BACK]/d->timings[TIMING]);
          if(d->require_timings == ESTIMATE_ALL
          || d->require_timings == ESTIMATE_FFT)
              printf(" (theoretical estimate)");
          printf("\n");
          printf("   decomp=%le (%lf)\n", d->timings[TIMING_DECOMP], 
                 d->timings[TIMING_DECOMP]/d->timings[TIMING]);
          printf("     comp=%le (%lf)\n", d->timings[TIMING_COMP], 
                 d->timings[TIMING_COMP]/d->timings[TIMING]);
#endif

        }

    static void ik_diff(data_struct* d, int dim) {
      fcs_int ind;
      fcs_int j[3];
      fcs_int* d_operator = NULL;
      /* direction in k space: */
      fcs_int dim_rs = (dim+d->ks_pnum)%3;

      P3M_DEBUG(printf( "  ik_diff() started...\n"));
      switch (dim) {
      case KX:
        d_operator = d->d_op[RX];
        break;
      case KY:
        d_operator = d->d_op[RY];
        break;
      case KZ:
        d_operator = d->d_op[RZ];
      }
    
      /* srqt(-1)*k differentiation */
      ind=0;
      for(j[0]=0; j[0]<d->fft.plan[3].new_grid[0]; j[0]++) {
        for(j[1]=0; j[1]<d->fft.plan[3].new_grid[1]; j[1]++) {
          for(j[2]=0; j[2]<d->fft.plan[3].new_grid[2]; j[2]++) {
            /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */
            d->rs_grid[ind] =
              -2.0*M_PI*(d->ks_grid[ind+1] * d_operator[ j[dim]+d->fft.plan[3].start[dim] ])
              / d->box_l[dim_rs];
            d->rs_grid[ind+1] =
              2.0*M_PI*d->ks_grid[ind] * d_operator[ j[dim]+d->fft.plan[3].start[dim] ]
              / d->box_l[dim_rs];
            ind+=2;
          }
        }
      }

      P3M_DEBUG(printf( "  ik_diff() finished.\n"));
      /* store the result in d->rs_grid */
    }

    static void spread_grid(data_struct* d, fcs_float* rs_grid) {
      int s_dir,r_dir,evenodd;
      MPI_Status status;
      fcs_float *tmp_ptr;
      P3M_DEBUG(printf( "  spread_grid() started...\n"));
  
      /* direction loop */
      for(s_dir=5; s_dir>=0; s_dir--) {
        if(s_dir%2==0) r_dir = s_dir+1;
        else           r_dir = s_dir-1;
        /* pack send block */ 
        if(d->sm.s_size[s_dir]>0) 
          fft_pack_block(rs_grid, d->send_grid, d->sm.r_ld[r_dir], d->sm.r_dim[r_dir], d->local_grid.dim, 1);
        /* communication */
        /** @todo Replace with MPI_Sendrecv */
        if (d->comm.node_neighbors[r_dir] != d->comm.rank) {
          for (evenodd=0; evenodd<2;evenodd++) {
            if ((d->comm.node_pos[r_dir/2]+evenodd)%2==0) {
              if (d->sm.r_size[r_dir]>0) 
                MPI_Send(d->send_grid, d->sm.r_size[r_dir], FCS_MPI_FLOAT, 
                         d->comm.node_neighbors[r_dir], REQ_P3M_SPREAD, d->comm.mpicomm);
            }
            else {
              if (d->sm.s_size[s_dir]>0) 
                MPI_Recv(d->recv_grid, d->sm.s_size[s_dir], FCS_MPI_FLOAT, 
                         d->comm.node_neighbors[s_dir], REQ_P3M_SPREAD, d->comm.mpicomm, &status); 	    
            }
          }
        }
        else {
          tmp_ptr = d->recv_grid;
          d->recv_grid = d->send_grid;
          d->send_grid = tmp_ptr;
        }
        /* unpack recv block */
        if(d->sm.s_size[s_dir]>0) {
          fft_unpack_block(d->recv_grid, rs_grid, d->sm.s_ld[s_dir], d->sm.s_dim[s_dir], d->local_grid.dim, 1); 
        }
      }

      P3M_DEBUG(printf( "  spread_grid() finished.\n"));
    }


    /** Compute the total energy of the system in kspace. No need to
        backtransform the FFT grid in this case! */
    static fcs_float compute_total_energy(data_struct* d) {
      fcs_float local_k_space_energy;
      fcs_float k_space_energy;

      P3M_DEBUG(printf( "  compute_total_energy() started...\n"));

      local_k_space_energy = 0.0;
      for (fcs_int i=0; i < d->fft.plan[3].new_size; i++)
        /* Use the energy optimized influence function */
        local_k_space_energy += d->g_energy[i] * ( SQR(d->rs_grid[2*i]) + SQR(d->rs_grid[2*i+1]) );

      MPI_Reduce(&local_k_space_energy, &k_space_energy, 1, FCS_MPI_FLOAT, 
                 MPI_SUM, 0, d->comm.mpicomm);
      fcs_float prefactor = 1.0 / (2.0 * d->box_l[0] * d->box_l[1] * d->box_l[2]);
      k_space_energy *= prefactor;

#ifdef P3M_INTERLACE
      /* In the case of interlacing we have calculated the sum of the
         shifted and unshifted charges, we have to take the average. */
      k_space_energy *= 0.5;
#endif

      /* self energy correction */
      k_space_energy -= d->sum_q2 * d->alpha * 0.5*FCS_2_SQRTPI;
      /* net charge correction */
      k_space_energy -= d->square_sum_q * M_PI * prefactor / SQR(d->alpha);

      P3M_DEBUG(printf( "  compute_total_energy() finished.\n"));
      return k_space_energy;
    }

    /* Backinterpolate the potentials obtained from k-space to the positions */
    static void 
    assign_potentials(data_struct* d, 
                      fcs_float *data,
                      fcs_int num_real_particles, 
                      fcs_float* positions, fcs_float* charges, 
                      fcs_int shifted,
                      fcs_float* potentials) {
      const fcs_int q2off = d->local_grid.q_2_off;
      const fcs_int q21off = d->local_grid.q_21_off;
      const fcs_float prefactor = 1.0 / (d->box_l[0] * d->box_l[1] * d->box_l[2]);
  
      P3M_DEBUG(printf( "  assign_potentials() started...\n"));
      /* Loop over all particles */
      for (fcs_int pid=0; pid < num_real_particles; pid++) {
        fcs_float potential = 0.0;
        fcs_int linind_grid = 
          get_ca_points(d, &positions[pid*3], shifted);

        /* Loop over all ca grid points nearby and compute charge assignment fraction */
        for (fcs_float *caf_x = d->cafx->begin(); caf_x < d->cafx->end(); caf_x++) {
          for (fcs_float *caf_y = d->cafy->begin(); caf_y < d->cafy->end(); caf_y++) {
            fcs_float caf_xy = *caf_x * *caf_y;
            for (fcs_float *caf_z = d->cafz->begin(); caf_z < d->cafz->end(); caf_z++) {
              potential += *caf_z * caf_xy * data[linind_grid];
              linind_grid++;
            }
            linind_grid += q2off;
          }
          linind_grid += q21off;
        }

        potential *= prefactor;
        /* self energy correction */
        potential -= charges[pid] * FCS_2_SQRTPI * d->alpha;
        /* net charge correction */
        /* potential -= fabs(charges[pid]) * PI * prefactor / SQR(d->alpha); */

        /* store the result */
        if (!shifted) {
          potentials[pid] = potential;
        } else {
          potentials[pid] = 0.5*(potentials[pid] + potential);
        }

      }
      P3M_DEBUG(printf( "  assign_potentials() finished.\n"));
    }

    /* Backinterpolate the forces obtained from k-space to the positions */
    static void 
    assign_fields_ik(data_struct* d, 
                     fcs_float *data,
                     fcs_int dim,
                     fcs_int num_real_particles,
                     fcs_float* positions,
                     fcs_int shifted,
                     fcs_float* fields) {
      const fcs_int q2off = d->local_grid.q_2_off;
      const fcs_int q21off = d->local_grid.q_21_off;
      const fcs_float prefactor = 1.0 / (2.0 * d->box_l[0] * d->box_l[1] * d->box_l[2]);
      const fcs_int dim_rs = (dim+d->ks_pnum) % 3;

      P3M_DEBUG(printf( "  assign_fields() started...\n"));
      /* Loop over all particles */
      for (fcs_int pid=0; pid < num_real_particles; pid++) {
        fcs_float field = 0.0;
        fcs_int linind_grid = 
          get_ca_points(d, &positions[3*pid], shifted);

        /* loop over the local grid, compute the field */
        for (fcs_float *caf_x = d->cafx->begin(); caf_x < d->cafx->end(); caf_x++) {
          for (fcs_float *caf_y = d->cafy->begin(); caf_y < d->cafy->end(); caf_y++) {
            fcs_float caf_xy = *caf_x * *caf_y;
            for (fcs_float *caf_z = d->cafz->begin(); caf_z < d->cafz->end(); caf_z++) {
              field -= *caf_z * caf_xy * data[linind_grid];
              linind_grid++;
            }
            linind_grid += q2off;
          }
          linind_grid += q21off;
        }

        field *= prefactor;

        /* store the result */
        if (!shifted)
          fields[3*pid + dim_rs] = field;
        else
          fields[3*pid + dim_rs] = 0.5 * (fields[3*pid + dim_rs] + field);
      }
      P3M_DEBUG(printf( "  assign_fields() finished.\n"));
    }

#ifdef P3M_AD
    /* Backinterpolate the forces obtained from k-space to the positions */
    static void 
    assign_fields_ad(data_struct* d,
                     fcs_float *data,
                     fcs_int num_real_particles, 
                     fcs_float* positions,
                     fcs_int shifted,
                     fcs_float* fields) {
      const fcs_int q2off = d->local_grid.q_2_off;
      const fcs_int q21off = d->local_grid.q_21_off;
      const fcs_float prefactor = 1.0 / (d->box_l[0] * d->box_l[1] * d->box_l[2]);
      const fcs_float l_x_inv = 1.0/d->box_l[0];
      const fcs_float l_y_inv = 1.0/d->box_l[1];
      const fcs_float l_z_inv = 1.0/d->box_l[2];
      const fcs_float grid[3] = 
        { (fcs_float)d->grid[0], 
          (fcs_float)d->grid[1], 
          (fcs_float)d->grid[2] };


      P3M_DEBUG(printf( "  assign_fields() [AD] started...\n"));
      /* Loop over all particles */
      for (fcs_int pid = 0; pid < num_real_particles; pid++) {
        fcs_float field[3] = { 0.0, 0.0, 0.0 };
        fcs_int linind_grid = 
          get_ca_points(d, &positions[pid*3], shifted);

        fcs_float *caf_x_d = d->cafx_d->begin();
        for (fcs_float *caf_x = d->cafx->begin(); caf_x < d->cafx->end(); caf_x++) {
          fcs_float *caf_y_d = d->cafy_d->begin();
          for (fcs_float *caf_y = d->cafy->begin(); caf_y < d->cafy->end(); caf_y++) {
            fcs_float *caf_z_d = d->cafz_d->begin();
            for (fcs_float *caf_z = d->cafz->begin(); caf_z < d->cafz->end(); caf_z++) {
              field[0] -= *caf_x_d * *caf_y * *caf_z * l_x_inv 
                * data[linind_grid] * grid[0];
              field[1] -= *caf_x * *caf_y_d * *caf_z * l_y_inv 
                * data[linind_grid] * grid[1];
              field[2] -= *caf_x * *caf_y * *caf_z_d * l_z_inv 
                * data[linind_grid] * grid[2];
              linind_grid++;
              caf_z_d++;
            }
            linind_grid += q2off;
            caf_y_d++;
          }
          linind_grid += q21off;
          caf_x_d++;
        }
        field[0] *= prefactor;
        field[1] *= prefactor;
        field[2] *= prefactor;

        if (!shifted) {
          fields[3*pid + 0] = field[0]; 
          fields[3*pid + 1] = field[1]; 
          fields[3*pid + 2] = field[2]; 
        } else {
          fields[3*pid + 0] = 0.5*(fields[3*pid + 0] + field[0]);
          fields[3*pid + 1] = 0.5*(fields[3*pid + 1] + field[1]);
          fields[3*pid + 2] = 0.5*(fields[3*pid + 2] + field[2]);
        }
      }
      P3M_DEBUG(printf( "  assign_fields() finished.\n"));
    }
#endif
  }
}
