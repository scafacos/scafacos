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

#include "FCSCommon.h"
#include "run.h"
#include "types.h"
#include "utils.h"
#include "caf.h"
#include "tune.h"
#include "fcs_p3m_p.h"
#include "common/gridsort/gridsort.h"
#include "common/near/near.h"
#include <stdlib.h>
#include <stdio.h>

#if FCS_ENABLE_INFO 
#define FCS_P3M_INIT_TIMING(comm)               \
  int tm_rank;                                  \
  MPI_Comm_rank(comm, &tm_rank);                \
  double tm_timer, tm_global_timer;
#define FCS_P3M_START_TIMING()                  \
  tm_timer = -MPI_Wtime();
#define FCS_P3M_FINISH_TIMING(comm, str)                                \
  tm_timer += MPI_Wtime();                                              \
  MPI_Reduce(&tm_timer, &tm_global_timer, 1,                            \
             MPI_DOUBLE, MPI_MAX, 0, comm);                             \
  if(!tm_rank) printf("P3M_TIMING: %s takes %e s\n", str, tm_global_timer);
#else
#define FCS_P3M_INIT_TIMING(comm)
#define FCS_P3M_START_TIMING()
#define FCS_P3M_FINISH_TIMING(comm, str)
#endif

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
/* domain decomposition */
static void
ifcs_p3m_domain_decompose(ifcs_p3m_data_struct *d, fcs_gridsort_t *gridsort,
                          fcs_int _num_particles, fcs_int _max_num_particles, 
                          fcs_float *_positions, fcs_float *_charges,
                          fcs_int *num_real_particles,
                          fcs_float **positions, fcs_float **charges,
                          fcs_gridsort_index_t **indices, 
                          fcs_int *num_ghost_particles,
                          fcs_float **ghost_positions, fcs_float **ghost_charges,
                          fcs_gridsort_index_t **ghost_indices
                          );

/* charge assignment */
static void ifcs_p3m_assign_charges(ifcs_p3m_data_struct* d,
				    fcs_float *data,
				    fcs_int num_particles,
				    fcs_float *positions, 
				    fcs_float *charges,
				    fcs_int shifted);
static void ifcs_p3m_assign_single_charge(ifcs_p3m_data_struct *d,
					  fcs_float *data,
					  fcs_int charge_id,
					  fcs_float real_pos[3],
					  fcs_float q,
					  fcs_int shifted);

/* collect grid from neighbor processes */
static void ifcs_p3m_gather_grid(ifcs_p3m_data_struct* d, fcs_float* rs_grid);
static void ifcs_p3m_add_block(fcs_float *in, fcs_float *out, int start[3], int size[3], int dim[3]);
/* spread grid to neighbor processors */
static void ifcs_p3m_spread_grid(ifcs_p3m_data_struct* d, fcs_float* rs_grid);

/* apply energy optimized influence function */
static void ifcs_p3m_apply_energy_influence_function(ifcs_p3m_data_struct* d);
/* apply force optimized influence function */
static void ifcs_p3m_apply_force_influence_function(ifcs_p3m_data_struct* d);
/* differentiate kspace in direction dim */
static void ifcs_p3m_ik_diff(ifcs_p3m_data_struct* d, int dim);

/* compute the total energy (in k-space, so no backtransform) */
static fcs_float ifcs_p3m_compute_total_energy(ifcs_p3m_data_struct* d);
/* assign the potentials to the positions */
static void 
ifcs_p3m_assign_potentials(ifcs_p3m_data_struct* d, fcs_float *data, 
                           fcs_int num_particles, 
                           fcs_float* positions, fcs_float* charges,
                           fcs_int shifted,
                           fcs_float* potentials);

#ifdef P3M_IK
/* assign the fields to the positions in dimension dim [IK]*/
static void 
ifcs_p3m_assign_fields_ik(ifcs_p3m_data_struct* d, 
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
ifcs_p3m_assign_fields_ad(ifcs_p3m_data_struct* d,
			  fcs_float *data,
			  fcs_int num_real_particles, 
			  fcs_float* positions,
			  fcs_int shifted,
			  fcs_float* fields);
#endif

/* callback function for near field computations */
static inline void 
ifcs_p3m_compute_near(const void *param, fcs_float dist, fcs_float *f, fcs_float *p)
{
  fcs_float alpha = *((fcs_float *) param);

  fcs_p3m_compute_near(alpha, dist, p, f);
}

/* callback function for performing a whole loop of near field computations (using ifcs_p3m_compute_near) */
FCS_NEAR_LOOP_FP(ifcs_p3m_compute_near_loop, ifcs_p3m_compute_near);

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
void ifcs_p3m_run(void* rd,
		  fcs_int _num_particles,
		  fcs_int _max_num_particles,
		  fcs_float *_positions, 
		  fcs_float *_charges,
		  fcs_float *_fields,
		  fcs_float *_potentials) {
  /* Here we assume, that the method is tuned and that all parameters are valid */

  P3M_INFO(printf( "ifcs_p3m_run() started...\n"));

  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;

  P3M_INFO(printf("    system parameters: box_l=(%"                     \
                  FCS_LMOD_FLOAT "f, %"                                 \
                  FCS_LMOD_FLOAT "f, %"                                 \
                  FCS_LMOD_FLOAT "f)\n",                                \
                  d->box_l[0], d->box_l[1], d->box_l[2]));
  P3M_INFO(printf(                                                      \
                  "    p3m params: r_cut=%" FCS_LMOD_FLOAT              \
                  "f, grid=%d, cao=%d,"                                 \
                  " alpha=%" FCS_LMOD_FLOAT "f, "                       \
                  "grid_off=(%" FCS_LMOD_FLOAT "f,%"                    \
                  FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f)\n",           \
                  d->r_cut, d->grid[0], d->cao, d->alpha,		\
                  d->grid_off[0], d->grid_off[1], d->grid_off[2]));

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

  FCS_P3M_INIT_TIMING(d->comm.mpicomm);

  ifcs_p3m_domain_decompose(d, &gridsort, 
                            _num_particles, _max_num_particles, _positions, _charges,
                            &num_real_particles,
                            &positions, &charges, &indices,
                            &num_ghost_particles,
                            &ghost_positions, &ghost_charges, &ghost_indices);

  /* allocate local fields and potentials */
  fcs_float *fields = NULL; 
  fcs_float *potentials = NULL; 
  if (_fields != NULL)
    fields = malloc(sizeof(fcs_float)*3*num_real_particles);
  if (_potentials != NULL || d->require_total_energy)
    potentials = malloc(sizeof(fcs_float)*num_real_particles);
  
#ifndef P3M_INTERLACE
  FCS_P3M_START_TIMING();
  /* charge assignment */
  ifcs_p3m_assign_charges(d, d->rs_grid, num_real_particles, 
                          positions, charges, 0);
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_charges");

  FCS_P3M_START_TIMING();
  /* gather the ca grid */
  ifcs_p3m_gather_grid(d, d->rs_grid);
  /* now d->rs_grid should contain the local ca grid */
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "gather_grid");
#else
  // First (unshifted) run
  FCS_P3M_START_TIMING();
  /* charge assignment */
  ifcs_p3m_assign_charges(d, d->fft.data_buf, num_real_particles, 
                          positions, charges, 0);
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_charges_1");

  FCS_P3M_START_TIMING();
  /* gather the ca grid */
  ifcs_p3m_gather_grid(d, d->fft.data_buf);
  /* now d->rs_grid should contain the local ca grid */
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "gather_grid_1");
  // Complexify
  for(fcs_int i=d->local_grid.size-1; i>=0; i--)
    d->rs_grid[2*i] = d->fft.data_buf[i];

  // Second (shifted) run
  FCS_P3M_START_TIMING();
  /* charge assignment */
  ifcs_p3m_assign_charges(d, d->fft.data_buf, num_real_particles, 
                          positions, charges, 1);
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_charges_2");

  FCS_P3M_START_TIMING();
  /* gather the ca grid */
  ifcs_p3m_gather_grid(d, d->fft.data_buf);
  /* now d->rs_grid should contain the local ca grid */
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "gather_grid_2");
  // Complexify
  for(fcs_int i=d->local_grid.size-1; i>=0; i--)
    d->rs_grid[2*i+1] = d->fft.data_buf[i];
#endif

  /* output charge assignment grid */
  /* fcs_float* ptr = d->rs_grid; */
  /* char filename[30]; */
  /* sprintf(filename, "charge-%d.vtk", d->comm.rank); */
  /* FILE* outfile = fopen(filename, "w"); */
  /* fprintf(outfile, "# vtk DataFile Version 2.0\ntest\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA %d\nSCALARS charge floats 1\nLOOKUP_TABLE default\n", */
  /* 	  d->local_grid.dim[0], d->local_grid.dim[1], d->local_grid.dim[2], */
  /* 	  d->local_grid.dim[0]* d->local_grid.dim[1]* d->local_grid.dim[2]); */

  /* for (fcs_int k=0; k < d->local_grid.dim[2]; k++) */
  /*   for (fcs_int j=0; j < d->local_grid.dim[1]; j++) */
  /*     for (fcs_int i=0; i < d->local_grid.dim[0]; i++) { */
  /* 	fcs_float charge = d->rs_grid[get_linear_index(i, j, k, d->local_grid.dim)]; */
  /* 	fprintf(outfile, "%e\n", charge); */
  /*     } */
  /* fclose(outfile); */

  /* forward transform */
  FCS_P3M_START_TIMING();
  P3M_DEBUG(printf( "  calling ifcs_fft_perform_forw()...\n"));
  ifcs_fft_perform_forw(&d->fft, &d->comm, d->rs_grid);
  P3M_DEBUG(printf( "  returned from ifcs_fft_perform_forw().\n"));
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "fft_perform_forw");
  
  if (d->require_total_energy || _potentials != NULL) {
    /* apply energy optimized influence function */
    FCS_P3M_START_TIMING();
    ifcs_p3m_apply_energy_influence_function(d);
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "apply_energy_influence_function");
    /* result is in d->ks_grid */
    
    /* compute total energy, but not potentials */
    FCS_P3M_START_TIMING();
    if (d->require_total_energy && potentials == NULL)
      d->total_energy = ifcs_p3m_compute_total_energy(d);
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "compute_total_energy");

    if (_potentials != NULL) {
      /* backtransform the grid */
      FCS_P3M_START_TIMING();
      P3M_DEBUG(printf( "  calling ifcs_fft_perform_back (potentials)...\n"));
      ifcs_fft_perform_back(&d->fft, &d->comm, d->ks_grid);
      P3M_DEBUG(printf( "  returned from ifcs_fft_perform_back.\n"));
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "fft_perform_back (potentials)");

#ifdef P3M_INTERLACE
      /** First (unshifted) run */
      for (fcs_int i=0; i<d->local_grid.size; i++) {
	d->fft.data_buf[i] = d->ks_grid[2*i];
      } 

      FCS_P3M_START_TIMING();
      ifcs_p3m_spread_grid(d, d->fft.data_buf);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid_1");

      FCS_P3M_START_TIMING();
      ifcs_p3m_assign_potentials(d, d->fft.data_buf,
                                 num_real_particles, positions, 
                                 charges, 0, potentials);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_potentials_1");

      /** Second (shifted) run */
      for (fcs_int i=0; i<d->local_grid.size; i++) {
        d->fft.data_buf[i] = d->ks_grid[2*i+1];
      }

      FCS_P3M_START_TIMING();
      ifcs_p3m_spread_grid(d, d->fft.data_buf);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid_2");

      FCS_P3M_START_TIMING();
      ifcs_p3m_assign_potentials(d, d->fft.data_buf,
                                 num_real_particles, positions,
                                 charges, 1, potentials);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_potentials_2");
#else
      /* redistribute energy grid */
      FCS_P3M_START_TIMING();
      ifcs_p3m_spread_grid(d, d->ks_grid);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid");
      /* compute potentials */
      FCS_P3M_START_TIMING();
      ifcs_p3m_assign_potentials(d, d->ks_grid,
                                 num_real_particles, positions, 
                                 charges, 0,
                                 potentials);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_potentials");
#endif
    }
  }

  if (_fields != NULL) {
    /* apply force optimized influence function */
    FCS_P3M_START_TIMING();
    ifcs_p3m_apply_force_influence_function(d);
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "apply_force_influence_function");
    
#ifdef P3M_AD
#ifdef P3M_INTERLACE
    /* backtransform the grid */
    FCS_P3M_START_TIMING();
    P3M_DEBUG(printf( "  calling ifcs_fft_perform_back...\n"));
    ifcs_fft_perform_back(&d->fft, &d->comm, d->ks_grid);
    P3M_DEBUG(printf( "  returned from ifcs_fft_perform_back.\n"));
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "fft_perform_back");
    
    /* First (unshifted) run */
    P3M_INFO(printf("  computing unshifted grid\n"));
    for (fcs_int i=0; i<d->local_grid.size; i++) {
      d->fft.data_buf[i] = d->ks_grid[2*i];
    } 
    
    FCS_P3M_START_TIMING();
    ifcs_p3m_spread_grid(d, d->fft.data_buf);
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid");
    
    FCS_P3M_START_TIMING();
    ifcs_p3m_assign_fields_ad(d, d->fft.data_buf, num_real_particles, 
                              positions, 0, fields);
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_fields");
    
    /* Second (shifted) run */
    P3M_INFO(printf("  computing shifted grid\n"));
    for (fcs_int i=0; i<d->local_grid.size; i++) {
      d->fft.data_buf[i] = d->ks_grid[2*i+1];
    }
    
    FCS_P3M_START_TIMING();
    ifcs_p3m_spread_grid(d, d->fft.data_buf);
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid");
    
    FCS_P3M_START_TIMING();
    ifcs_p3m_assign_fields_ad(d, d->fft.data_buf, num_real_particles, 
                              positions, 1, fields);
    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign_fields");
#endif /* P3M_INTERLACE */
#endif /* P3M_AD */
  
#ifdef P3M_IK
    /* result is in d->ks_grid */
    for (int dim = 0; dim < 3; dim++) {
      /* differentiate in direction dim */
      /* result is stored in d->rs_grid */
      FCS_P3M_START_TIMING();
      ifcs_p3m_ik_diff(d, dim);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "ik_diff");
      /* backtransform the grid */
      FCS_P3M_START_TIMING();
      P3M_DEBUG(printf( "  calling ifcs_fft_perform_back (field dim=%d)...\n", dim));
      ifcs_fft_perform_back(&d->fft, &d->comm, d->rs_grid);
      P3M_DEBUG(printf( "  returned from ifcs_fft_perform_back.\n"));
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "fft_perform_back");
      
#ifdef P3M_INTERLACE
      /** First (unshifted) run */
      for(fcs_int i=0;i<d->local_grid.size;i++) {
        d->fft.data_buf[i] = d->rs_grid[2*i];
      } 
      
      FCS_P3M_START_TIMING();
      ifcs_p3m_spread_grid(d, d->fft.data_buf);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid");
      
      FCS_P3M_START_TIMING();
      ifcs_p3m_assign_fields_ik(d, d->fft.data_buf, dim, num_real_particles, 
                                positions, 0, fields);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign fields");
      
      /** Second (shifted) run */
      for(fcs_int i=0;i<d->local_grid.size;i++) {
        d->fft.data_buf[i] = d->rs_grid[2*i+1];
    } 
      
      FCS_P3M_START_TIMING();
      ifcs_p3m_spread_grid(d, d->fft.data_buf);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid");
      
      FCS_P3M_START_TIMING();
      ifcs_p3m_assign_fields_ik(d, d->fft.data_buf, dim, num_real_particles, 
                                positions, 1, fields);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign fields");
#else
      /* redistribute force grid */
      FCS_P3M_START_TIMING();
      ifcs_p3m_spread_grid(d, d->rs_grid);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "spread_grid");
      
      /* Output force grid */
      /* fcs_int i,j,k; */
      /* fcs_float* ptr = d->rs_grid; */
      /* char filename[30]; */
      /* sprintf(filename, "field-%d.vtk", dim); */
      /* FILE* outfile = fopen(filename, "w"); */
      /* fprintf(outfile, "# vtk DataFile Version 2.0\ntest\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA %d\nSCALARS OutArray  floats 1\nLOOKUP_TABLE default\n", */
      /* 	      d->local_grid.dim[0], d->local_grid.dim[1], d->local_grid.dim[2], */
      /* 	      d->local_grid.dim[0]* d->local_grid.dim[1]* d->local_grid.dim[2]); */
      
      /* for (k=0; k < d->local_grid.dim[2]; k++) */
      /* 	for (j=0; j < d->local_grid.dim[1]; j++) */
      /* 	  for (i=0; i < d->local_grid.dim[0]; i++) { */
      /* 	    fcs_float charge = d->rs_grid[get_linear_index(i, j, k, d->local_grid.dim)]; */
      /* 	    fprintf(outfile, "%e\n", charge); */
      /* 	  } */
      /* fclose(outfile); */
      
      FCS_P3M_START_TIMING();
      ifcs_p3m_assign_fields_ik(d, d->rs_grid, dim, num_real_particles, 
                                positions, 0, fields);
      FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "assign fields");
#endif /* P3M_INTERLACE */
    }
#endif /* P3M_IK */
  }   
 
  if (d->near_field_flag) {
    FCS_P3M_START_TIMING();
    /* compute near field */
    fcs_near_t near;
    fcs_float alpha = d->alpha;
  
    fcs_near_create(&near);
    /*  fcs_near_set_field_potential(&near, ifcs_p3m_compute_near);*/
    fcs_near_set_loop(&near, ifcs_p3m_compute_near_loop);

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

    FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "near field calculation");
  }

  /* sort particles back */
  P3M_DEBUG(printf( "  calling fcs_gridsort_sort_backward()...\n"));
  FCS_P3M_START_TIMING();
  fcs_gridsort_sort_backward(&gridsort,
                             fields, potentials,
                             _fields, _potentials, 1,
                             d->comm.mpicomm);
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "backward sort");
  P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_backward().\n"));
  
  fcs_gridsort_free(&gridsort);
  fcs_gridsort_destroy(&gridsort);

  if (fields != NULL) free(fields);
  if (potentials != NULL) free(potentials);

  P3M_INFO(printf( "ifcs_p3m_run() finished.\n"));
}

/***************************************************/
/* RUN COMPONENTS */
static void
ifcs_p3m_domain_decompose(ifcs_p3m_data_struct *d, fcs_gridsort_t *gridsort,
                          fcs_int _num_particles, fcs_int _max_num_particles, 
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

  FCS_P3M_INIT_TIMING(d->comm.mpicomm);

  fcs_gridsort_create(gridsort);
  
  fcs_gridsort_set_system(gridsort, box_base, box_a, box_b, box_c, NULL);
  fcs_gridsort_set_particles(gridsort, _num_particles, _max_num_particles, _positions, _charges);

  P3M_DEBUG(printf( "  calling fcs_gridsort_sort_forward()...\n"));
  FCS_P3M_START_TIMING();
  /* @todo: Set skin to r_cut only, when near field is wanted! */
  fcs_gridsort_sort_forward(gridsort,
        		    (d->near_field_flag ? d->r_cut : 0.0),
        		    d->comm.mpicomm);
  FCS_P3M_FINISH_TIMING(d->comm.mpicomm, "forward sort");
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
/** Assign the charges to the grid */
static void 
ifcs_p3m_assign_charges(ifcs_p3m_data_struct* d,
			fcs_float *data,
			fcs_int num_real_particles,
			fcs_float *positions, 
			fcs_float *charges,
			fcs_int shifted) {

  P3M_DEBUG(printf( "  ifcs_p3m_assign_charges() started...\n"));
  /* init local charge grid */
  for (fcs_int i=0; i<d->local_grid.size; i++) data[i] = 0.0;

  /* now assign the charges */
  for (fcs_int pid=0; pid < num_real_particles; pid++)
    ifcs_p3m_assign_single_charge(d, data, pid, &positions[pid*3], 
                                  charges[pid], shifted);

  P3M_DEBUG(printf( "  ifcs_p3m_assign_charges() finished...\n"));
}

/** Compute the data of the charge assignment grid points.

    The function returns the linear index of the top left grid point
    in the charge assignment grid that corresponds to real_pos. When
    "shifted" is set, it uses the shifted position for interlacing.
    After the call, "caf_ptr" contain pointers into fcs_float arrays
    that contain the value of the charge assignment fraction (caf) for
    x,y,z.
 */
static fcs_int 
ifcs_get_ca_points(ifcs_p3m_data_struct *d, 
                   fcs_float real_pos[3], 
                   fcs_int shifted, 
                   fcs_float *caf_ptr[6]) {
  const fcs_int cao = d->cao;
  const fcs_int interpol = d->n_interpol;

  /* linear index of the grid point */
  fcs_int linind = 0;

  for (fcs_int dim=0; dim<3; dim++) {
    /* position in normalized coordinates in [0,1] */
    fcs_float pos    = (real_pos[dim] - d->local_grid.ld_pos[dim]) * d->ai[dim];
    /* shift position to the corner of the charge assignment area */
    pos -= d->pos_shift;
    /* if using the interlaced grid, shift it more */
    if (shifted) pos -= 0.5;
    /* nearest grid point in the ca grid */
    fcs_int grid_ind  = (fcs_int)floor(pos);
    /* linear index of grid point */
    linind = grid_ind + linind*d->local_grid.dim[dim];

    if (interpol > 0) {
      /* index of nearest point in the interpolation grid */
      fcs_int intind = (fcs_int) ((pos-grid_ind) * (2*d->n_interpol + 1));
      /* generate pointer into interpolation grids */
      caf_ptr[dim] = &d->int_caf[intind*cao];
#ifdef P3M_AD
      caf_ptr[3+dim] = &d->int_caf_d[intind*cao];
#endif
    } else {
      /* distance between position and nearest grid point */
      fcs_float dist = (pos-grid_ind)-0.5;
      /* pointer to the interpolation grid, which is now used to
         store the precomputed values of caf. */
      caf_ptr[dim] = &d->int_caf[dim*cao];
      /* fill array */
      for (fcs_int i = 0; i < cao; i++)
        caf_ptr[dim][i] = ifcs_p3m_caf(i, dist, cao);
#ifdef P3M_AD
      caf_ptr[3+dim] = &d->int_caf_d[dim*cao];
      for (fcs_int i = 0; i < cao; i++)
        caf_ptr[3+dim][i] = ifcs_p3m_caf_d(i, dist, cao);
#endif
    }

#ifdef ADDITIONAL_CHECKS
    if (real_pos[dim] < d->comm.my_left[dim] 
        || real_pos[dim] > d->comm.my_right[dim]) {
      printf("%d: dim %d: position not in domain! (%" FCS_LMOD_FLOAT    \
             "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f)\n", 
             d->comm.rank, dim, real_pos[0], real_pos[1], real_pos[2]);
    }
#endif
  }

  return linind;
}

/** Assign a single charge to the grid. If pid <=0, the charge
    fractions are not stored. This can be used to add virtual charges
    to the system, where no field is computed. This can be used for e.g. ELC or ICC*. */
static void 
ifcs_p3m_assign_single_charge(ifcs_p3m_data_struct *d,
			      fcs_float *data,
			      fcs_int pid,
			      fcs_float real_pos[3],
			      fcs_float q,
			      fcs_int shifted) {
  const fcs_int cao = d->cao;
  const fcs_int q2off = d->local_grid.q_2_off;
  const fcs_int q21off = d->local_grid.q_21_off;
  
  fcs_float *caf_begin[6];
  fcs_int linind_grid = 
    ifcs_get_ca_points(d, real_pos, shifted, caf_begin);
  fcs_float *caf_end[3]  = { caf_begin[0]+cao,
                             caf_begin[1]+cao,
                             caf_begin[2]+cao };

  /* Loop over all ca grid points nearby and compute charge assignment fraction */
  for (fcs_float *caf_x = caf_begin[0]; caf_x != caf_end[0]; caf_x++) {
    for (fcs_float *caf_y = caf_begin[1]; caf_y != caf_end[1]; caf_y++) {
      fcs_float caf_xy = *caf_x * *caf_y;
      for (fcs_float *caf_z = caf_begin[2]; caf_z != caf_end[2]; caf_z++) {
        /* add it to the grid */
        data[linind_grid] += q * caf_xy * *caf_z;
        linind_grid++;
      }
      linind_grid += q2off;
    }
    linind_grid += q21off;
  } 
}

/* Gather information for FFT grid inside the nodes domain (inner local grid) */
static void ifcs_p3m_gather_grid(ifcs_p3m_data_struct* d, fcs_float* rs_grid) {
  MPI_Status status;
  fcs_float *tmp_ptr;

  P3M_DEBUG(printf( "  ifcs_p3m_gather_grid() started...\n"));

  /* direction loop */
  for(fcs_int s_dir=0; s_dir<6; s_dir++) {
    fcs_int r_dir;
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /* pack send block */
    if (d->sm.s_size[s_dir]>0)
      ifcs_fft_pack_block(rs_grid, d->send_grid, d->sm.s_ld[s_dir], 
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
      ifcs_p3m_add_block(d->recv_grid, rs_grid, d->sm.r_ld[r_dir], 
			 d->sm.r_dim[r_dir], d->local_grid.dim);
    }
  }

  P3M_DEBUG(printf( "  ifcs_p3m_gather_grid() finished.\n"));
}

static void ifcs_p3m_add_block(fcs_float *in, fcs_float *out, int start[3], int size[3], int dim[3]) {
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
static void ifcs_p3m_apply_energy_influence_function(ifcs_p3m_data_struct* d) {
  P3M_DEBUG(printf( "  ifcs_p3m_apply_energy_influence_function() started...\n"));
  const fcs_int size = d->fft.plan[3].new_size;
  for (fcs_int i=0; i < size; i++) {
    d->ks_grid[2*i] = d->g_energy[i] * d->rs_grid[2*i]; 
    d->ks_grid[2*i+1] = d->g_energy[i] * d->rs_grid[2*i+1]; 
  }
  P3M_DEBUG(printf( "  ifcs_p3m_apply_energy_influence_function() finished.\n"));
}

/* apply the influence function */
static void ifcs_p3m_apply_force_influence_function(ifcs_p3m_data_struct* d) {
  P3M_DEBUG(printf( "  ifcs_p3m_apply_force_influence_function() started...\n"));
  const fcs_int size = d->fft.plan[3].new_size;
  for (fcs_int i=0; i < size; i++) {
    d->ks_grid[2*i] = d->g_force[i] * d->rs_grid[2*i]; 
    d->ks_grid[2*i+1] = d->g_force[i] * d->rs_grid[2*i+1]; 
  }
  P3M_DEBUG(printf( "  ifcs_p3m_apply_force_influence_function() finished.\n"));
}

static void ifcs_p3m_ik_diff(ifcs_p3m_data_struct* d, int dim) {
  fcs_int ind;
  fcs_int j[3];
  fcs_float* d_operator = NULL;
  /* direction in k space: */
  fcs_int dim_rs = (dim+d->ks_pnum)%3;

  P3M_DEBUG(printf( "  ifcs_p3m_ik_diff() started...\n"));
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
	  -2.0*FCS_PI*(d->ks_grid[ind+1] * d_operator[ j[dim]+d->fft.plan[3].start[dim] ])
	  / d->box_l[dim_rs];
	d->rs_grid[ind+1] =
	  2.0*FCS_PI*d->ks_grid[ind] * d_operator[ j[dim]+d->fft.plan[3].start[dim] ]
	  / d->box_l[dim_rs];
	ind+=2;
      }
    }
  }

  P3M_DEBUG(printf( "  ifcs_p3m_ik_diff() finished.\n"));
  /* store the result in d->rs_grid */
}

 static void ifcs_p3m_spread_grid(ifcs_p3m_data_struct* d, fcs_float* rs_grid) {
   int s_dir,r_dir,evenodd;
   MPI_Status status;
   fcs_float *tmp_ptr;
  P3M_DEBUG(printf( "  ifcs_p3m_spread_grid() started...\n"));
  
   /* direction loop */
   for(s_dir=5; s_dir>=0; s_dir--) {
     if(s_dir%2==0) r_dir = s_dir+1;
     else           r_dir = s_dir-1;
     /* pack send block */ 
     if(d->sm.s_size[s_dir]>0) 
       ifcs_fft_pack_block(rs_grid, d->send_grid, d->sm.r_ld[r_dir], d->sm.r_dim[r_dir], d->local_grid.dim, 1);
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
       ifcs_fft_unpack_block(d->recv_grid, rs_grid, d->sm.s_ld[s_dir], d->sm.s_dim[s_dir], d->local_grid.dim, 1); 
     }
   }

  P3M_DEBUG(printf( "  ifcs_p3m_spread_grid() finished.\n"));
}


/** Compute the total energy of the system in kspace. No need to
    backtransform the FFT grid in this case! */
static fcs_float ifcs_p3m_compute_total_energy(ifcs_p3m_data_struct* d) {
  fcs_float local_k_space_energy;
  fcs_float k_space_energy;

  P3M_DEBUG(printf( "  ifcs_p3m_compute_total_energy() started...\n"));

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
  k_space_energy -= d->square_sum_q * FCS_PI * prefactor / SQR(d->alpha);

  P3M_DEBUG(printf( "  ifcs_p3m_compute_total_energy() finished.\n"));
  return k_space_energy;
}

/* Backinterpolate the potentials obtained from k-space to the positions */
static void 
ifcs_p3m_assign_potentials(ifcs_p3m_data_struct* d, 
			   fcs_float *data,
                           fcs_int num_real_particles, 
                           fcs_float* positions, fcs_float* charges, 
                           fcs_int shifted,
                           fcs_float* potentials) {
  const fcs_int cao = d->cao;
  const fcs_int q2off = d->local_grid.q_2_off;
  const fcs_int q21off = d->local_grid.q_21_off;
  const fcs_float prefactor = 1.0 / (d->box_l[0] * d->box_l[1] * d->box_l[2]);
  
  P3M_DEBUG(printf( "  ifcs_p3m_assign_potentials() started...\n"));
  /* Loop over all particles */
  for (fcs_int pid=0; pid < num_real_particles; pid++) {
    fcs_float potential = 0.0;

    fcs_float *caf_begin[6];
    fcs_int linind_grid = 
      ifcs_get_ca_points(d, &positions[pid*3], shifted, caf_begin);
    fcs_float *caf_end[3]  = { caf_begin[0]+cao,
                               caf_begin[1]+cao,
                               caf_begin[2]+cao };
    
    /* Loop over all ca grid points nearby and compute charge assignment fraction */
    for (fcs_float *caf_x = caf_begin[0]; caf_x != caf_end[0]; caf_x++) {
      for (fcs_float *caf_y = caf_begin[1]; caf_y != caf_end[1]; caf_y++) {
        fcs_float caf_xy = *caf_x * *caf_y;
        for (fcs_float *caf_z = caf_begin[2]; caf_z != caf_end[2]; caf_z++) {
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
  P3M_DEBUG(printf( "  ifcs_p3m_assign_potentials() finished.\n"));
}

/* Backinterpolate the forces obtained from k-space to the positions */
static void 
ifcs_p3m_assign_fields_ik(ifcs_p3m_data_struct* d, 
                          fcs_float *data,
                          fcs_int dim,
                          fcs_int num_real_particles,
                          fcs_float* positions,
                          fcs_int shifted,
                          fcs_float* fields) {
  const fcs_int cao = d->cao;
  const fcs_int q2off = d->local_grid.q_2_off;
  const fcs_int q21off = d->local_grid.q_21_off;
  const fcs_float prefactor = 1.0 / (2.0 * d->box_l[0] * d->box_l[1] * d->box_l[2]);
  const fcs_int dim_rs = (dim+d->ks_pnum) % 3;

  P3M_DEBUG(printf( "  ifcs_p3m_assign_fields() started...\n"));
  /* Loop over all particles */
  for (fcs_int pid=0; pid < num_real_particles; pid++) {
    fcs_float field = 0.0;
    
    fcs_float *caf_begin[6];
    fcs_int linind_grid = 
      ifcs_get_ca_points(d, &positions[3*pid], shifted, caf_begin);
    fcs_float *caf_end[3]  = { caf_begin[0]+cao,
                               caf_begin[1]+cao,
                               caf_begin[2]+cao };
    
    /* loop over the local grid, compute the field */
    for (fcs_float *caf_x = caf_begin[0]; caf_x != caf_end[0]; caf_x++) {
      for (fcs_float *caf_y = caf_begin[1]; caf_y != caf_end[1]; caf_y++) {
        fcs_float caf_xy = *caf_x * *caf_y;
        for (fcs_float *caf_z = caf_begin[2]; caf_z != caf_end[2]; caf_z++) {
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
  P3M_DEBUG(printf( "  ifcs_p3m_assign_fields() finished.\n"));
}

#ifdef P3M_AD
/* Backinterpolate the forces obtained from k-space to the positions */
static void 
ifcs_p3m_assign_fields_ad(ifcs_p3m_data_struct* d,
			  fcs_float *data,
			  fcs_int num_real_particles, 
			  fcs_float* positions,
                          fcs_int shifted,
			  fcs_float* fields) {
  const fcs_int cao = d->cao;
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

  fcs_float field[3] = { 0.0, 0.0, 0.0 };

  P3M_DEBUG(printf( "  ifcs_p3m_assign_fields() [AD] started...\n"));
  /* Loop over all particles */
  for (fcs_int pid=0; pid < num_real_particles; pid++) {
    field[0] = field[1] = field[2] = 0.0;

    fcs_float *caf_begin[6];
    fcs_int linind_grid = 
      ifcs_get_ca_points(d, &positions[pid*3], shifted, caf_begin);
    fcs_float *caf_end[3]  = { caf_begin[0]+cao,
                               caf_begin[1]+cao,
                               caf_begin[2]+cao };

    fcs_float *caf_x_d = caf_begin[3];
    for (fcs_float *caf_x = caf_begin[0]; caf_x != caf_end[0]; caf_x++) {
      fcs_float *caf_y_d = caf_begin[4];
      for (fcs_float *caf_y = caf_begin[1]; caf_y != caf_end[1]; caf_y++) {
        fcs_float *caf_z_d = caf_begin[5];
        for (fcs_float *caf_z = caf_begin[2]; caf_z != caf_end[2]; caf_z++) {
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
  P3M_DEBUG(printf( "  ifcs_p3m_assign_fields() finished.\n"));
}
#endif
