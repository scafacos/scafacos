 /*
  Copyright (C) 2014 Olaf Lenz, Gabriel Sichardt
  Copyright (C) 2013 Olaf Lenz, Florian Weik, Gabriel Sichardt
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
#include "run.hpp"
#include "types.hpp"
#include "utils.hpp"
#include "caf.hpp"
#include "tune.hpp"
#include "fcs_p3m_p.h"
#include "common/gridsort/gridsort.h"
#include "common/near/near.h"
#include <cassert>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <stdio.h>
#include "prepare.hpp"

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

/* sum and collect timing information.*/
static void collect_print_timings(ifcs_p3m_data_struct *d);
static fcs_int triclinic_check(fcs_float* tricl, fcs_int number);

static fcs_float * position_store(fcs_float positions, fcs_int number);

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
//todo: move methods somewhere else (bottom, utils, triclinic-utils,...)
/*add up the measured timings and collect information from all procs.*/
static void collect_print_timings(ifcs_p3m_data_struct *d) {
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
            d->timings[TIMING_FAR] / d->timings[TIMING]);
    printf("     near=%le (%lf)\n", d->timings[TIMING_NEAR],
            d->timings[TIMING_NEAR] / d->timings[TIMING]);
    printf("       ca=%le (%lf)\n", d->timings[TIMING_CA],
            d->timings[TIMING_CA] / d->timings[TIMING]);
    printf("      pot=%le (%lf)\n", d->timings[TIMING_POTENTIALS],
            d->timings[TIMING_POTENTIALS] / d->timings[TIMING]);
    printf("   fields=%le (%lf)\n", d->timings[TIMING_FIELDS],
            d->timings[TIMING_FIELDS] / d->timings[TIMING]);
    printf("   gather=%le (%lf)\n", d->timings[TIMING_GATHER],
            d->timings[TIMING_GATHER] / d->timings[TIMING]);
    printf("   spread=%le (%lf)\n", d->timings[TIMING_SPREAD],
            d->timings[TIMING_SPREAD] / d->timings[TIMING]);
    printf("  forward=%le (%lf)\n", d->timings[TIMING_FORWARD],
            d->timings[TIMING_FORWARD] / d->timings[TIMING]);
    printf("     back=%le (%lf)\n", d->timings[TIMING_BACK],
            d->timings[TIMING_BACK] / d->timings[TIMING]);
    printf("   decomp=%le (%lf)\n", d->timings[TIMING_DECOMP],
            d->timings[TIMING_DECOMP] / d->timings[TIMING]);
    printf("     comp=%le (%lf)\n", d->timings[TIMING_COMP],
            d->timings[TIMING_COMP] / d->timings[TIMING]);
#endif

}
//todo: method description
//todo: does not work with [][] (arbitrary dimensions).
static void print_matrix(fcs_float  matrix[][3],fcs_int dimhor, fcs_int dimver){
    fcs_int counterhor, counterver;
    printf("printmatrix\n");
    for(counterhor=0;counterhor<dimhor;counterhor++){
        for(counterver=0;counterver<dimver;counterver++){
            printf("%f ", matrix[counterhor][counterver]);
        }printf("\n");
    }
}
//todo: method description
static fcs_float* to_triclinic(ifcs_p3m_data_struct *d, fcs_float *vectors, fcs_int number){
    if(d->cosy_flag==triclinic){
        printf("to_triclinic detected triclinic flag! \n");
        if(triclinic_check(vectors, number)==0)
        return vectors;
        else
            printf("but the vector cannot be triclinic. trying to convert!\n");
    }
    
    fcs_float* tricl;
    tricl=static_cast<fcs_float*>(malloc(3*number*sizeof(fcs_float)));
 //   print_matrix(d->box_matrix,3,3);
  int i;
  for (i = 0; i < number; i++) {
        printf("to_tricl input vectors: %f %f %f \n", i, vectors[3 * i], vectors[3 * i + 1], vectors[3 * i + 2]);
    }
        for (i = 0; i < number; i++) {
            tricl[3 * i] = 1 / d->box_matrix[0][0] * vectors[3 * i] - d->box_matrix[1][0] / (d->box_matrix[0][0] * d->box_matrix[1][1]) * vectors[3 * i + 1] + (d->box_matrix[1][0]*d->box_matrix[2][1]-d->box_matrix[1][1]*d->box_matrix[2][0])/(d->box_matrix[0][0]*d->box_matrix[1][1]*d->box_matrix[2][2])*vectors[3 * i + 2];
            tricl[3 * i + 1] = 1 / d->box_matrix[1][1] * vectors[3 * i + 1] - d->box_matrix[2][1] / (d->box_matrix[1][1] * d->box_matrix[2][2]) * vectors[3 * i + 2];
            tricl[3 * i + 2] = 1 / d->box_matrix[2][2] * vectors[3 * i + 2];
        }
    
  //todo: move this to the calling method
   // d->box_l[0]=d->box_l[1]=d->box_l[2]=1.0;

#ifdef ADDITIONAL_CHECKS
  if(triclinic_check(tricl, number)==1)MPI_Abort(MPI_COMM_WORLD,1);
#endif
  //todo: remove printing as soon as all works
    printf("to_tricl: vectors are now in triclinic coordinates: \n");


    for (i = 0; i < number; i++) {
        printf("tricl: %f %f %f \n", i, tricl[3 * i], tricl[3 * i + 1], tricl[3 * i + 2]);
    }
    d->cosy_flag=triclinic;
    return tricl;
}
/*checks if the given vector tricl can be triclinic. simply checks whether all positions are below 1 and above 0.
 * 
 * @param tricl fcs_float pointer to an array with 3*number floats where 3 in a row are the position of one particle
 * 
 * @param number the number of particles, i.e. 1/3 of the number fcs_floats in tricl
 * 
 * @return 1 if the positions cannot be triclinic, 0 if the positions might be triclinic
 *
 */
static fcs_int triclinic_check(fcs_float* tricl, fcs_int number) {
    fcs_int i;
    fcs_int res = 0;
    for (i = 0; i < number; i++) {
        if (tricl[3 * i] > 1 || tricl[3 * i] < 0 || tricl[3 * i + 1] > 1 || tricl[3 * i + 1] < 0 || tricl[3 * i + 2] > 1 || tricl[3 * i + 2] < 0) {
            printf("ERROR: this cannot be a triclinic position: particle %d, coordinate: %f %f %f \n", i, tricl[3 * i], tricl[3 * i + 1], tricl[3 * i + 2]);
            res = 1;
        }
    }
    return res;
}
/* convert the given triclinic vectors to cartesian vectors. This method is not
 * usable for the force conversion!
 * 
 * @param d the datastruct that contains system information like the box matrix.
 * 
 * @param vectors a pointer to a fcs_float array containing vectors where indexes
 * 3*n , 3*n+1 and 3*n+2 (n is an integer) belong to the same vector and are
 * sorted according to the box vector order
 * 
 * @param number the number of vectors contained in vectors. i.e. vectors
 * contains 3*number fcs_floats.
 * 
 * @return a pointer to a fcs_float array containing the cartesian vectors.
 */
static fcs_float* to_cartesian(ifcs_p3m_data_struct *d,fcs_float* vectors, fcs_int number){
    if(d->cosy_flag==cartesian) return vectors;
    fcs_int counter;
    fcs_float* cart = static_cast<fcs_float*>(malloc(3*number*sizeof(fcs_float)));
    for(counter=0;counter<number; counter++){
        cart[3*counter]=vectors[3*counter]*d->box_matrix[0][0]+vectors[3*counter+1]*d->box_matrix[1][0]+vectors[3*counter+2]*d->box_matrix[2][0];
        cart[3*counter+1]=vectors[3*counter+1]*d->box_matrix[1][1]+vectors[3*counter+2]*d->box_matrix[2][1];
        cart[3*counter+2]=vectors[3*counter+2]*d->box_matrix[2][2];
    }
    d->cosy_flag=cartesian;
    return cart;
}

static void print_vector(fcs_float* vector, fcs_int number){
    fcs_int cnt;
    for(cnt=0; cnt<number; cnt++){
        printf("vector %d : %f %f %f\n",cnt,vector[3*cnt],vector[3*cnt+1],vector[3*cnt+2]);
    }
}
static void print_array(fcs_float* array, fcs_int number){
    fcs_int cnt;
    for(cnt=0; cnt<number; cnt++){
        printf("array[%d] = %f, ",cnt,array[cnt]);
    }
    printf("\n");
}
static fcs_float* cartesian_field(ifcs_p3m_data_struct *d, fcs_float* triclinic_field, fcs_int number) {
    fcs_float* cart = static_cast<fcs_float*> (malloc(3 * number * sizeof (fcs_float)));
    
    fcs_int part_no;
    
    // the next lines do: //UEBERPRUEFUNG NOETIG!!!!!
//      for (part_no = 0; part_no < number; part_no++){
//          
//          cart[3*part_no]=2*triclinic_field[3*part_no+2]+2*triclinic_field[3*part_no+1]+4*triclinic_field[3*part_no];
//          cart[3*part_no+1]=2*triclinic_field[3*part_no+2]+2*triclinic_field[3*part_no+1];
//          cart[3*part_no+2]=2*triclinic_field[3*part_no+2];
//          
//      }
//    
    //the next lines do F = M^-1^T F where M transforms from tric to cart

    
    for (part_no = 0; part_no < number; part_no++) {
        cart[3 * part_no + 2] = (d->box_matrix[1][0] * d->box_matrix[2][1] - d->box_matrix[1][1] * d->box_matrix[2][0]) / (d->box_matrix[0][0] * d->box_matrix[1][1] * d->box_matrix[2][2]) * triclinic_field[3 * part_no]-(d->box_matrix[2][1]) / (d->box_matrix[1][1] * d->box_matrix[2][2]) * triclinic_field[3 * part_no + 1] + 1 / d->box_matrix[2][2] * triclinic_field[3 * part_no + 2];
        cart[3 * part_no + 1] = -(d->box_matrix[1][0]) / (d->box_matrix[0][0] * d->box_matrix[1][1]) * triclinic_field[3 * part_no] + 1 / d->box_matrix[1][1] * triclinic_field[3 * part_no + 1];
        cart[3 * part_no] = 1 / (d->box_matrix[0][0]) * triclinic_field[3 * part_no];
    }

    
    //the following lines do: F=CT F where CT transforms from cart to tric.
//    for (part_no = 0; part_no < number; part_no++) {
//        cart[3 * part_no ] = 1 / d->box_matrix[0][0] * triclinic_field[3 * part_no] - d->box_matrix[1][0] / (d->box_matrix[0][0] * d->box_matrix[1][1]) * triclinic_field[3 * part_no + 1]+((d->box_matrix[1][0] * d->box_matrix[2][1])-(d->box_matrix[1][1] * d->box_matrix[2][0])) / (d->box_matrix[0][0] * d->box_matrix[1][1] * d->box_matrix[2][2]) * triclinic_field[3 * part_no + 2];
//        cart[3 * part_no + 1] = 1 / d->box_matrix[1][1] * triclinic_field[3 * part_no + 1]-(d->box_matrix[2][1]) / (d->box_matrix[1][1] * d->box_matrix[2][2]) * triclinic_field[3 * part_no + 2];
//        cart[3 * part_no + 2] = 1 / d->box_matrix[2][2] * triclinic_field[3 * part_no + 2];
//    }
    
    printf("field transformed\n");
    return cart;
}


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
#define START(ID)                                    \
  if (d->require_timings) d->timings[ID] += -MPI_Wtime();
#define STOP(ID)                                     \
  if (d->require_timings) d->timings[ID] += MPI_Wtime();
#define STOPSTART(ID1, ID2)                            \
  if (d->require_timings) {                                       \
    d->timings[ID1] += MPI_Wtime();                               \
    d->timings[ID2] += -MPI_Wtime();                              \
  }                                                               \

void ifcs_p3m_circumvent_tuning(void* rd){
ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
d->alpha = 1.494582;
d->cao = 7;
d->grid[0]=46; d->grid[1]=32; d->grid[2]=40;
d->r_cut=2.828427;
ifcs_p3m_prepare(d,5);
}

#if defined(P3M_INTERLACE) && defined(P3M_AD)
void ifcs_p3m_run(void* rd,
		  fcs_int _num_particles,
		  fcs_int _max_num_particles,
		  fcs_float *_positions, 
		  fcs_float *_charges,
		  fcs_float *_fields,
		  fcs_float *_potentials) {
  /* Here we assume, that the method is already tuned and that all
     parameters are valid */
  P3M_INFO(printf( "ifcs_p3m_run() started...\n"));
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;
  fcs_float* _positions_triclinic = static_cast<fcs_float*>(malloc(3*_num_particles*sizeof(fcs_float)));
  _positions_triclinic=to_triclinic(d,_positions, _num_particles);
    if(triclinic_check(_positions_triclinic, _num_particles)==1){
        //ugly way of getting the right cosy
     while(triclinic_check(_positions_triclinic, _num_particles)==1){   
         d->cosy_flag==cartesian;
         _positions_triclinic=to_triclinic(d,_positions_triclinic, _num_particles);
     }
    }
  /* reset all timers */
  if (d->require_timings) {
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
  printf("positions triclinic before decompose\n");
  print_vector(_positions_triclinic, _num_particles);
  START(TIMING_DECOMP)
  ifcs_p3m_domain_decompose(d, &gridsort, 
                            _num_particles, _max_num_particles, _positions_triclinic, _charges,
                            &num_real_particles,
                            &positions, &charges, &indices,
                            &num_ghost_particles,
                            &ghost_positions, &ghost_charges, &ghost_indices);
//  printf("positions after decompose:\n");
//  print_vector(positions,num_real_particles);
  /* allocate local fields and potentials */
  fcs_float *fields = NULL; 
  fcs_float *potentials = NULL; 
  fcs_float *nearfields = NULL; 
  fcs_float *nearpotentials = NULL; 
  
  if (_fields != NULL){
    fields = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*3*num_real_particles));
    nearfields = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*3*num_real_particles));
  }
  if (_potentials != NULL || d->require_total_energy){
    potentials = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*num_real_particles));
    nearpotentials = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*num_real_particles));
  }
 
  STOPSTART(TIMING_DECOMP, TIMING_CA);

  /* charge assignment */
  ifcs_p3m_assign_charges(d, d->fft.data_buf, num_real_particles, 
                          positions, charges, 0);
  STOPSTART(TIMING_CA, TIMING_GATHER);
  /* gather the ca grid */
  ifcs_p3m_gather_grid(d, d->fft.data_buf);

  // Complexify
  for (fcs_int i = d->local_grid.size-1; i >= 0; i--)
    d->rs_grid[2*i] = d->fft.data_buf[i];

  STOPSTART(TIMING_GATHER, TIMING_CA);

  // Second (shifted) run
  /* charge assignment */
  ifcs_p3m_assign_charges(d, d->fft.data_buf, num_real_particles, 
                          positions, charges, 1);

  STOPSTART(TIMING_CA, TIMING_GATHER);

  /* gather the ca grid */
  ifcs_p3m_gather_grid(d, d->fft.data_buf);
  /* now d->rs_grid should contain the local ca grid */
  
  // Complexify
  for (fcs_int i = d->local_grid.size-1; i >= 0; i--){    
    d->rs_grid[2*i+1] = d->fft.data_buf[i];
//     printf("rsgrid %f at i=%d\n",d->rs_grid[2*i+1],i);
  }
  STOP(TIMING_GATHER);

  /* forward transform */
  START(TIMING_FORWARD);
  P3M_DEBUG(printf( "  calling ifcs_fft_perform_forw()...\n"));
  ifcs_fft_perform_forw(&d->fft, &d->comm, d->rs_grid);
  P3M_DEBUG(printf( "  returned from ifcs_fft_perform_forw().\n"));
  STOP(TIMING_FORWARD);
  
  /********************************************/
  /* POTENTIAL COMPUTATION */
  /********************************************/
  if (d->require_total_energy || _potentials != NULL) {
    /* apply energy optimized influence function */
    START(TIMING_INFLUENCE)
    ifcs_p3m_apply_energy_influence_function(d);
    /* result is in d->ks_grid */
    STOP(TIMING_INFLUENCE)
            
//for (fcs_int i = d->local_grid.size-1; i >= 0; i--){    
//     printf("ksgrid %f at i=%d\n",d->ks_grid[2*i+1],i);
//  }
    /* compute total energy, but not potentials */
    if (d->require_total_energy && potentials == NULL) {
      START(TIMING_POTENTIALS)
      d->total_energy = ifcs_p3m_compute_total_energy(d);
      STOP(TIMING_POTENTIALS)
    }

    if (_potentials != NULL) {
      /* backtransform the grid */
      P3M_DEBUG(printf( "  calling ifcs_fft_perform_back (potentials)...\n"));
      START(TIMING_BACK)
      ifcs_fft_perform_back(&d->fft, &d->comm, d->ks_grid);
      STOP(TIMING_BACK)
      P3M_DEBUG(printf( "  returned from ifcs_fft_perform_back.\n"));

      /** First (unshifted) run */
      START(TIMING_SPREAD)
      for (fcs_int i=0; i<d->local_grid.size; i++) {
	d->fft.data_buf[i] = d->ks_grid[2*i];
      } 

      ifcs_p3m_spread_grid(d, d->fft.data_buf);

      STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS)

      ifcs_p3m_assign_potentials(d, d->fft.data_buf,
                                 num_real_particles, positions, 
                                 charges, 0, potentials);

      STOPSTART(TIMING_POTENTIALS, TIMING_SPREAD)

      /** Second (shifted) run */
      for (fcs_int i=0; i<d->local_grid.size; i++) {
        d->fft.data_buf[i] = d->ks_grid[2*i+1];
      }
      ifcs_p3m_spread_grid(d, d->fft.data_buf);

      STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS)

      ifcs_p3m_assign_potentials(d, d->fft.data_buf,
                                 num_real_particles, positions,
                                 charges, 1, potentials);
//printf("potentials far:\n");print_array(potentials, num_real_particles);
      STOP(TIMING_POTENTIALS)
    }
  }
  
  /********************************************/
  /* FIELD COMPUTATION */
  /********************************************/
  if (_fields != NULL) {
    /* apply force optimized influence function */
    START(TIMING_INFLUENCE);
    ifcs_p3m_apply_force_influence_function(d);
    STOP(TIMING_INFLUENCE);
    
    /* backtransform the grid */
    START(TIMING_BACK);
    P3M_DEBUG(printf( "  calling ifcs_fft_perform_back...\n"));
    ifcs_fft_perform_back(&d->fft, &d->comm, d->ks_grid);
    P3M_DEBUG(printf( "  returned from ifcs_fft_perform_back.\n"));

    STOPSTART(TIMING_BACK, TIMING_SPREAD)
    
    /* First (unshifted) run */
    P3M_INFO(printf("  computing unshifted grid\n"));
    for (fcs_int i=0; i<d->local_grid.size; i++) {
      d->fft.data_buf[i] = d->ks_grid[2*i];
    } 
    
    ifcs_p3m_spread_grid(d, d->fft.data_buf);

    STOPSTART(TIMING_SPREAD, TIMING_FIELDS)

    ifcs_p3m_assign_fields_ad(d, d->fft.data_buf, num_real_particles, 
                              positions, 0, fields);

    STOPSTART(TIMING_FIELDS, TIMING_SPREAD)
    
    /* Second (shifted) run */
    P3M_INFO(printf("  computing shifted grid\n"));
    for (fcs_int i=0; i<d->local_grid.size; i++) {
      d->fft.data_buf[i] = d->ks_grid[2*i+1];
    }
    
    ifcs_p3m_spread_grid(d, d->fft.data_buf);
    
    STOPSTART(TIMING_SPREAD, TIMING_FIELDS)

    ifcs_p3m_assign_fields_ad(d, d->fft.data_buf, num_real_particles, 
                              positions, 1, fields);
//printf("field in tricl now.\n");
//    print_vector(fields, num_real_particles);
//    fields=cartesian_field(d,fields, num_real_particles);
    printf("field in cart after long range:\n");
    print_vector(fields, num_real_particles);  
    
    STOP(TIMING_FIELDS)
  } 
          
  START(TIMING_COMP)
  /* sort particles back */
  P3M_DEBUG(printf( "  calling fcs_gridsort_sort_backward()...\n"));
  fcs_gridsort_sort_backward(&gridsort,
                             fields, potentials,
                             _fields, _potentials, 1,
                             d->comm.mpicomm);
  P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_backward().\n"));
  printf("positions triclinic after sort back (far): \n"); print_vector(positions, num_real_particles);
   printf("_fields after sort back (far): \n"); print_vector(_fields, num_real_particles); 
  fcs_gridsort_free(&gridsort);
  fcs_gridsort_destroy(&gridsort);

  if (fields != NULL) free(fields);
  if (potentials != NULL) free(potentials);
  
  
  STOP(TIMING_COMP)
              fcs_float* _fields_near = static_cast<fcs_float*>(malloc(3*num_real_particles*sizeof(fcs_float)));
    fcs_float* _potentials_near = static_cast<fcs_float*>(malloc(num_real_particles*sizeof(fcs_float)));
  if (d->near_field_flag) {
      d->cosy_flag=cartesian;
        fcs_int num_real_particles_near;
  fcs_int num_ghost_particles_near;
  fcs_float *positions_near, *ghost_positions_near;
  fcs_float *charges_near, *ghost_charges_near;
  fcs_gridsort_index_t *indices_near, *ghost_indices_near;
  fcs_gridsort_t gridsort_near;
      
       ifcs_p3m_domain_decompose(d, &gridsort_near, 
                            _num_particles, _max_num_particles, _positions, _charges,
                            &num_real_particles_near,
                            &positions_near, &charges_near, &indices_near,
                            &num_ghost_particles_near,
                            &ghost_positions_near, &ghost_charges_near, &ghost_indices_near);
      
      
  /* start near timer */
    START(TIMING_NEAR)  
            
            printf("pos after decomp start of short range:\n");
            print_vector(positions_near, num_real_particles);
        
    /* compute near field */
    fcs_near_t near;
    fcs_float alpha = d->alpha;

    fcs_near_create(&near);
    /*  fcs_near_set_field_potential(&near, ifcs_p3m_compute_near);*/
    fcs_near_set_loop(&near, ifcs_p3m_compute_near_loop);

    fcs_float box_base[3] = {0.0, 0.0, 0.0 };
//    fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
//    fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
//    fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
    
    fcs_near_set_system(&near, box_base, d->box_matrix[0], d->box_matrix[1], d->box_matrix[2], NULL);
    
//    fcs_near_set_particles(&near, num_real_particles, num_real_particles,
//                           positions, charges, indices,
//                           (_fields != NULL) ? fields : NULL, 
//                           (_potentials != NULL) ? potentials : NULL);
     fcs_near_set_particles(&near, num_real_particles, num_real_particles,
                           positions_near, charges_near, indices_near,
                           (_fields != NULL) ? nearfields : NULL, 
                           (_potentials != NULL) ? nearpotentials : NULL);
    fcs_near_set_ghosts(&near, num_ghost_particles_near,
                        ghost_positions_near, ghost_charges_near, ghost_indices_near);
    printf("ghost number near: %d \n", num_ghost_particles_near);
    P3M_DEBUG(printf( "  calling fcs_near_compute()...\n"));
    fcs_near_compute(&near, d->r_cut, &alpha, d->comm.mpicomm);
    P3M_DEBUG(printf( "  returning from fcs_near_compute().\n"));
    
//    P3M_DEBUG(printf("pos after compute near (cart):\n");print_vector(positions,num_real_particles);
    printf("nearfield after compute near (cart):\n");print_vector(nearfields,num_real_particles);
//    printf("nearpots: \n"); print_array(nearpotentials, num_real_particles);)
    fcs_near_destroy(&near);
//    printf("long range fields after compute near: \n");print_vector(fields, num_real_particles);
//    printf("positions after compute near: \n");print_vector(positions, num_real_particles);
    STOP(TIMING_NEAR)
            
                     START(TIMING_COMP)
  /* sort particles back */
  P3M_DEBUG(printf( "  calling fcs_gridsort_sort_backward()...\n"));
    fcs_float* _fields_near = static_cast<fcs_float*>(malloc(3*num_real_particles*sizeof(fcs_float)));
    fcs_float* _potentials_near = static_cast<fcs_float*>(malloc(num_real_particles*sizeof(fcs_float)));
  fcs_gridsort_sort_backward(&gridsort_near,
                             nearfields, nearpotentials,
                             _fields_near, _potentials_near, 1,
                             d->comm.mpicomm);
  P3M_DEBUG(printf( "  returning from fcs_gridsort_sort_backward().\n"));
  printf("positions cart after sort back near: \n"); print_vector(positions, num_real_particles);
  printf("fields cart after sort back near: \n"); print_vector(_fields_near, num_real_particles);
  fcs_gridsort_free(&gridsort);
  fcs_gridsort_destroy(&gridsort);

  
  
  STOP(TIMING_COMP)  

if (nearfields != NULL) free(nearfields);
  if (nearpotentials != NULL) free(nearpotentials);
            
            
  }//end of near field
          
     fcs_int another_counter;
for(another_counter =0 ; another_counter<<num_real_particles; another_counter++){
    potentials[another_counter]+=nearpotentials[another_counter];
    fields[3*another_counter]+=nearfields[3*another_counter];
    fields[3*another_counter+1]+=nearfields[3*another_counter+1];
    fields[3*another_counter+2]+=nearfields[3*another_counter+2];
}
          
 
  /* collect timings from the different nodes */
  if (d->require_timings) 
        collect_print_timings(d);
  d->cosy_flag==cartesian;
  P3M_INFO(printf( "ifcs_p3m_run() finished.\n"));
}

#elif !defined(P3M_INTERLACE) && defined(P3M_IK)
void ifcs_p3m_run(void* rd,
		  fcs_int _num_particles,
		  fcs_int _max_num_particles,
		  fcs_float *_positions, 
		  fcs_float *_charges,
		  fcs_float *_fields,
		  fcs_float *_potentials) {
  /* Here we assume, that the method is already tuned and that all
     parameters are valid */
  P3M_INFO(printf( "ifcs_p3m_run() started...\n"));
  ifcs_p3m_data_struct *d = (ifcs_p3m_data_struct*)rd;

  /* reset all timers */
  if (d->require_timings) {
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
    fields = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*3*num_real_particles));
  if (_potentials != NULL || d->require_total_energy)
    potentials = static_cast<fcs_float*>(malloc(sizeof(fcs_float)*num_real_particles));
  
  STOPSTART(TIMING_DECOMP, TIMING_CA);

  /* charge assignment */
  ifcs_p3m_assign_charges(d, d->rs_grid, num_real_particles, 
                          positions, charges, 0);
  STOPSTART(TIMING_CA, TIMING_GATHER);
  /* gather the ca grid */
  ifcs_p3m_gather_grid(d, d->rs_grid);
  /* now d->rs_grid should contain the local ca grid */
  STOP(TIMING_GATHER);

  /* forward transform */
  START(TIMING_FORWARD);
  P3M_DEBUG(printf( "  calling ifcs_fft_perform_forw()...\n"));
  ifcs_fft_perform_forw(&d->fft, &d->comm, d->rs_grid);
  P3M_DEBUG(printf( "  returned from ifcs_fft_perform_forw().\n"));
  STOP(TIMING_FORWARD);
  
  /********************************************/
  /* POTENTIAL COMPUTATION */
  /********************************************/
  if (d->require_total_energy || _potentials != NULL) {
    /* apply energy optimized influence function */
    START(TIMING_INFLUENCE)
    ifcs_p3m_apply_energy_influence_function(d);
    /* result is in d->ks_grid */
    STOP(TIMING_INFLUENCE)

    /* compute total energy, but not potentials */
    if (d->require_total_energy && potentials == NULL) {
      START(TIMING_POTENTIALS)
      d->total_energy = ifcs_p3m_compute_total_energy(d);
      STOP(TIMING_POTENTIALS)
    }

    if (_potentials != NULL) {
      /* backtransform the grid */
      P3M_DEBUG(printf( "  calling ifcs_fft_perform_back (potentials)...\n"));
      START(TIMING_BACK)
      ifcs_fft_perform_back(&d->fft, &d->comm, d->ks_grid);
      STOP(TIMING_BACK)
      P3M_DEBUG(printf( "  returned from ifcs_fft_perform_back.\n"));

      /* redistribute energy grid */
      START(TIMING_SPREAD)
      ifcs_p3m_spread_grid(d, d->ks_grid);

      STOPSTART(TIMING_SPREAD, TIMING_POTENTIALS)

      /* compute potentials */
      ifcs_p3m_assign_potentials(d, d->ks_grid,
                                 num_real_particles, positions, 
                                 charges, 0,
                                 potentials);

      STOP(TIMING_POTENTIALS)
    }
  }

  /********************************************/
  /* FIELD COMPUTATION */
  /********************************************/
  if (_fields != NULL) {
    /* apply force optimized influence function */
    START(TIMING_INFLUENCE);
    ifcs_p3m_apply_force_influence_function(d);
    STOP(TIMING_INFLUENCE);
    
    /* result is in d->ks_grid */
    for (int dim = 0; dim < 3; dim++) {
      /* differentiate in direction dim */
      /* result is stored in d->rs_grid */
      START(TIMING_FIELDS);
      ifcs_p3m_ik_diff(d, dim);
      
      STOPSTART(TIMING_FIELDS, TIMING_BACK);

      /* backtransform the grid */
      P3M_DEBUG(printf( "  calling ifcs_fft_perform_back (field dim=%d)...\n", dim));
      ifcs_fft_perform_back(&d->fft, &d->comm, d->rs_grid);
      P3M_DEBUG(printf( "  returned from ifcs_fft_perform_back.\n"));
      STOP(TIMING_BACK);
      
      /* redistribute force grid */
      START(TIMING_SPREAD);
      ifcs_p3m_spread_grid(d, d->rs_grid);
      STOPSTART(TIMING_SPREAD, TIMING_FIELDS);
      ifcs_p3m_assign_fields_ik(d, d->rs_grid, dim, num_real_particles, 
                                positions, 0, fields);
      //TODO: transform back for triclinic needs to be put here.
      //
//      int part_no;
//      for (part_no = 0; part_no <  num_real_particles; part_no++){
//          fields[3*part_no]=d->box_matrix[2][0]*fields[3*part_no+2]+d->box_matrix[1][0]*fields[3*part_no+1]+d->box_matrix[0][0]*fields[3*part_no];
//          fields[3*part_no+1]=d->box_matrix[2][1]*fields[3*part_no+2]+d->box_matrix[1][1]*fields[3*part_no+1];
//          fields[3*part_no+2]=d->box_matrix[2][2]*fields[3*part_no+2];
//          
//      }
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
    /*  fcs_near_set_field_potential(&near, ifcs_p3m_compute_near);*/
    fcs_near_set_loop(&near, ifcs_p3m_compute_near_loop);

    fcs_float box_base[3] = {0.0, 0.0, 0.0 };
    
/*
    fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
    fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
    fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
*/
    fcs_near_set_system(&near, box_base, d->box_vector_a, d->box_vector_b, d->box_vector_c, NULL);

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

  /* collect timings from the different nodes */
  if (d->require_timings) 
      collect_print_timings(d);
  

  P3M_INFO(printf( "ifcs_p3m_run() finished.\n"));
}
#endif

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
  fcs_float box_base[3] = {0.0, 0.0, 0.0};
    //  fcs_float box_a[3] = {d->box_l[0], 0.0, 0.0 };
    //  fcs_float box_b[3] = {0.0, d->box_l[1], 0.0 };
    //  fcs_float box_c[3] = {0.0, 0.0, d->box_l[2] };
    fcs_float box_a[3] = {1.0, 0.0, 0.0};
    fcs_float box_b[3] = {0.0, 1.0, 0.0};
    fcs_float box_c[3] = {0.0, 0.0, 1.0 };
  fcs_int num_particles;
  
  printf("positions in decompose:\n");print_vector(_positions, _num_particles);///triclinics
  fcs_gridsort_create(gridsort);
  
  if(d->cosy_flag==triclinic){
      printf("triclinic decomposition. 100 010 001\n");
      fcs_gridsort_set_system(gridsort, box_base, box_a, box_b, box_c, NULL);
  }
  else{
      printf("cartesian decomposition.400 220 222\n");
  fcs_gridsort_set_system(gridsort, box_base, d->box_matrix[0], d->box_matrix[1], d->box_matrix[2], NULL);
  }
  fcs_gridsort_set_particles(gridsort, _num_particles, _max_num_particles, _positions, _charges);

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
ifcs_get_ca_points(ifcs_p3m_data_struct *d, 
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
ifcs_p3m_assign_charges(ifcs_p3m_data_struct* d,
			fcs_float *data,
			fcs_int num_real_particles,
			fcs_float *positions, 
			fcs_float *charges,
			fcs_int shifted) {
  P3M_DEBUG(printf( "  ifcs_p3m_assign_charges() started...\n"));

  const fcs_int q2off = d->local_grid.q_2_off;
  const fcs_int q21off = d->local_grid.q_21_off;

  /* init local charge grid */
  for (fcs_int i=0; i<d->local_grid.size; i++) data[i] = 0.0;

  /* now assign the charges */
  for (fcs_int pid=0; pid < num_real_particles; pid++) {
    const fcs_float q = charges[pid];
    fcs_int linind_grid = ifcs_get_ca_points(d, &positions[pid*3], shifted);

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

  P3M_DEBUG(printf( "  ifcs_p3m_assign_charges() finished...\n"));
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
  fcs_int* d_operator = NULL;
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
  const fcs_int q2off = d->local_grid.q_2_off;
  const fcs_int q21off = d->local_grid.q_21_off;
  const fcs_float prefactor = 1.0 / (d->box_l[0] * d->box_l[1] * d->box_l[2]);
  
  P3M_DEBUG(printf( "  ifcs_p3m_assign_potentials() started...\n"));
  /* Loop over all particles */
  for (fcs_int pid=0; pid < num_real_particles; pid++) {
    fcs_float potential = 0.0;
    fcs_int linind_grid = 
      ifcs_get_ca_points(d, &positions[pid*3], shifted);

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
  const fcs_int q2off = d->local_grid.q_2_off;
  const fcs_int q21off = d->local_grid.q_21_off;
  const fcs_float prefactor = 1.0 / (2.0 * d->box_l[0] * d->box_l[1] * d->box_l[2]);
  const fcs_int dim_rs = (dim+d->ks_pnum) % 3;

  P3M_DEBUG(printf( "  ifcs_p3m_assign_fields() started...\n"));
  /* Loop over all particles */
  for (fcs_int pid=0; pid < num_real_particles; pid++) {
    fcs_float field = 0.0;
    fcs_int linind_grid = 
      ifcs_get_ca_points(d, &positions[3*pid], shifted);

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


  P3M_DEBUG(printf( "  ifcs_p3m_assign_fields() [AD] started...\n"));
  /* Loop over all particles */
  for (fcs_int pid = 0; pid < num_real_particles; pid++) {
    fcs_float field[3] = { 0.0, 0.0, 0.0 };
    fcs_int linind_grid = 
      ifcs_get_ca_points(d, &positions[pid*3], shifted);

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
  P3M_DEBUG(printf( "  ifcs_p3m_assign_fields() finished.\n"));
}
#endif
