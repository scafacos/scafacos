/*
  Copyright (C) 2011-2012 Rene Halver

  This file is part of ScaFaCoS.

  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.

  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/



#ifndef FCS_FMM_INCLUDED
#define FCS_FMM_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "fcs_fmm_p.h"
#include "FCSResult.h"
#include "FCSInterface.h"
#include "fmm/sl_fmm/mpi_fmm_resort.h"


/**
 * @file fcs_fmm.h
 * @brief file containing all fmm specific functions
 * @author Rene Halver
 */


typedef struct fcs_fmm_parameters_t
{
  /* list of parameters used by fmm */
  /* what error type should be used? */
  /* 0 -> relative error 10^(-3) */
  /* 1 -> relative error <fmm_deltaE> */
  /* 2 -> absolute error <fmm_deltaE> */
  fcs_int absrel;
  /* size of chosen error [10^(-1) .. 10^(-14)] */
  fcs_float tolerance_value;
  /* dipole correction */
  fcs_int dipole_correction;
  /* maximum tree depth [0 .. 19] */
  fcs_int maxdepth;
  /* limit for unrolled functions [0 .. 50] */
  fcs_int limit;
  /* status of the load balancing */
  fcs_int balance;
  /* using Coulomb potential */
  fcs_int potential;
  /* radius for cusp potential */
  fcs_float cusp_radius;
  /* load balancing vector */
  fcs_int define_loadvector;
  /* internal fmm tuning for inhomogenous systems */
  long long system;

  /* storage space for the virial */
  fcs_float virial[9];

  /* size and memory pointer for wigner */
  long long wignersize;
  void *wignerptr;
  
  /* resort parameters */
  fcs_float max_particle_move;
  fcs_int resort;
  fcs_fmm_resort_t fmm_resort;

} fcs_fmm_parameters_t;


/**
 * @brief function to check if obligatory fmm information is set
 * @param handle FCS-object to be checked if obligatory
 * information needed by fmm is set (no plausibility check!)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_check(FCS handle, fcs_int local_particles);

/**
 * @brief initialization routine for the basic parameters needed by fmm
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_init(FCS handle );

/**
 * @brief tuning routine for the basic parameters needed by fmm
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions,  fcs_float *charges);

/**
 * @brief run method for fmm
 * @param handle the FCS-object into which the method specific parameters
 * can be entered
 * @param local_particles actual number of particles on process
 * @param local_max_particles size of allocated arrays
 * @param positons fcs_float* list of positions of particles in form
 *        (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)
 * @param charges fcs_float* list of charges
 * @param output FCSOutput* pointer that contains a FCSOutput-object with the
 * results after the run
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
                             fcs_float *positions, fcs_float *charges,
                             fcs_float *field, fcs_float *potentials);

/**
 * @brief a method to calculate fmm near-field potentials with the chosen method (without charges of particles)
 *        e.g.: instead of q1q2/r the function return 1/r
 * @param handle pointer to a FCS on which the FCS with the information
 * are saved
 * @param local_particles the amount of particles on local process
 * @param r distance for which the potential
 * @return fcs_float containing the potential value for the given distance (without charges, s. above)
 */
extern FCSResult fcs_fmm_near_field_potential (FCS handle, fcs_float);

/**
 * @brief clean-up method for fmm
 * @param handle the FCS-object, which contains the parameters
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_destroy(FCS handle);

/**
 * @brief function to activate computation of the virial 
 * @param handle FCS-object that contains the parameter
 * @param flag whether or not to compute the virial in the next call
 * to fcs_run().
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_require_virial(FCS handle, fcs_int flag);

/**
 * @brief function to fetch the virial
 * @param handle FCS-object that contains the parameter
 * @param virial pointer to the array where the virial is returned.
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_virial(FCS handle, fcs_float *virial);


FCSResult fcs_fmm_set_max_particle_move(FCS handle, fcs_float max_particle_move);
FCSResult fcs_fmm_set_resort(FCS handle, fcs_int resort);
FCSResult fcs_fmm_get_resort(FCS handle, fcs_int *resort);
FCSResult fcs_fmm_get_resort_availability(FCS handle, fcs_int *availability);
FCSResult fcs_fmm_get_resort_particles(FCS handle, fcs_int *resort_particles);
FCSResult fcs_fmm_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);
FCSResult fcs_fmm_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);
FCSResult fcs_fmm_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm);


#endif
