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



#ifndef FCS_INTERFACE_INCLUDED
#define FCS_INTERFACE_INCLUDED

#include <mpi.h>

#include "FCSInterface_p.h"
#include "FCSResult.h"

#ifdef FCS_ENABLE_DIRECT
#include "fcs_direct.h"
#endif
#ifdef FCS_ENABLE_EWALD
#include "fcs_ewald.h"
#endif
#ifdef FCS_ENABLE_FMM
#include "fcs_fmm.h"
#endif
#ifdef FCS_ENABLE_MEMD
#include "fcs_memd.h"
#endif
#ifdef FCS_ENABLE_MMM1D
#include "fcs_mmm1d.h"
#endif
#ifdef FCS_ENABLE_MMM2D
#include "fcs_mmm2d.h"
#endif
#ifdef FCS_ENABLE_P2NFFT
#include "fcs_p2nfft.h"
#endif
#ifdef FCS_ENABLE_P3M
#include "fcs_p3m.h"
#endif
#ifdef FCS_ENABLE_PEPC
#include "fcs_pepc.h"
#endif
#ifdef FCS_ENABLE_PP3MG
#include "fcs_pp3mg.h"
#endif
#ifdef FCS_ENABLE_VMG
#include "fcs_vmg.h"
#endif
#ifdef FCS_ENABLE_WOLF
#include "fcs_wolf.h"
#endif


/**
 * @file FCSInterface.h
 * @brief supplying the C-interface routine for the ScaFaCoS library
 * @author Rene Halver, Olaf Lenz
 */

/* definition of FCS_t structure */
typedef struct FCS_t
{
  /* list of obligatory parameters for all runs (set in fcs_init) */

  /* the chosen method */
  fcs_int method;
  /* communicator to be used in methods */
  MPI_Comm communicator;
  /* dimension of the system */
  fcs_int dimension;
  /* flag marking near field calculation
     0 = done by calling routine
     1 = done by library routine*/
  fcs_int near_field_flag;
  /* the three vectors spanning the system box */
  fcs_float box_a[3];
  fcs_float box_b[3];
  fcs_float box_c[3];
  /* offset of each particle [for boxes from -L/2:L/2] */
  fcs_float offset[3];
  /* total number of particles */
  fcs_int total_particles;
  /* total number of charges */
  fcs_int total_charges;
  /* periodicity, fcs_int[3] to set in which direction the */
  /* system has to be periodic, as a result the sum over   */
  /* the array gives the dimension of periodicity */
  /* 0 = system is not periodic in given direction */
  /* 1 = system is periodic in given direction */
  /* [0] -> x, [1] -> y, [2] -> z */
  fcs_int periodicity[3];
  
  /* structures containing the method-specific parameters */
#ifdef FCS_ENABLE_DIRECT
  fcs_direct_parameters direct_param;
#endif
#ifdef FCS_ENABLE_FMM
  fcs_fmm_parameters fmm_param;
#endif
#ifdef FCS_ENABLE_MEMD
  fcs_memd_parameters memd_param;
#endif
#ifdef FCS_ENABLE_P2NFFT
  fcs_p2nfft_parameters p2nfft_param;
#endif
#ifdef FCS_ENABLE_PEPC
  fcs_pepc_parameters pepc_param;
#endif
#ifdef FCS_ENABLE_PP3MG
  fcs_pp3mg_parameters pp3mg_param;
#endif
#ifdef FCS_ENABLE_VMG
  fcs_vmg_parameters vmg_param;
#endif
#ifdef FCS_ENABLE_WOLF
  fcs_wolf_parameters wolf_param;
#endif

  /* current instance of the method */
  void *method_context;

  fcs_int values_changed;

  FCSResult (*set_max_particle_move)(FCS handle, fcs_float max_particle_move);
  FCSResult (*set_resort)(FCS handle, fcs_int resort);
  FCSResult (*get_resort)(FCS handle, fcs_int *resort);
  FCSResult (*get_resort_availability)(FCS handle, fcs_int *availability);
  FCSResult (*get_resort_particles)(FCS handle, fcs_int *resort_particles);
  FCSResult (*resort_ints)(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm);
  FCSResult (*resort_floats)(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm);
  FCSResult (*resort_bytes)(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm);

} FCS_t;

/**
 * @brief getter function for values changed
 * @param handle FCS-object to be changed
 * @return fcs_int indicating if one of the parameters in
 * handle has been changed since the last call to fcs_tune
 */
fcs_int fcs_get_values_changed(FCS handle);

/**
 * @brief setter function for values changed
 * @param handle FCS-object to be changed
 * @param changed indicating if one of the parameters in 
 * handle has been changed since the last call to fcs_tune
 * @return FCSResult-object containing the return state
 * or NULL if successful
 */
FCSResult fcs_set_values_changed(FCS handle, fcs_int changed);


#endif
