/*
  Copyright (C) 2011, 2012, 2013 Rene Halver, Michael Hofmann

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


/**
 * @file FCSInterface.h
 * @brief supplying the C-interface routine for the ScaFaCoS library
 * @author Rene Halver, Olaf Lenz
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


#define FCS_MAX_METHOD_NAME_LENGTH  32


/*
 * @brief data structure that is used for storing the parameters of an FCS solver
 */
typedef struct _FCS_t
{
  /* list of obligatory parameters for all runs (set in fcs_init) */

  /* numerical identifier of the solver method */
  fcs_int method;
  /* name of the solver method */
  char method_name[FCS_MAX_METHOD_NAME_LENGTH];
  /* MPI communicator to be used for the parallel execution */
  MPI_Comm communicator;
  /* dimensions of the system */
  fcs_int dimensions;
  /* whether near-field computations should be performed
     0 = done by calling routine
     1 = done by library routine*/
  fcs_int near_field_flag;
  /* the three base vectors of the system box */
  fcs_float box_a[3];
  fcs_float box_b[3];
  fcs_float box_c[3];
  /* origin vector of the system box */
  fcs_float box_origin[3];
  /* total number of particles in the system */
  fcs_int total_particles;
  /* periodicity of the system in each dimension (value 0: open, value 1: periodic) */
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
 * @brief function to set the method context information
 * @param handle FCS-object representing an FCS solver
 * @param pointer to the method context informations
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_method_context(FCS handle, void *method_context);

/**
 * @brief function to return the method context information
 * @param handle FCS handle representing an FCS solver object
 * @return pointer to the method context informations
 */
void *fcs_get_method_context(FCS handle);

/**
 * @brief function to set whether parameter values of the FCS solver have changed
 * @param handle FCS-object representing an FCS solver
 * @param values_changed whether parameter values of the FCS solver have changed
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_values_changed(FCS handle, fcs_int values_changed);

/**
 * @brief function to return whether parameter values of the FCS solver have changed
 * @param handle FCS-object representing an FCS solver
 * @return whether parameter values of the FCS solver have changed
 */
fcs_int fcs_get_values_changed(FCS handle);

/**
 * @brief Fortran wrapper function ot initialize an FCS solver method
 * @param handle FCS-object representing an FCS solver
 * @param method string for selecting the solver method
 * @param communicator MPI communicator to be used for parallel execution
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_init_f(FCS *handle, const char *method_name, MPI_Fint communicator);


#endif
