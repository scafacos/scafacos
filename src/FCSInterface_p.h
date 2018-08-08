/*
  Copyright (C) 2011, 2012, 2013 Rene Halver, Olaf Lenz, Michael Hofmann

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



#ifndef FCS_INTERFACE_P_INCLUDED
#define FCS_INTERFACE_P_INCLUDED

#include <mpi.h>

#include "FCSDefinitions.h"
#include "FCSResult_p.h"


/**
 * @file FCSInterface_p.h
 * @brief public interface definitions for the main solver-independent
 * functionality of the ScaFaCoS library 
 * @author Rene Halver, Olaf Lenz, Michael Hofmann
 */


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief FCS-object representing an FCS solver
 */
typedef struct _FCS_t *FCS;

#define FCS_NULL  NULL

/**
 * @brief function to initialize an FCS solver
 * @param new_handle pointer to an FCS-object that will represent the FCS solver
 * @param method_name string for selecting the solver method
 * @param communicator MPI communicator to be used for the parallel execution
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_init(FCS *new_handle, const char* method_name, MPI_Comm communicator);

/**
 * @brief function to destroy an FCS solver
 * @param handle FCS-object representing an FCS solver
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_destroy(FCS handle);

/**
 * @brief function to return the numerical identifier of the solver method
 * @param handle FCS-object representing an FCS solver
 * @return numerical identifier of the solver method
 */
fcs_int fcs_get_method(FCS handle);

/**
 * @brief function to return the name of the solver method
 * @param handle FCS-object representing an FCS solver
 * @return name of the solver method
 */
const char *fcs_get_method_name(FCS handle);

/**
 * @brief function to return the MPI communicator used for the parallel execution
 * @param handle FCS-object representing an FCS solver
 * @return MPI communicator used for the parallel execution
 */
MPI_Comm fcs_get_communicator(FCS handle);

/**
 * @brief function to set all obligatory parameters for an FCS solver
 * @param handle FCS-object representing an FCS solver
 * @param near_field_flag whether near-field computations should be performed
 * by the solver method (value 1) or not (value 0)
 * @param box_a first base vector of the system box
 * @param box_b second base vector of the system box
 * @param box_c third base vector of the system box
 * @param box_origin origin vector of the system box
 * @param periodicity periodicity of the system in each dimension (value 0: open, value 1: periodic)
 * @param total_particles total number of particles in the system
 * @return FCSResult-object containing the return state
 */
FCSResult 
fcs_set_common(FCS handle, fcs_int near_field_flag, 
               const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c, 
               const fcs_float *box_origin, 
               const fcs_int *periodicity, fcs_int total_particles);

/**
 * @brief function to set the dimensions of the system
 * @param handle FCS-object representing an FCS solver
 * @param dim dimensions of the system (values 1-3, default: 3)
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_dimensions(FCS handle, fcs_int dim);

/**
 * @brief function to return the dimensions of the system
 * @param handle FCS handle representing an FCS solver object
 * @return dim dimensions of the system
 */
fcs_int fcs_get_dimensions(FCS handle);

/**
 * @brief function to set the near-field flag
 * @param handle FCS-object representing an FCS solver
 * @param near_field_flag whether near-field computations should be performed
 * by the solver method (value 1) or not (value 0)
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_near_field_flag(FCS handle, fcs_int near_field_flag);

/**
 * @brief function to return the near-field flag
 * @param handle FCS-object representing an FCS solver
 * @return value of near-field flag (see ::fcs_set_near_field_flag)
 */
fcs_int fcs_get_near_field_flag(FCS handle);
/**
 * @brief function to set the redistribution of particles before the computations (and back after the computations)
 * @param handle FCS-object representing an FCS solver
 * @param redistribute whether the particles should remain where they are (value 0) or be equally distributed among all processes (value 1)
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_redistribute(FCS handle, fcs_int redistribute);

/**
 * @brief function to return the particle redistribution setting
 * @param handle FCS-object representing an FCS solver
 * @return value of the particle redistribution setting (see ::fcs_set_redistribute)
 */
fcs_int fcs_get_redistribute(FCS handle);

/**
 * @brief function to set the first base vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @param box_a first base vector of the system box
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_box_a(FCS handle, const fcs_float *box_a);

/**
 * @brief function to return the first base vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @return first base vector of the system box
 */
const fcs_float *fcs_get_box_a(FCS handle);

/**
 * @brief function to set the second base vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @param box_b second base vector of the system box
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_box_b(FCS handle, const fcs_float *box_b);

/**
 * @brief function to return the second base vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @return second base vector of the system box
 */
const fcs_float *fcs_get_box_b(FCS handle);

/**
 * @brief function to set the third base vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @param box_c third base vector of the system box
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_box_c(FCS handle, const fcs_float *box_c);

/**
 * @brief function to return the third base vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @return third base vector of the system box
 */
const fcs_float *fcs_get_box_c(FCS handle);

/**
 * @brief function to set the origin vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @param box_origin origin vector of the system box
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_box_origin(FCS handle, const fcs_float *box_origin);

/**
 * @brief function function to return the origin vector of the system box
 * @param handle FCS-object representing an FCS solver
 * @return origin vector of the system box
 */
const fcs_float *fcs_get_box_origin(FCS handle);

/**
 * @brief function to set the periodicity of the system
 * @param handle FCS-object representing an FCS solver
 * @param periodicity periodicity of the system in each dimension (value 0: open, value 1: periodic)
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_periodicity(FCS handle, const fcs_int *periodicity);

/**
 * @brief function to return the periodicity of the system
 * @param handle FCS-object representing an FCS solver
 * @return periodicity of the system (see ::fcs_set_periodicity)
 */
const fcs_int *fcs_get_periodicity(FCS handle);

/**
 * @brief function to set the total number of particles in the system
 * @param handle FCS-object representing an FCS solver
 * @param total_particles total number of particles in the system
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_total_particles(FCS handle, fcs_int total_particles);

/**
 * @brief function to return the total number of particles in the system
 * @param handle FCS-object representing an FCS solver
 * @return total number of particles in the system
 */
fcs_int fcs_get_total_particles(FCS handle);

/**
 * @brief function to set the maximum number of particles that can be stored in the specified local particle data arrays
 * @param handle FCS-object representing an FCS solver
 * @param max_local_particles maximum number of particles that can be stored locally
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_max_local_particles(FCS handle, fcs_int max_local_particles);

/**
 * @brief function to return the maximum number of particles that can be stored in the specified local particle data arrays
 * @param handle FCS-object representing an FCS solver
 * @return maximum number of particles that can be stored locally
 */
fcs_int fcs_get_max_local_particles(FCS handle);

/**
 * @brief function to set the error tolerance of the FCS solver
 * @param handle FCS-object representing an FCS solver
 * @param tolerance_type constant to select the type of the error tolerance
 * value (see FCSDefinitions.h)
 * @param tolerance error tolerance value
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance);

/**
 * @brief function to return the error tolerance of the FCS solver
 * @param handle FCS-object representing an FCS solver
 * @param tolerance_type constant to select the type of the error tolerance
 * value (see FCSDefinitions.h)
 * @param tolerance error tolerance value
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance);

/**
 * @brief function to set a user-defined cutoff radius for the near-field
 * computations of the FCS solver (if supported)
 * @param handle FCS handle representing an FCS solver object
 * @param r_cut cutoff radius for the near-field computations 
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_r_cut(FCS handle, fcs_float r_cut);

/**
 * @brief function to disable a user-defined cutoff radius for the near-field
 * computations of the FCS solver (if supported)
 * @param handle FCS handle representing an FCS solver object
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_unset_r_cut(FCS handle);

/**
 * @brief function to return the user-defined cutoff radius for the near-field
 * computations of the FCS solver (if supported)
 * @param handle FCS handle representing an FCS solver object
 * @param r_cut cutoff radius for the near-field computations
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_r_cut(FCS handle, fcs_float *r_cut);

/**
 * @brief function to set the parameters of the FCS solver based on a parameter string
 * @param handle FCS-object representing an FCS solver
 * @param parameters char* parameter string
 * (format: "<1st parameter name>,<1st comma-separated value(s)>,<2nd parameter name>,<2nd comma-separated value(s)>,...")
 * @param continue_on_errors whether to continue if setting a parameter fails
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_parameters(FCS handle, const char *parameters, fcs_int continue_on_errors);

/**
 * @brief function to print the parameters of an FCS solver to stdout
 * @param handle FCS-object representing an FCS solver
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_print_parameters(FCS handle);

/**
 * @brief function to tune method specific parameters depending on the particles
 * @param handle FCS-object representing an FCS solver
 * @param local_particles local number of particles
 * @param positions positions of the local particles used for the tuning
 * @param charges charges of the local particles used for the tuning
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_tune(FCS handle, fcs_int local_particles,
  fcs_float *positions, fcs_float *charges);

/**
 * @brief function to run the solver method
 * @param handle FCS-object representing an FCS solver
 * @param local_particles local number of particles
 * @param positions positions of the local particles
 * @param charges charges of the local particles
 * @param field calculated field values of the local particles
 * @param potentials calculated potential values of the local particles
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_run(FCS handle, fcs_int local_particles,
  fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials);

/**
 * @brief function to compute the correction to the field and total energy when
 * periodic boundary conditions with a finite dielectric constant of
 * the surrounding medium epsilon are used
 * @param handle FCS-object representing an FCS solver
 * @param local_particles local number of particles
 * @param positions positions of the local particles. These positions should NOT be
 * folded into the system box.
 * @param charges charges of the local particles
 * @param epsilon value of the dielectric constant of the surrounding
 * medium. For metallic boundary conditions, use -1.0 (or simply do
 * not use this function!).
 * @param[out] field_correction required correction to the field values caused
 * by the dipole term.
 * @param[out] energy_correction required correction to the total energy caused
 * by the dipole term.
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_compute_dipole_correction(FCS handle, fcs_int local_particles,
  fcs_float* positions, fcs_float *charges, fcs_float epsilon,
  fcs_float *field_correction, fcs_float *energy_correction);

/**
 * @brief function to return whether the solver method supports the delegation of
 * near-field computations to an external application
 * @param handle FCS-object representing an FCS solver
 * @param near_field_delegation whether the delegation of near field computation is supported
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_near_field_delegation(FCS handle, fcs_int *near_field_delegation);

/**
 * @brief function to compute the near-field components of the potential and the
 * field for the solver method
 * @param handle FCS-object representing an FCS solver
 * @param dist distance between interacting particles
 * @param potential near-field component of the potential
 * @param field near-field component of the field
 *
 * Note: The computed near-field components of the potential and the field do
 * NOT include the charge values of the charges of the interacting particles.
 *
 * Note: For performance reasons, this function does not perform any error checking.
 * It is the responsibility of the user of this function to use it only (1) with
 * a valid FCS-object, (2) when the FCS solver supports the delegation of
 * near-field computations, and (3) all necessary parameters of the FCS solver are
 * set (e.g., with fcs_set_common, fcs_set_... setters, or after ::fcs_run).
 * Furthermore, it is the responsibility of the user of this function to use only
 * dist values between 0.0 and the user-defined cutoff radius of the near-field
 * computations (see ::fcs_set_r_cut).
 */
FCSResult fcs_compute_near(FCS handle, fcs_float dist, fcs_float *potential, fcs_float *field);

/**
 * @brief function to compute the near-field component of the potential for the solver method
 * @param handle FCS-object representing an FCS solver
 * @param dist distance between interacting particles
 * @param potential near-field component of the potential
 * @return FCSResult-object containing the return state
 *
 * Note: See ::fcs_compute_near for details.
 */
FCSResult fcs_compute_near_potential(FCS handle, fcs_float dist, fcs_float *potential);

/**
 * @brief function to compute the near-field component of the field for the solver method
 * @param handle FCS-object representing an FCS solver
 * @param dist distance between interacting particles
 * @param field near-field component of the field
 * @return FCSResult-object containing the return state
 *
 * Note: See ::fcs_compute_near for details.
 */
FCSResult fcs_compute_near_field(FCS handle, fcs_float dist, fcs_float *field);

/**
 * @brief function to set whether the virial should be computed
 * @param handle FCS-object representing an FCS solver
 * @param compute_virial whether the virial should be computed
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_compute_virial(FCS handle, fcs_int compute_virial);

/**
 * @brief function to return whether the virial should be computed
 * @param handle FCS-object representing an FCS solver
 * @param compute_virial whether the virial should be computed
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_compute_virial(FCS handle, fcs_int *compute_virial);

/**
 * @brief function to return the comuputed virial
 * @param handle FCS-object representing an FCS solver
 * @param virial computed virial
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_virial(FCS handle, fcs_float *virial);

/**
 * @brief function to set the maximum distance the particles have moved since the call of ::fcs_run
 * @param handle FCS-object representing an FCS solver
 * @param max_particle_move maximum distance the particles have moved 
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_max_particle_move(FCS handle, fcs_float max_particle_move);

/**
 * @brief function to set whether resort support is requested (default is no resort support)
 *   if resort support is requested (and supported by the solver) then the solver tries to retain its sorted particle data order,
 *   i.e., the distribution and order of the given particles is changed by fcs_run such that the position, charge, field, and potential values correspond to the sorted particle order of the solver,
 *   this can only be performed successfully if the local sizes of the particle arrays (specified by the local_max_particles of fcs_run) on all processes are large enough to store the sorted particle order,
 *   after performing fcs_run, function fcs_get_resort_availability can be used to determine whether resort support is available or not,
 *   if resort support is available then fcs_get_resort_particles returns the new local number of particles and fcs_resort_[ints,floats,bytes] can be used to bring additional particle data into the new sorted order
 * @param handle FCS-object representing an FCS solver
 * @param resort whether resort support is requested
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_set_resort(FCS handle, fcs_int resort);

/**
 * @brief function to return whether resort support is requested
 * @param handle FCS-object representing an FCS solver
 * @param resort whether resort support is requested
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_resort(FCS handle, fcs_int *resort);

/**
 * @brief function to return whether resort support is available after executing ::fcs_run
 * @param handle FCS-object representing an FCS solver
 * @param availability whether resort support is available
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_resort_availability(FCS handle, fcs_int *availability);

/**
 * @brief function to return the new local number of particles
 * @param handle FCS-object representing an FCS solver
 * @param resort_particles new local number of particles
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_get_resort_particles(FCS handle, fcs_int *resort_particles);

/**
 * @brief function to sort additional integer particle data into the new sorted particle order
 * @param handle FCS-object representing an FCS solver
 * @param src array of integer values in unsorted (original) order
 * @param dst array to store the sorted integer values
 * @param n number of integer values for each particle
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n);

/**
 * @brief function to sort additional float particle data into the new sorted particle order
 * @param handle FCS-object representing an FCS solver
 * @param src array of float values in unsorted (original) order
 * @param dst array to store the sorted float values
 * @param n fcs_int number of float values for each particle
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n);

/**
 * @brief function to sort additional byte particle data into the new sorted particle order
 * @param handle FCS-object representing an FCS solver
 * @param src array of byte values in unsorted (original) order
 * @param dst array to store the sorted byte values
 * @param n fcs_int number of byte values for each particle
 * @return FCSResult-object containing the return state
 */
FCSResult fcs_resort_bytes(FCS handle, void *src, void *dst, fcs_int n);


#ifdef FCS_ENABLE_DEPRECATED
#define fcs_set_dimension    fcs_set_dimensions
#define fcs_get_dimension    fcs_get_dimensions
#define fcs_set_offset       fcs_set_box_origin
#define fcs_get_offset       fcs_get_box_origin
#define fcs_printHandle      fcs_print_parameters
#define fcs_method_has_near  fcs_delegate_near_field
#define fcs_parser           fcs_set_parameters
#define fcs_require_virial   fcs_set_compute_virial
#endif


#ifdef __cplusplus
}
#endif


#endif
