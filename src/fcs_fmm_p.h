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



#ifndef FCS_FMM_P_INCLUDED
#define FCS_FMM_P_INCLUDED

#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file fcs_fmm_p.h
 * @brief file containing all fmm specific functions (public version)
 * @author Rene Halver
 */

typedef struct fcs_fmm_parameters_t *fcs_fmm_parameters;

#define FCS_FMM_COULOMB 64
#define FCS_FMM_CUSP 65
#define FCS_FMM_NO_DIPOLE_CORRECTION -1
#define FCS_FMM_STANDARD_DIPOLE_CORRECTION 0


#define FCS_FMM_ACTIVE_DIPOLE_CORRECTION 1
#define FCS_FMM_STANDARD_ERROR 0
#define FCS_FMM_CUSTOM_ABSOLUTE 1
#define FCS_FMM_CUSTOM_RELATIVE 2

#define FCS_FMM_INHOMOGENOUS_SYSTEM 1LL
#define FCS_FMM_HOMOGENOUS_SYSTEM 0LL

/**
 * @brief function to set the optional fmm absrel parameter
 * @param handle FCS-object that is modified
 * @param choice choice for fmm which error boundary should be used:
 *        FCS_FMM_STANDARD_ERROR  - relative error with 10^(-3)
 *        FCS_FMM_CUSTOM_RELATIVE - relative error with deltaE
 *        FCS_FMM_CUSTOM_ABSOLUTE - absolute error with deltaE
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_absrel(FCS handle, fcs_int choice);

/**
 * @brief function to get the optional fmm absrel parameter
 * @param handle FCS-object that contains the parameter
 * @param absrel pointer to fcs_int
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_absrel(FCS handle, fcs_int *absrel);

/**
 * @brief function to set the optional fmm deltaE parameter, which gives the energy tolerance for the FMM
 * @param handle FCS-object that is modified
 * @param tolerance_value the error boundary to be used by fmm (only applicable
 *        if absrel is either 1 or 2)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_tolerance_energy(FCS handle, fcs_float tolerance_value);

/**
 * @brief function to get the optional fmm energy tolerance (deltaE) parameter
 * @param handle FCS-object that contains the parameter
 * @param tolerance_value pointer to a fcs_float variable (only applicable if absrel is either 1 or 2)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_tolerance_energy(FCS handle, fcs_float *tolerance_value);

/**
 * @brief function to set the optional fmm energy tolerance (deltaE) parameter
 * @param handle FCS-object that is modified
 * @param dipole_correction chooses which form of dipole correction fmm should use:
 *        FCS_FMM_NO_DIPOLE_CORRECTION        no dipole correction
 *        FCS_FMM_STANDARD_DIPOLE_CORRECTION  standard dipole correction
 *        FCS_FMM_ACTIVE_DIPOLE_CORRECTION    dipole correction activated
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_dipole_correction(FCS handle, fcs_int dipole_correction);

/**
 * @brief function to get the optional fmm dipole correction parameter
 * @param handle FCS-object that contains the parameter
 * @param dipole_correction pointer to fcs_int variable where the
 * function returns the dipole correction
 * @return the parameter dipole correction
 */
extern FCSResult fcs_fmm_get_dipole_correction(FCS handle, fcs_int *dipole_correction);

/**
 * @brief function to set the fmm parameter parameter
 * @param handle FCS-object that is modified
 * @param coulomb fcs_int indicates if coulomb potential is used
 *        FCS_FMM_COULOMB if coulomb potential is not used
 *        FCS_FMM_CUSP    if cusp potential is used
 *        (exclusive to FFM_cusp!)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_potential(FCS handle, fcs_int coulomb);

/**
 * @brief function to get the fmm coulomb parameter
 * @param handle FCS-object that contains the parameter
 * @param coulomb pointer to fcs_int variable where the function
 * returns the potential parameter
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_potential(FCS handle, fcs_int *coulomb);

/**
 * @brief function to set the fmm internal tuning parameter
 * @param handle FCS-object that is modified
 * @param system long long defining the kind of system to tune for
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_internal_tuning(FCS handle, long long system);

/**
 * @brief function to get the fmm internal tuning parameter
 * @param handle FCS-object that contains the parameter 
 * @param system pointer to long long variable stating 
 * which kind of system is used
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_internal_tuning(FCS handle, long long *system);

/**
 * @brief function to set the fmm cusp radius parameter
 * @param handle FCS-object that is modified
 * @param radius fcs_float the radius for the cusp potential
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_cusp_radius(FCS handle, fcs_float radius);

/**
 * @brief function to get the fmm cusp radius parameter
 * @param handle FCS-object that contains the parameter 
 * @param coulomb pointer to fcs_float variable where the function
 * returns the cusp radius
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_cusp_radius(FCS handle, fcs_float *cusp_radius);

/**
 * @brief function to set the optional fmm absrel parameter
 * @param handle FCS-object that is modified
 * @param depth fcs_int containing the maximum tree depth for the FMM (0 to 19)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_maxdepth(FCS handle, long long depth);

/**
 * @brief function to get the optional fmm absrel parameter
 * @param handle FCS-object that contains the parameter
 * @param depth  fcs_int containing the set maximum tree depth 
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_maxdepth(FCS handle, long long *depth);

/**
 * @brief function to set the optional fmm absrel parameter
 * @param handle FCS-object that is modified
 * @param limit fcs_int the limit for unrolled functions within the FMM (0 to 50)
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_unroll_limit(FCS handle, long long limit);

/**
 * @brief function to get the optional fmm absrel parameter
 * @param handle FCS-object that contains the parameter
 * @param depth  fcs_int containing the limit for the unrolled functions
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_unroll_limit(FCS handle, long long *limit);

/**
 * @brief function to set the optional fmm absrel parameter
 * @param handle FCS-object that is modified
 * @param depth fcs_int activates load balancing routines (0 (deactivated) or 1 (activated))
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_balanceload(FCS handle, long long load);

/**
 * @brief function to get the optional fmm absrel parameter
 * @param handle FCS-object that contains the parameter
 * @param depth  fcs_int containing the load balancing status
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_get_balanceload(FCS handle, long long *load);

/**
 * @brief function to set the optional fmm loadvector parameter
 * @param handle FCS-object that contains the parameter
 * @param define_loadvector  fcs_int containing the load balancing initialization status
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_set_define_loadvector(FCS handle, long long  define_loadvector);
extern FCSResult fcs_fmm_get_define_loadvector(FCS handle, long long *define_loadvector);

/*
 * @brief combined setter function for all fmm related parameters
 * @param handle FCS-object that is modified
 * @param choice fcs_int choice for fmm which error boundary should be used:
 *        0 - relative error with 10^(-3)
 *        1 - relative error with energy tolerance (deltaE)
 *        2 - absolute error with energy tolerance (deltaE)
 * @param tolerance_value [deltaE] fcs_float error boundary to be used by fmm (only applicable
 *        if absrel is either 1 or 2)
 * @param dipole_correction fcs_int chooses which form of dipole correction fmm should use:
 *        -1 = no dipole correction
 *         0 = standard dipole correction
 *         1 = dipole correction activated
 * @param potential fcs_int chooses which potential to use
 *         1 = coulomb potential (ignores radius, see below)
 *         2 = cusp potential (requires radius)
 * @param radius fcs_float sets the radius for the cusp potential, if chosen
 * @return FCSResult-object containing the return state
 */
extern FCSResult fcs_fmm_setup(FCS handle, fcs_int absrel, fcs_float tolerance_value, fcs_int dipole_correction, long long system, long long maxdepth, long long unroll_limit, long long load/*, fcs_int potential, fcs_float radius*/);

void fcs_fmm_setup_f(void *handle, fcs_int absrel, fcs_float tolerance_value, fcs_int dipole_correction, fcs_int *return_value);

#ifdef __cplusplus
}
#endif

#endif
