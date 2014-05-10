/*
  Copyright (C) 2011,2012 Rene Halver, Olaf Lenz

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

#ifndef _FCS_P3M_P_H
#define _FCS_P3M_P_H

#include "FCSDefinitions.h"
#include "FCSResult_p.h"
#include "FCSInterface_p.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file fcs_p3m_p.h
 * @brief file containing all p3m specific functions
 * @author Rene Halver, Olaf Lenz
 */

/* TODO: test whether to use it or not */
#define FCS_P3M_USE_ERFC_APPROXIMATION 1
  /* #define FCS_P3M_USE_ERFC_APPROXIMATION 0 */

/**
 * @brief The struct that keeps the parameters for the near field
 * component of the method.
 *
 * The struct can be obtained via fcs_p3m_get_near_parameters()
 * and is used in fcs_p3m_compute_near_potential(),
 * fcs_p3m_compute_near_field() and fcs_p3m_compute_near().
 */
/* This struct is defined open so that an MD implementation can use
   the parameters to perform the near field computation in its own
   code. */
typedef struct {fcs_float alpha; fcs_float potentialOffset;} fcs_p3m_near_parameters_t;

/**
 * @brief extracts the parameters required to compute the near-field
 * component of p3m from the method handle
 * @param handle the FCS-object, which contains the parameters
 * @param near_params the fcs_p3m_near_parameters_t object that will be
 * updated with the near-field component parameters of the method.
 * @return FCSResult-object containing the return state
 *
 * The struct can be obtained via fcs_p3m_get_near_parameters()
 * and is used in fcs_p3m_compute_near_potential(),
 * fcs_p3m_compute_near_field() and fcs_p3m_compute_near().
 * Note that the set of near field parameters might change whenever
 * any of the parameters of the method changes. Therefore make sure to
 * update the near field parameters after each parameter change.
 */
FCSResult 
fcs_p3m_get_near_parameters(FCS handle,
			    fcs_p3m_near_parameters_t *near_params);

/**
 * @brief compute the near-field component of the potential of p3m
 * @param params the struct that contains the parameters for the near
 * field computation
 * @param dist the distance of both charges
 * @return the value of the potential. Note that the function will
 * return 1/dist, not q1q2/dist.
 * 
 * Note that it is in the responsibility of the user of this function
 * to make sure that the function will get values of r between 0.0 and
 * the cutoff range. Values outside of this range might result in
 * undefined behavior.
 */
/* This function is defined inline for maximal performance! */
static inline fcs_float 
fcs_p3m_compute_near_potential(fcs_p3m_near_parameters_t params, fcs_float dist) {
  const fcs_float alpha = params.alpha;
  const fcs_float potentialOffset = params.potentialOffset;  
  fcs_float adist = alpha * dist;

#if FCS_P3M_USE_ERFC_APPROXIMATION

  /* approximate \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
     Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
     (9. ed.), chapter 7 */
  fcs_float t = 1.0 / (1.0 + 0.3275911 * adist);
  fcs_float erfc_part_ri = exp(-adist*adist) * 
    (t * (0.254829592 + 
	  t * (-0.284496736 + 
	       t * (1.421413741 + 
		    t * (-1.453152027 + 
			 t * 1.061405429))))) 
    / dist;
  
#else

  /* use erf instead of erfc to fix ICC performance problems */
  fcs_float erfc_part_ri = (1.0 - erf(adist)) / dist; 

#endif

  return erfc_part_ri-potentialOffset;

}

/**
 * @brief compute the near-field component of the field of p3m
 * @param params the struct that contains the parameters for the near
 * field computation
 * @param dist the distance of both charges
 * @return the value of the field. Note that the function will
 * return 1/dist, not q1q2/dist
 * 
 * Note that it is in the responsibility of the user of this function
 * to make sure that the function will get values of r between 0.0 and
 * the cutoff range. Values outside of this range might result in
 * undefined behavior.
 */
static inline fcs_float 
fcs_p3m_compute_near_field(fcs_p3m_near_parameters_t params, fcs_float dist) {
  const fcs_float alpha = params.alpha;
  fcs_float adist = alpha * dist;

#if FCS_P3M_USE_ERFC_APPROXIMATION

  /* approximate \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
     Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
     (9. ed.), chapter 7 */
  fcs_float t = 1.0 / (1.0 + 0.3275911 * adist);
  fcs_float erfc_part_ri = exp(-adist*adist) * 
    (t * (0.254829592 + 
	  t * (-0.284496736 + 
	       t * (1.421413741 + 
		    t * (-1.453152027 + 
			 t * 1.061405429))))) 
    / dist;
  /* 1/sqrt(PI) = 0.56418958354775627928034964498 */
  return -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist)) 
    / dist;
  
#else

  fcs_float erfc_part_ri = (1.0 - erf(adist)) / dist; /* use erf instead of erfc to fix ICC performance problems */
  /* 1/sqrt(PI) = 0.56418958354775627928034964498 */
  return -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist))
    / dist;

#endif

}


/**
 * @brief compute the near-field component of the potential and the field of p3m
 * @param params the struct that contains the parameters for the near
 * field computation
 * @param dist the distance of both charges
 * @param potential pointer to the fcs_float variable where the
 * potential will be written to
 * @param field pointer to the fcs_float variable where the magnitude
 * of the field will be written to
 * @return FCSResult-object containing the return state
 * 
 * Note that it is in the responsibility of the user of this function
 * to make sure that the function will get values of r between 0.0 and
 * the cutoff range. Values outside of this range might result in
 * undefined behavior.
 */
  static inline void
fcs_p3m_compute_near(fcs_p3m_near_parameters_t params, fcs_float dist, 
		     fcs_float *potential, fcs_float *field) {
  const fcs_float alpha = params.alpha;
  const fcs_float potentialOffset = params.potentialOffset;
  fcs_float adist = alpha * dist;

#if FCS_P3M_USE_ERFC_APPROXIMATION

  /* approximate \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
     Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
     (9. ed.), chapter 7 */
  fcs_float t = 1.0 / (1.0 + 0.3275911 * adist);
  fcs_float erfc_part_ri = exp(-adist*adist) * 
    (t * (0.254829592 + 
	  t * (-0.284496736 + 
	       t * (1.421413741 + 
		    t * (-1.453152027 + 
			 t * 1.061405429))))) 
    / dist;

  *potential = erfc_part_ri-potentialOffset;
  *field = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist))
    / dist;
  
#else

  fcs_float erfc_part_ri = (1.0 - erf(adist)) / dist; /* use erf instead of erfc to fix ICC performance problems */
  *potential = erfc_part_ri-potentialOffset;
  *field = -(erfc_part_ri + 2.0*alpha*0.56418958354775627928034964498*exp(-adist*adist)) 
    / dist;

#endif

}

FCSResult fcs_p3m_set_r_cut(FCS handle, fcs_float r_cut);
FCSResult fcs_p3m_set_r_cut_tune(FCS handle);
FCSResult fcs_p3m_get_r_cut(FCS handle, fcs_float *r_cut);

FCSResult fcs_p3m_set_alpha(FCS handle, fcs_float alpha);
FCSResult fcs_p3m_set_alpha_tune(FCS handle);
FCSResult fcs_p3m_get_alpha(FCS handle, fcs_float *alpha);

FCSResult fcs_p3m_set_grid(FCS handle, fcs_int grid);
FCSResult fcs_p3m_set_grid_tune(FCS handle);
FCSResult fcs_p3m_get_grid(FCS handle, fcs_int *grid);

FCSResult fcs_p3m_set_cao(FCS handle, fcs_int cao);
FCSResult fcs_p3m_set_cao_tune(FCS handle);
FCSResult fcs_p3m_get_cao(FCS handle, fcs_int *cao);

FCSResult fcs_p3m_require_total_energy(FCS handle, fcs_int total_energy);
FCSResult fcs_p3m_get_total_energy(FCS handle, fcs_float *total_energy);

FCSResult fcs_p3m_set_tolerance_field(FCS handle, fcs_float tolerance_field);
/* FORTRAN wrapper */
void fcs_p3m_set_tolerance_field_f(void *handle, fcs_float tolerance_field, fcs_int *return_value);
FCSResult fcs_p3m_set_tolerance_field_tune(FCS handle);
FCSResult fcs_p3m_get_tolerance_field(FCS handle, fcs_float *tolerance_field);

FCSResult fcs_p3m_set_potential_shift(FCS handle, fcs_int flag);
FCSResult fcs_p3m_get_potential_shift(FCS handle, fcs_int *flag);

  /*
FCSResult fcs_p3m_set_variant(FCS handle, fcs_int variant);
FCSResult fcs_p3m_get_variant(FCS handle, fcs_int *variant);
  */

#ifdef __cplusplus
}
#endif

#endif
