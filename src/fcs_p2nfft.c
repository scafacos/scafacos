/*
  Copyright (C) 2011-2013 Michael Pippig
  Copyright (C) 2011-2012 Rene Halver
  Copyright (C) 2011 Sebastian Banert

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



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fcs_p2nfft.h"
#include "../lib/p2nfft/p2nfft.h"

static FCSResult ifcs_p2nfft_check(
    FCS handle, const char* fnc_name);

/* define the checks that are executed at the beginning of each wrapper  */
#define FCS_P2NFFT_INTERFACE_CHECK(FUNCNAME) \
  char* fnc_name =  #FUNCNAME; \
  FCSResult result = ifcs_p2nfft_check(handle, fnc_name); \
  if (result != NULL) return result;

/* Enable definition of multiple wrappers which point to the same internal function.
 * We use this trick to get an P3M-compliant interface but also remain have our own nomenclature. */
#undef FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS
#define FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS 1

/* clear macro definitions of headers (redefine as wrappers later) */
#undef FCS_P2NFFT_INTERFACE_WRAPPER_0
#undef FCS_P2NFFT_INTERFACE_WRAPPER_1
#undef FCS_P2NFFT_INTERFACE_WRAPPER_2
#undef FCS_P2NFFT_INTERFACE_WRAPPER_3

#undef FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS
#define FCS_P2NFFT_INTERFACE_WITH_REDIRECTIONS 1

/* define macros for definition of wrappers with redirection:
 * User sees the function name fcs_p2nfft_ ## NAME
 * that is a redirection to ifcs_p2nnft_ ## INAME */
#define FCS_P2NFFT_INTERFACE_WRAPPER_0(NAME, INAME) \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle) { \
    FCS_P2NFFT_INTERFACE_CHECK(fcs_p2nfft_ ## NAME) \
    return ifcs_p2nfft_ ## INAME(handle->method_context, fnc_name); }
#define FCS_P2NFFT_INTERFACE_WRAPPER_1(NAME, INAME, TYPE1, ARG1) \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1) { \
    FCS_P2NFFT_INTERFACE_CHECK(fcs_p2nfft_ ## NAME) \
    return ifcs_p2nfft_ ## INAME(handle->method_context, fnc_name, ARG1); }
#define FCS_P2NFFT_INTERFACE_WRAPPER_2(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2)   \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1, TYPE2 ARG2) { \
    FCS_P2NFFT_INTERFACE_CHECK(fcs_p2nfft_ ## NAME) \
    return ifcs_p2nfft_ ## INAME(handle->method_context, fnc_name, ARG1, ARG2); }
#define FCS_P2NFFT_INTERFACE_WRAPPER_3(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3)   \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3) { \
    FCS_P2NFFT_INTERFACE_CHECK(fcs_p2nfft_ ## NAME) \
    return ifcs_p2nfft_ ## INAME(handle->method_context, fnc_name, ARG1, ARG2, ARG3); }

/************************************************************
 *     Getter and Setter for P2NFFT, PNFFT, PFFT Parameters 
 ************************************************************/
/* call macros that create the wrappers */
#include "fcs_p2nfft_wrappers.h"


/************************************************************
 *     P2NFFT Interface 
 ************************************************************/
/* initialization function for basic p2nfft parameters */
extern FCSResult fcs_p2nfft_init(
    FCS handle
    )
{
  char* fnc_name =  "fcs_p2nfft_init";
  FCSResult result;

  result = ifcs_p2nfft_check(handle, fnc_name);
  if (result != NULL)
    return result;

  ifcs_p2nfft_init(&(handle->method_context), handle->communicator);

  handle->set_max_particle_move = fcs_p2nfft_set_max_particle_move;
  handle->set_resort = fcs_p2nfft_set_resort;
  handle->get_resort = fcs_p2nfft_get_resort;
  handle->get_resort_availability = fcs_p2nfft_get_resort_availability;
  handle->get_resort_particles = fcs_p2nfft_get_resort_particles;
  handle->resort_ints = fcs_p2nfft_resort_ints;
  handle->resort_floats = fcs_p2nfft_resort_floats;
  handle->resort_bytes = fcs_p2nfft_resort_bytes;

  return NULL;
}

static fcs_int periodic_dims(
    fcs_int *periodicity
    )
{
  fcs_int count = 0;
  for(fcs_int t=0; t<3; t++)
    if(periodicity[t])
      count++;

  return count;
}

static fcs_int nonperiodic_box_lengths_are_equal(
    fcs_float a, fcs_float b, fcs_float c, fcs_int *periodicity
    )
{
  if(!periodicity[0] && !periodicity[1])
    if(!fcs_float_is_equal(a,b))
      return 0;

  if(!periodicity[0] && !periodicity[2])
    if(!fcs_float_is_equal(a,c))
      return 0;

  if(!periodicity[1] && !periodicity[2])
    if(!fcs_float_is_equal(b,c))
      return 0;

  return 1;
}


/* internal p2nfft-specific tuning function */
extern FCSResult fcs_p2nfft_tune(
    FCS handle, fcs_int local_particles, fcs_int local_max_particles,
    fcs_float *positions, fcs_float *charges
    )
{
  char* fnc_name = "fcs_p2nfft_tune";
  FCSResult result;
  fcs_float box_l[3];

  result = ifcs_p2nfft_check(handle, fnc_name);
  if (result != NULL)
    return result;

  /* Check for periodicity */
  fcs_int *periodicity = fcs_get_periodicity(handle);

  /* Check for correct box parameters */
  fcs_float *a = fcs_get_box_a(handle);
  fcs_float *b = fcs_get_box_b(handle);
  fcs_float *c = fcs_get_box_c(handle);
  if (!fcs_uses_principal_axes(a, b, c))
    return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
        "The p2nfft method needs a rectangular box with box vectors parallel to the principal axes.");
  if(periodic_dims(periodicity) == 1)
    if(!nonperiodic_box_lengths_are_equal(a[0], b[1], c[2], periodicity))
      return fcsResult_create(FCS_LOGICAL_ERROR, fnc_name,
          "The p2nfft method currently depends on equal nonperiodic box lengths with 1d-periodic boundary conditions.");
   
  /* Get box size */ 
  box_l[0] = fcs_norm(a);
  box_l[1] = fcs_norm(b);
  box_l[2] = fcs_norm(c);

  /* Call the p2nfft solver's tuning routine */
  result = ifcs_p2nfft_tune(handle->method_context, periodicity,
      local_particles, positions, charges, box_l,
      fcs_get_near_field_flag(handle));

  return result;
}

/* internal p2nfft-specific run function */
extern FCSResult fcs_p2nfft_run(
    FCS handle, fcs_int local_particles, fcs_int local_max_particles,
    fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials
    )
{
  char* fnc_name =  "fcs_p2nfft_run";
  FCSResult result;

  result = ifcs_p2nfft_check(handle, fnc_name);
  if (result != NULL)
    return result;

  result = fcs_p2nfft_tune(handle, local_particles, local_max_particles, positions, charges);
  if (result != NULL)
    return result;
  
  result = ifcs_p2nfft_run(handle->method_context, local_particles, local_max_particles,
      positions, charges, potentials, field);

  return result;
}

/* clean-up function for p2nfft */
extern FCSResult fcs_p2nfft_destroy(
    FCS handle
    )
{
  ifcs_p2nfft_destroy(handle->method_context);
  fcs_set_method_context(handle, NULL);
  return NULL;
}


/*********************************************************************
 * method to check if p2nfft parameters are entered into checked FCS 
 *********************************************************************/
extern FCSResult fcs_p2nfft_check(
    FCS handle
    )
{
  char* fnc_name = "fcs_p2nfft_check";
  return ifcs_p2nfft_check(handle, fnc_name);
}

static FCSResult ifcs_p2nfft_check(
    FCS handle, const char* fnc_name
    )
{
  if (handle == NULL) 
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "Supplied handle must not be a null pointer.");

  if (fcs_get_method(handle) != FCS_P2NFFT)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. Please choose \"p2nfft\".");

  return NULL;
}

/************************************************************
 *     Nearfield computation 
 ************************************************************/

fcs_float fcs_p2nfft_compute_near_potential(
    FCS handle, fcs_float dist
    )
{
  /* Since this function is not fast at all,
   * we also have the time to add some extra safety in the case dist==0.0 */
  if(fcs_float_is_zero(dist))
    return ifcs_p2nfft_compute_self_potential(handle->method_context);
  else
    return ifcs_p2nfft_compute_near_potential(handle->method_context, dist);
}

fcs_float fcs_p2nfft_compute_near_field(
    FCS handle, fcs_float dist
    )
{
  /* Since this function is not fast at all,
   * we also have the time to add some extra safety in the case dist==0.0 */
  if(fcs_float_is_zero(dist))
    return 0.0;
  else
    return ifcs_p2nfft_compute_near_field(handle->method_context, dist);
}

void fcs_p2nfft_compute_near(
    FCS handle, fcs_float dist,
    fcs_float *potential, fcs_float *field
    )
{
  ifcs_p2nfft_compute_near(handle->method_context, dist, field, potential);
}

/************************************************************
 *     Resort support
 ************************************************************/

FCSResult fcs_p2nfft_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  ifcs_p2nfft_set_max_particle_move(handle->method_context, max_particle_move);

  return NULL;
}


FCSResult fcs_p2nfft_set_resort(FCS handle, fcs_int resort)
{
  ifcs_p2nfft_set_resort(handle->method_context, resort);

  return NULL;
}


FCSResult fcs_p2nfft_get_resort(FCS handle, fcs_int *resort)
{
  ifcs_p2nfft_get_resort(handle->method_context, resort);

  return NULL;
}


FCSResult fcs_p2nfft_get_resort_availability(FCS handle, fcs_int *availability)
{
  ifcs_p2nfft_get_resort_availability(handle->method_context, availability);

  return NULL;
}


FCSResult fcs_p2nfft_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  ifcs_p2nfft_get_resort_particles(handle->method_context, resort_particles);

  return NULL;
}


FCSResult fcs_p2nfft_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_resort_ints(handle->method_context, src, dst, n, comm);
  
  return NULL;
}


FCSResult fcs_p2nfft_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_resort_floats(handle->method_context, src, dst, n, comm);
  
  return NULL;
}


FCSResult fcs_p2nfft_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_resort_bytes(handle->method_context, src, dst, n, comm);
  
  return NULL;
}
