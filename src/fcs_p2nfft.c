/*
  Copyright (C) 2011-2013 Michael Pippig
  Copyright (C) 2011-2012 Rene Halver
  Copyright (C) 2011 Sebastian Banert
  Copyright (C) 2016 Michael Hofmann

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

#define P2NFFT_CHECK_RETURN_RESULT(_h_, _f_)  do { \
  CHECK_HANDLE_RETURN_RESULT(_h_, _f_); \
  CHECK_METHOD_RETURN_RESULT(_h_, _f_, FCS_METHOD_P2NFFT, "p2nfft"); \
  } while (0)

#define P2NFFT_CHECK_RETURN_VAL(_h_, _f_, _v_)  do { \
  CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_); \
  CHECK_METHOD_RETURN_VAL(_h_, _f_, FCS_METHOD_P2NFFT, "p2nfft", _v_); \
  } while (0)

/* Enable definition of multiple wrappers which point to the same internal function.
 * We use this trick in order to get an P3M-compliant interface but also remain our own nomenclature. */
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
    P2NFFT_CHECK_RETURN_RESULT(handle, __func__); \
    return ifcs_p2nfft_ ## INAME(handle->method_context, __func__); }
#define FCS_P2NFFT_INTERFACE_WRAPPER_1(NAME, INAME, TYPE1, ARG1) \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1) { \
    P2NFFT_CHECK_RETURN_RESULT(handle, __func__); \
    return ifcs_p2nfft_ ## INAME(handle->method_context, __func__, ARG1); }
#define FCS_P2NFFT_INTERFACE_WRAPPER_2(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2)   \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1, TYPE2 ARG2) { \
    P2NFFT_CHECK_RETURN_RESULT(handle, __func__); \
    return ifcs_p2nfft_ ## INAME(handle->method_context, __func__, ARG1, ARG2); }
#define FCS_P2NFFT_INTERFACE_WRAPPER_3(NAME, INAME, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3)   \
  FCSResult fcs_p2nfft_ ## NAME(FCS handle, TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3) { \
    P2NFFT_CHECK_RETURN_RESULT(handle, __func__); \
    return ifcs_p2nfft_ ## INAME(handle->method_context, __func__, ARG1, ARG2, ARG3); }

/************************************************************
 *     Getter and Setter for P2NFFT, PNFFT, PFFT Parameters 
 ************************************************************/
/* call macros that create the wrappers */
#include "fcs_p2nfft_wrappers.h"


/************************************************************
 *     P2NFFT Interface 
 ************************************************************/
/* initialization function for basic p2nfft parameters */
FCSResult fcs_p2nfft_init(
    FCS handle
    )
{
  P2NFFT_CHECK_RETURN_RESULT(handle, __func__);

  handle->shift_positions = 0;

  handle->destroy = fcs_p2nfft_destroy;
  handle->set_r_cut = fcs_p2nfft_set_r_cut;
  handle->unset_r_cut = fcs_p2nfft_set_r_cut_tune;
  handle->get_r_cut = fcs_p2nfft_get_r_cut;
  handle->set_tolerance = fcs_p2nfft_set_tolerance;
  handle->get_tolerance = fcs_p2nfft_get_tolerance;
  handle->set_parameter = fcs_p2nfft_set_parameter;
  handle->print_parameters = fcs_p2nfft_print_parameters;
  handle->tune = fcs_p2nfft_tune;
  handle->run = fcs_p2nfft_run;
  handle->set_compute_virial = fcs_p2nfft_require_virial;
  handle->get_virial = fcs_p2nfft_get_virial;
  handle->set_max_particle_move = fcs_p2nfft_set_max_particle_move;
  handle->set_resort = fcs_p2nfft_set_resort;
  handle->get_resort = fcs_p2nfft_get_resort;
  handle->get_resort_availability = fcs_p2nfft_get_resort_availability;
  handle->get_resort_particles = fcs_p2nfft_get_resort_particles;
  handle->resort_ints = fcs_p2nfft_resort_ints;
  handle->resort_floats = fcs_p2nfft_resort_floats;
  handle->resort_bytes = fcs_p2nfft_resort_bytes;

  ifcs_p2nfft_init(&(handle->method_context), handle->communicator);

  return FCS_RESULT_SUCCESS;
}

static fcs_int periodic_dims(
    const fcs_int *periodicity
    )
{
  fcs_int count = 0;
  for(fcs_int t=0; t<3; t++)
    if(periodicity[t])
      count++;

  return count;
}

static fcs_int nonperiodic_box_lengths_are_equal(
    fcs_float a, fcs_float b, fcs_float c, const fcs_int *periodicity
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

static fcs_int is_principal_axis(const fcs_float *a)
{
  return (1 == (!fcs_float_is_zero(a[0]) +  !fcs_float_is_zero(a[1]) + !fcs_float_is_zero(a[2])) );
}

static fcs_int nonperiodic_axes_are_principal(
    const fcs_float *a, const fcs_float *b, const fcs_float *c, const fcs_int *periodicity
    )
{
  if( !periodicity[0] && !is_principal_axis(a) ) return 0;
  if( !periodicity[1] && !is_principal_axis(b) ) return 0;
  if( !periodicity[2] && !is_principal_axis(c) ) return 0;
  return 1;
}

static fcs_int nonperiodic_axes_are_orthogonal_to_all_other_axes(
    const fcs_float *a, const fcs_float *b, const fcs_float *c, const fcs_int *periodicity
    )
{
  /* We use the fact that nonperiodic axes are principal axes */
  if( !periodicity[0] ) if( !fcs_float_is_zero(b[0]) || !fcs_float_is_zero(c[0]) ) return 0;
  if( !periodicity[1] ) if( !fcs_float_is_zero(a[1]) || !fcs_float_is_zero(c[1]) ) return 0;
  if( !periodicity[2] ) if( !fcs_float_is_zero(a[2]) || !fcs_float_is_zero(b[2]) ) return 0;
  return 1;
}



/* internal p2nfft-specific tuning function */
FCSResult fcs_p2nfft_tune(
    FCS handle, fcs_int local_particles,
    fcs_float *positions, fcs_float *charges
    )
{
  FCSResult result;

  P2NFFT_CHECK_RETURN_RESULT(handle, __func__);

  /* Check for periodicity */
  const fcs_int *periodicity = fcs_get_periodicity(handle);

  /* Check for correct box parameters */
  const fcs_float *a = fcs_get_box_a(handle);
  const fcs_float *b = fcs_get_box_b(handle);
  const fcs_float *c = fcs_get_box_c(handle);

  if( !nonperiodic_axes_are_principal(a, b, c, periodicity) )
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
        "Principal box vectors are required for all nonperiodic dims.");

  if( !nonperiodic_axes_are_orthogonal_to_all_other_axes(a, b, c, periodicity) )
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
        "Nonperiodic dims require box vector that are orthognonal to all box vectors of periodic dims.");

  if(periodic_dims(periodicity) == 1)
    if(!nonperiodic_box_lengths_are_equal(fcs_norm(a), fcs_norm(b), fcs_norm(c), periodicity))
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
          "The p2nfft method currently depends on equal nonperiodic box lengths with 1d-periodic boundary conditions.");

  /* Call the p2nfft solver's tuning routine */
  result = ifcs_p2nfft_tune(handle->method_context, periodicity,
      local_particles, positions, charges, a, b, c, fcs_get_offset(handle), 
      fcs_get_near_field_flag(handle));

  return result;
}

/* internal p2nfft-specific run function */
FCSResult fcs_p2nfft_run(
    FCS handle, fcs_int local_particles,
    fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials
    )
{
  FCSResult result;

  P2NFFT_CHECK_RETURN_RESULT(handle, __func__);

  result = fcs_p2nfft_tune(handle, local_particles, positions, charges);
  CHECK_RESULT_RETURN(result);
  
  fcs_int max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  result = ifcs_p2nfft_run(handle->method_context, local_particles, max_local_particles,
      positions, charges, potentials, field);

  return result;
}

/* clean-up function for p2nfft */
FCSResult fcs_p2nfft_destroy(
    FCS handle
    )
{
  ifcs_p2nfft_destroy(handle->method_context);
  fcs_set_method_context(handle, NULL);
  return FCS_RESULT_SUCCESS;
}


/*********************************************************************
 * method to check if p2nfft parameters are entered into checked FCS 
 *********************************************************************/
FCSResult fcs_p2nfft_check(
    FCS handle
    )
{
  P2NFFT_CHECK_RETURN_RESULT(handle, __func__);

  return FCS_RESULT_SUCCESS;
}


/************************************************************
 *     General P2NFFT parameter handling
 ************************************************************/
FCSResult fcs_p2nfft_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  char *param = *current;
  char *cur = *next;

  *matched = 0;

  /* P2NFFT specific parameters */
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_r_cut",            p2nfft_set_r_cut,                     FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_epsI",             p2nfft_set_epsI,                      FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC3_GOTO_NEXT("p2nfft_grid",             p2nfft_set_grid,                      FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC3_GOTO_NEXT("p2nfft_oversampled_grid", p2nfft_set_oversampled_grid,          FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_alpha",            p2nfft_set_alpha,                     FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_cao",              p2nfft_set_cao,                       FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_near",         p2nfft_set_reg_near,                  FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_near_name",    p2nfft_set_reg_near_by_name,          FCS_PARSE_VAL(fcs_p_char_t));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_far",          p2nfft_set_reg_far,                   FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_far_name",     p2nfft_set_reg_far_by_name,           FCS_PARSE_VAL(fcs_p_char_t));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_kernel",       p2nfft_set_reg_kernel,                FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_kernel_name",  p2nfft_set_reg_kernel_by_name,        FCS_PARSE_VAL(fcs_p_char_t));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_epsB",             p2nfft_set_epsB,                      FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_p",                p2nfft_set_p,                         FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_c",                p2nfft_set_c,                         FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_require_virial",   p2nfft_require_virial,                FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_intpol_order",     p2nfft_set_interpolation_order,       FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_k_cut",            p2nfft_set_k_cut,                     FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_ignore_tolerance", p2nfft_set_ignore_tolerance,          FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_ignore_field",     p2nfft_set_ignore_field,              FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_verbose_tuning",   p2nfft_set_verbose_tuning,            FCS_PARSE_VAL(fcs_int));

  /* PNFFT specific parameters */
  FCS_PARSE_IF_PARAM_THEN_FUNC3_GOTO_NEXT("pnfft_N",                 p2nfft_set_pnfft_N,                   FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC3_GOTO_NEXT("pnfft_n",                 p2nfft_set_pnfft_n,                   FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int), FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_window",            p2nfft_set_pnfft_window,              FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_window_name",       p2nfft_set_pnfft_window_by_name,      FCS_PARSE_VAL(fcs_p_char_t));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_m",                 p2nfft_set_pnfft_m,                   FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC3_GOTO_NEXT("pnfft_b",                 p2nfft_set_pnfft_b,                   FCS_PARSE_VAL(fcs_float), FCS_PARSE_VAL(fcs_float), FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_intpol_order",      p2nfft_set_pnfft_interpolation_order, FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_direct",            p2nfft_set_pnfft_direct,              FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_phi_hat",       p2nfft_set_pnfft_pre_phi_hat,         FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_fg_psi",            p2nfft_set_pnfft_fg_psi,              FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_fft_in_place",      p2nfft_set_pnfft_fft_in_place,        FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_sort_nodes",        p2nfft_set_pnfft_sort_nodes,          FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_interlaced",        p2nfft_set_pnfft_interlaced,          FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_grad_ik",           p2nfft_set_pnfft_grad_ik,             FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_psi",           p2nfft_set_pnfft_pre_psi,             FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_fg_psi",        p2nfft_set_pnfft_pre_fg_psi,          FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_full_psi",      p2nfft_set_pnfft_pre_full_psi,        FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_full_fg_psi",   p2nfft_set_pnfft_pre_full_fg_psi,     FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_real_f",            p2nfft_set_pnfft_real_f,              FCS_PARSE_VAL(fcs_int));

  /* PFFT specific parameters */
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_patience",           p2nfft_set_pfft_patience,             FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_patience_name",      p2nfft_set_pfft_patience_by_name,     FCS_PARSE_VAL(fcs_p_char_t));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_tune",               p2nfft_set_pfft_tune,                 FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_preserve_input",     p2nfft_set_pfft_preserve_input,       FCS_PARSE_VAL(fcs_int));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_print_parameters(FCS handle)
{
  fcs_int tolerance_type;
  fcs_float tolerance;
  fcs_get_tolerance(handle, &tolerance_type, &tolerance);
  if(tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
    printf("p2nfft: tolerance_type = FCS_TOLERANCE_TYPE_FIELD, tolerance = %" FCS_LMOD_FLOAT "e\n", tolerance);
  else if(tolerance_type == FCS_TOLERANCE_TYPE_POTENTIAL)
    printf("p2nfft: tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL, tolerance = %" FCS_LMOD_FLOAT "e\n", tolerance);

  return FCS_RESULT_SUCCESS;
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

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_set_resort(FCS handle, fcs_int resort)
{
  ifcs_p2nfft_set_resort(handle->method_context, resort);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_get_resort(FCS handle, fcs_int *resort)
{
  ifcs_p2nfft_get_resort(handle->method_context, resort);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_get_resort_availability(FCS handle, fcs_int *availability)
{
  ifcs_p2nfft_get_resort_availability(handle->method_context, availability);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  ifcs_p2nfft_get_resort_particles(handle->method_context, resort_particles);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_resort_ints(handle->method_context, src, dst, n, comm);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_resort_floats(handle->method_context, src, dst, n, comm);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_p2nfft_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_p2nfft_resort_bytes(handle->method_context, src, dst, n, comm);

  return FCS_RESULT_SUCCESS;
}
