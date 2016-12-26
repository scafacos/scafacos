/*
  Copyright (C) 2011, 2012, 2013, 2016 Michael Hofmann

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "wolf/wolf.h"
#include "fcs_wolf.h"


#define WOLF_CHECK_RETURN_RESULT(_h_, _f_)  do { \
  CHECK_HANDLE_RETURN_RESULT(_h_, _f_); \
  CHECK_METHOD_RETURN_RESULT(_h_, _f_, FCS_METHOD_WOLF, "wolf"); \
  } while (0)

#define WOLF_CHECK_RETURN_VAL(_h_, _f_, _v_)  do { \
  CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_); \
  CHECK_METHOD_RETURN_VAL(_h_, _f_, FCS_METHOD_WOLF, "wolf", _v_); \
  } while (0)


FCSResult fcs_wolf_init(FCS handle)
{
  fcs_wolf_context_t *ctx;
  const fcs_float default_cutoff = 0.0;
  const fcs_float default_alpha = 0.0;

  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  handle->wolf_param = malloc(sizeof(*handle->wolf_param));

  ctx = malloc(sizeof(fcs_wolf_context_t));

  ctx->comm = fcs_get_communicator(handle);

  MPI_Comm_size(ctx->comm, &ctx->comm_size);
  MPI_Comm_rank(ctx->comm, &ctx->comm_rank);

  fcs_set_method_context(handle, ctx);

  ifcs_wolf_create(&handle->wolf_param->wolf);

  fcs_wolf_set_cutoff(handle, default_cutoff);
  fcs_wolf_set_alpha(handle, default_alpha);

/*  handle->wolf_param->metallic_boundary_conditions = 1;*/

  handle->shift_positions = 0;

  handle->destroy = fcs_wolf_destroy;
  handle->set_parameter = fcs_wolf_set_parameter;
  handle->print_parameters = fcs_wolf_print_parameters;
  handle->tune = fcs_wolf_tune;
  handle->run = fcs_wolf_run;
/*  handle->set_compute_virial = fcs_wolf_require_virial;
  handle->get_virial = fcs_wolf_get_virial;*/

  handle->set_max_particle_move = fcs_wolf_set_max_particle_move;
  handle->set_resort = fcs_wolf_set_resort;
  handle->get_resort = fcs_wolf_get_resort;
  handle->get_resort_availability = fcs_wolf_get_resort_availability;
  handle->get_resort_particles = fcs_wolf_get_resort_particles;
  handle->resort_ints = fcs_wolf_resort_ints;
  handle->resort_floats = fcs_wolf_resort_floats;
  handle->resort_bytes = fcs_wolf_resort_bytes;

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_destroy(FCS handle)
{
  fcs_wolf_context_t *ctx;

  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_destroy(&handle->wolf_param->wolf);

  ctx = fcs_get_method_context(handle);

  free(ctx);

  fcs_set_method_context(handle, NULL);

  free(handle->wolf_param);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_check(FCS handle)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

	return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_tune(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_run(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  fcs_wolf_context_t *ctx;

/*  fcs_int i;
  fcs_float field_correction[3];*/

  const fcs_float *box_base, *box_a, *box_b, *box_c;
  const fcs_int *periodicity;
  fcs_int max_local_particles;

  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ctx = fcs_get_method_context(handle);

  box_base = fcs_get_box_origin(handle);
  box_a = fcs_get_box_a(handle);
  box_b = fcs_get_box_b(handle);
  box_c = fcs_get_box_c(handle);
  periodicity = fcs_get_periodicity(handle);

  max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  ifcs_wolf_set_system(&handle->wolf_param->wolf, box_base, box_a, box_b, box_c, periodicity);

  ifcs_wolf_set_particles(&handle->wolf_param->wolf, local_particles, max_local_particles, positions, charges, field, potentials);

  ifcs_wolf_run(&handle->wolf_param->wolf, ctx->comm);

/*  if (handle->wolf_param->metallic_boundary_conditions && (periodicity[0] || periodicity[1] || periodicity[2]))
  {
    fcs_compute_dipole_correction(handle, local_particles, positions, charges, 0.0, field_correction, NULL);

    for (i = 0; i < local_particles; ++i)
    {
      field[i * 3 + 0] -= field_correction[0];
      field[i * 3 + 1] -= field_correction[1];
      field[i * 3 + 2] -= field_correction[2];
    }
  }*/

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_require_virial(FCS handle, fcs_int compute_virial)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_virial(FCS handle, fcs_float *virial)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  fcs_int i;

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  for (i = 0; i < 9; ++i) virial[i] = handle->wolf_param->wolf.virial[i];

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("wolf_cutoff", wolf_set_cutoff, FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("wolf_alpha", wolf_set_alpha, FCS_PARSE_VAL(fcs_float));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_print_parameters(FCS handle)
{
  fcs_float cutoff, alpha;

  FCS_DEBUG_FUNC_INTRO(__func__);

  fcs_wolf_get_cutoff(handle, &cutoff);
  fcs_wolf_get_alpha(handle, &alpha);

  printf("wolf cutoff: %" FCS_LMOD_FLOAT "f\n", cutoff);
  printf("wolf alpha: %" FCS_LMOD_FLOAT "f\n", alpha);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_cutoff(FCS handle, fcs_float cutoff)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_set_cutoff(&handle->wolf_param->wolf, cutoff);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_cutoff(FCS handle, fcs_float *cutoff)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_get_cutoff(&handle->wolf_param->wolf, cutoff);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_alpha(FCS handle, fcs_float alpha)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_set_alpha(&handle->wolf_param->wolf, alpha);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_alpha(FCS handle, fcs_float *alpha)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_get_alpha(&handle->wolf_param->wolf, alpha);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_set_max_particle_move(&handle->wolf_param->wolf, max_particle_move);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_resort(FCS handle, fcs_int resort)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_set_resort(&handle->wolf_param->wolf, resort);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_resort(FCS handle, fcs_int *resort)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_get_resort(&handle->wolf_param->wolf, resort);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_resort_availability(FCS handle, fcs_int *availability)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_get_resort_availability(&handle->wolf_param->wolf, availability);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_get_resort_particles(&handle->wolf_param->wolf, resort_particles);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_resort_ints(&handle->wolf_param->wolf, src, dst, n, comm);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_resort_floats(&handle->wolf_param->wolf, src, dst, n, comm);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  ifcs_wolf_resort_bytes(&handle->wolf_param->wolf, src, dst, n, comm);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_setup(FCS handle, fcs_float cutoff, fcs_float alpha)
{
  FCSResult result;

  FCS_DEBUG_FUNC_INTRO(__func__);

  WOLF_CHECK_RETURN_RESULT(handle, __func__);

  result = fcs_wolf_set_cutoff(handle, cutoff);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_wolf_set_alpha(handle, alpha);
  if (result != FCS_RESULT_SUCCESS) return result;

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


void fcs_wolf_setup_f(void *handle, fcs_float cutoff, fcs_float alpha, fcs_int *return_value)
{
  FCSResult result;

  FCS_DEBUG_FUNC_INTRO(__func__);

  result = fcs_wolf_setup((FCS) handle, cutoff, alpha);

  if (result == FCS_RESULT_SUCCESS)
    *return_value = 0;
  else
    *return_value = fcs_result_get_return_code(result);

  FCS_DEBUG_FUNC_OUTRO(__func__, result);
}
