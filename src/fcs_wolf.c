/*
  Copyright (C) 2011,2012,2013 Michael Hofmann

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


#ifdef FCS_ENABLE_DEBUG
# define DEBUG_MOP(_mop_)  do { _mop_; } while (0)
#else
# define DEBUG_MOP(_mop_)  do {  } while (0)
#endif


#define WOLF_HANDLE_CHECK(_h_, _f_) do { \
  if ((_h_) == FCS_NULL) \
    return fcsResult_create(FCS_NULL_ARGUMENT, (_f_), "null pointer supplied as handle"); \
  if ((_h_)->method != FCS_WOLF) \
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, (_f_), "handle does not represent method \"wolf\""); \
} while (0)


FCSResult fcs_wolf_init(FCS handle)
{
  static const char func_name[] = "fcs_wolf_init";

  fcs_wolf_context_t *ctx;
  const fcs_float default_cutoff = 0.0;
  const fcs_float default_alpha = 0.0;


  DEBUG_MOP(printf("fcs_wolf_init\n"));

  WOLF_HANDLE_CHECK(handle, func_name);

  ctx = malloc(sizeof(fcs_wolf_context_t));

  ctx->comm = fcs_get_communicator(handle);

  MPI_Comm_size(ctx->comm, &ctx->comm_size);
  MPI_Comm_rank(ctx->comm, &ctx->comm_rank);

  fcs_set_method_context(handle, ctx);

  ifcs_wolf_create(&handle->wolf_param->wolf);

  fcs_wolf_set_cutoff(handle, default_cutoff);
  fcs_wolf_set_alpha(handle, default_alpha);

/*  handle->wolf_param->metallic_boundary_conditions = 1;*/

  handle->set_max_particle_move = fcs_wolf_set_max_particle_move;
  handle->set_resort = fcs_wolf_set_resort;
  handle->get_resort = fcs_wolf_get_resort;
  handle->get_resort_availability = fcs_wolf_get_resort_availability;
  handle->get_resort_particles = fcs_wolf_get_resort_particles;
  handle->resort_ints = fcs_wolf_resort_ints;
  handle->resort_floats = fcs_wolf_resort_floats;
  handle->resort_bytes = fcs_wolf_resort_bytes;

  DEBUG_MOP(printf("fcs_wolf_init: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_destroy(FCS handle)
{
  static const char func_name[] = "fcs_wolf_destroy";

  fcs_wolf_context_t *ctx;


  DEBUG_MOP(printf("fcs_wolf_destroy\n"));

  WOLF_HANDLE_CHECK(handle, func_name);

  ifcs_wolf_destroy(&handle->wolf_param->wolf);

  ctx = fcs_get_method_context(handle);

  free(ctx);

  fcs_set_method_context(handle, NULL);

  DEBUG_MOP(printf("fcs_wolf_destroy: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_check(FCS handle)
{
  static const char func_name[] = "fcs_wolf_check";

  DEBUG_MOP(printf("fcs_wolf_check\n"));

  WOLF_HANDLE_CHECK(handle, func_name);

  DEBUG_MOP(printf("fcs_wolf_check: done\n"));

	return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges)
{
  static const char func_name[] = "fcs_wolf_tune";

  DEBUG_MOP(printf("fcs_wolf_tune\n"));

  WOLF_HANDLE_CHECK(handle, func_name);

  DEBUG_MOP(printf("fcs_wolf_tune: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  static const char func_name[] = "fcs_wolf_run";

  fcs_wolf_context_t *ctx;

/*  fcs_int i;
  fcs_float field_correction[3];*/

  fcs_float *box_base, *box_a, *box_b, *box_c;
  fcs_int *periodicity;


  DEBUG_MOP(printf("fcs_wolf_run\n"));

  WOLF_HANDLE_CHECK(handle, func_name);

  ctx = fcs_get_method_context(handle);

  box_base = fcs_get_offset(handle);
  box_a = fcs_get_box_a(handle);
  box_b = fcs_get_box_b(handle);
  box_c = fcs_get_box_c(handle);
  periodicity = fcs_get_periodicity(handle);

  ifcs_wolf_set_system(&handle->wolf_param->wolf, box_base, box_a, box_b, box_c, periodicity);

  ifcs_wolf_set_particles(&handle->wolf_param->wolf, local_particles, local_max_particles, positions, charges, field, potentials);

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

  DEBUG_MOP(printf("fcs_wolf_run: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_require_virial(FCS handle, fcs_int compute_virial)
{
  static const char func_name[] = "fcs_wolf_require_virial";

  WOLF_HANDLE_CHECK(handle, func_name);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_virial(FCS handle, fcs_float *virial)
{
/*  static const char func_name[] = "fcs_wolf_get_virial";

  fcs_int i;

  WOLF_HANDLE_CHECK(handle, func_name);

  for (i = 0; i < 9; ++i) virial[i] = handle->wolf_param->wolf.virial[i];*/

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_cutoff(FCS handle, fcs_float cutoff)
{
  static const char func_name[] = "fcs_wolf_set_cutoff";

  WOLF_HANDLE_CHECK(handle, func_name);

  ifcs_wolf_set_cutoff(&handle->wolf_param->wolf, cutoff);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_cutoff(FCS handle, fcs_float *cutoff)
{
  static const char func_name[] = "fcs_wolf_get_cutoff";

  WOLF_HANDLE_CHECK(handle, func_name);

  ifcs_wolf_get_cutoff(&handle->wolf_param->wolf, cutoff);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_alpha(FCS handle, fcs_float alpha)
{
  static const char func_name[] = "fcs_wolf_set_alpha";

  WOLF_HANDLE_CHECK(handle, func_name);

  ifcs_wolf_set_alpha(&handle->wolf_param->wolf, alpha);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_alpha(FCS handle, fcs_float *alpha)
{
  static const char func_name[] = "fcs_wolf_get_alpha";

  WOLF_HANDLE_CHECK(handle, func_name);

  ifcs_wolf_get_alpha(&handle->wolf_param->wolf, alpha);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  ifcs_wolf_set_max_particle_move(&handle->wolf_param->wolf, max_particle_move);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_set_resort(FCS handle, fcs_int resort)
{
  ifcs_wolf_set_resort(&handle->wolf_param->wolf, resort);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_resort(FCS handle, fcs_int *resort)
{
  ifcs_wolf_get_resort(&handle->wolf_param->wolf, resort);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_resort_availability(FCS handle, fcs_int *availability)
{
  ifcs_wolf_get_resort_availability(&handle->wolf_param->wolf, availability);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  ifcs_wolf_get_resort_particles(&handle->wolf_param->wolf, resort_particles);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_wolf_resort_ints(&handle->wolf_param->wolf, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_wolf_resort_floats(&handle->wolf_param->wolf, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  ifcs_wolf_resort_bytes(&handle->wolf_param->wolf, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_wolf_setup(FCS handle, fcs_float cutoff, fcs_float alpha)
{
  static const char func_name[] = "fcs_wolf_setup";

  FCSResult result;

  WOLF_HANDLE_CHECK(handle, func_name);

  result = fcs_wolf_set_cutoff(handle, cutoff);
  if (result != FCS_RESULT_SUCCESS) return result;

  result = fcs_wolf_set_alpha(handle, alpha);
  if (result != FCS_RESULT_SUCCESS) return result;

  return FCS_RESULT_SUCCESS;
}


void fcs_wolf_setup_f(void *handle, fcs_float cutoff, fcs_float alpha, fcs_int *return_value)
{
  FCSResult result;

  result = fcs_wolf_setup((FCS) handle, cutoff, alpha);

  if (result == FCS_RESULT_SUCCESS)
    *return_value = 0;
  else
    *return_value = fcsResult_getReturnCode(result);
}
