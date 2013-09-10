/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann, Rene Halver

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

#include "direct/directc.h"
#include "fcs_direct.h"


#ifdef FCS_ENABLE_DEBUG
# define DEBUG_MOP(_mop_)  do { _mop_; } while (0)
#else
# define DEBUG_MOP(_mop_)  do {  } while (0)
#endif


#define DIRECT_HANDLE_CHECK(_h_, _f_) do { \
  if ((_h_) == FCS_NULL) \
    return fcsResult_create(FCS_NULL_ARGUMENT, (_f_), "null pointer supplied as handle"); \
  if ((_h_)->method != FCS_DIRECT) \
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, (_f_), "handle does not represent method \"direct\""); \
} while (0)


FCSResult fcs_direct_init(FCS handle)
{
  static const char func_name[] = "fcs_direct_init";

  fcs_direct_context_t *ctx;
  fcs_float default_cutoff = 0.0;
  fcs_int default_periodic_images[3] = { 1, 1, 1 };


  DEBUG_MOP(printf("fcs_direct_init\n"));

  DIRECT_HANDLE_CHECK(handle, func_name);

  ctx = malloc(sizeof(fcs_direct_context_t));

  ctx->comm = fcs_get_communicator(handle);

  MPI_Comm_size(ctx->comm, &ctx->comm_size);
  MPI_Comm_rank(ctx->comm, &ctx->comm_rank);

  fcs_set_method_context(handle, ctx);

  fcs_directc_create(&handle->direct_param->directc);

  fcs_direct_set_cutoff(handle, default_cutoff);

  fcs_direct_set_periodic_images(handle, default_periodic_images);
  
  fcs_direct_set_cutoff_with_near(handle, FCS_FALSE);

  handle->direct_param->metallic_boundary_conditions = 1;

  handle->set_max_particle_move = fcs_direct_set_max_particle_move;
  handle->set_resort = fcs_direct_set_resort;
  handle->get_resort = fcs_direct_get_resort;
  handle->get_resort_availability = fcs_direct_get_resort_availability;
  handle->get_resort_particles = fcs_direct_get_resort_particles;
  handle->resort_ints = fcs_direct_resort_ints;
  handle->resort_floats = fcs_direct_resort_floats;
  handle->resort_bytes = fcs_direct_resort_bytes;

  DEBUG_MOP(printf("fcs_direct_init: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_destroy(FCS handle)
{
  static const char func_name[] = "fcs_direct_destroy";

  fcs_direct_context_t *ctx;


  DEBUG_MOP(printf("fcs_direct_destroy\n"));

  DIRECT_HANDLE_CHECK(handle, func_name);

  fcs_directc_destroy(&handle->direct_param->directc);

  ctx = fcs_get_method_context(handle);

  free(ctx);

  fcs_set_method_context(handle, NULL);

  DEBUG_MOP(printf("fcs_direct_destroy: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_check(FCS handle)
{
  static const char func_name[] = "fcs_direct_check";

  DEBUG_MOP(printf("fcs_direct_check\n"));

  DIRECT_HANDLE_CHECK(handle, func_name);

  DEBUG_MOP(printf("fcs_direct_check: done\n"));

	return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges)
{
  static const char func_name[] = "fcs_direct_tune";

  DEBUG_MOP(printf("fcs_direct_tune\n"));

  DIRECT_HANDLE_CHECK(handle, func_name);

  DEBUG_MOP(printf("fcs_direct_tune: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  static const char func_name[] = "fcs_direct_run";

  fcs_direct_context_t *ctx;

  fcs_int i;

  fcs_float field_correction[3];

  fcs_float *box_base, *box_a, *box_b, *box_c;
  fcs_int *periodicity;


  DEBUG_MOP(printf("fcs_direct_run\n"));

  DIRECT_HANDLE_CHECK(handle, func_name);

  ctx = fcs_get_method_context(handle);

  box_base = fcs_get_offset(handle);
  box_a = fcs_get_box_a(handle);
  box_b = fcs_get_box_b(handle);
  box_c = fcs_get_box_c(handle);
  periodicity = fcs_get_periodicity(handle);

  fcs_directc_set_system(&handle->direct_param->directc, box_base, box_a, box_b, box_c, periodicity);

  fcs_directc_set_particles(&handle->direct_param->directc, local_particles, local_max_particles, positions, charges, field, potentials);

  fcs_directc_run(&handle->direct_param->directc, ctx->comm);

  if (handle->direct_param->metallic_boundary_conditions && (periodicity[0] || periodicity[1] || periodicity[2]))
  {
    fcs_compute_dipole_correction(handle, local_particles, positions, charges, 0.0, field_correction, NULL);

    for (i = 0; i < local_particles; ++i)
    {
      field[i * 3 + 0] -= field_correction[0];
      field[i * 3 + 1] -= field_correction[1];
      field[i * 3 + 2] -= field_correction[2];
    }
  }

  DEBUG_MOP(printf("fcs_direct_run: done\n"));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_require_virial(FCS handle, fcs_int compute_virial)
{
  static const char func_name[] = "fcs_direct_require_virial";

  DIRECT_HANDLE_CHECK(handle, func_name);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_virial(FCS handle, fcs_float *virial)
{
  static const char func_name[] = "fcs_direct_get_virial";

  fcs_int i;

  DIRECT_HANDLE_CHECK(handle, func_name);

  for (i = 0; i < 9; ++i) virial[i] = handle->direct_param->directc.virial[i];

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_cutoff(FCS handle, fcs_float cutoff)
{
  static const char func_name[] = "fcs_direct_set_cutoff";

  DIRECT_HANDLE_CHECK(handle, func_name);

  fcs_directc_set_cutoff(&handle->direct_param->directc, cutoff);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_cutoff(FCS handle, fcs_float *cutoff)
{
  static const char func_name[] = "fcs_direct_get_cutoff";

  DIRECT_HANDLE_CHECK(handle, func_name);

  fcs_directc_get_cutoff(&handle->direct_param->directc, cutoff);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_cutoff_with_near(FCS handle, fcs_bool cutoff_with_near)
{
  static const char func_name[] = "fcs_direct_set_cutoff_with_near";

  DIRECT_HANDLE_CHECK(handle, func_name);

  fcs_directc_set_cutoff_with_near(&handle->direct_param->directc, FCS_IS_TRUE(cutoff_with_near));

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_cutoff_with_near(FCS handle, fcs_bool *cutoff_with_near)
{
  static const char func_name[] = "fcs_direct_get_cutoff_with_near";

  DIRECT_HANDLE_CHECK(handle, func_name);

  *cutoff_with_near = handle->direct_param->directc.cutoff_with_near;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_periodic_images(FCS handle, fcs_int *periodic_images)
{
  static const char func_name[] = "fcs_direct_set_periodic_images";

  DIRECT_HANDLE_CHECK(handle, func_name);

  fcs_directc_set_periodic_images(&handle->direct_param->directc, periodic_images);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_periodic_images(FCS handle, fcs_int *periodic_images)
{
  static const char func_name[] = "fcs_direct_get_periodic_images";

  DIRECT_HANDLE_CHECK(handle, func_name);

  periodic_images[0] = handle->direct_param->directc.periodic_images[0];
  periodic_images[1] = handle->direct_param->directc.periodic_images[1];
  periodic_images[2] = handle->direct_param->directc.periodic_images[2];

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_in_particles(FCS handle, fcs_int nin_particles, fcs_float *in_positions, fcs_float *in_charges)
{
  static const char func_name[] = "fcs_directc_set_in_particles";

  DIRECT_HANDLE_CHECK(handle, func_name);

  fcs_directc_set_in_particles(&handle->direct_param->directc, nin_particles, in_positions, in_charges);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  fcs_directc_set_max_particle_move(&handle->direct_param->directc, max_particle_move);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_resort(FCS handle, fcs_int resort)
{
  fcs_directc_set_resort(&handle->direct_param->directc, resort);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_resort(FCS handle, fcs_int *resort)
{
  fcs_directc_get_resort(&handle->direct_param->directc, resort);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_resort_availability(FCS handle, fcs_int *availability)
{
  fcs_directc_get_resort_availability(&handle->direct_param->directc, availability);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  fcs_directc_get_resort_particles(&handle->direct_param->directc, resort_particles);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  fcs_directc_resort_ints(&handle->direct_param->directc, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  fcs_directc_resort_floats(&handle->direct_param->directc, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  fcs_directc_resort_bytes(&handle->direct_param->directc, src, dst, n, comm);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_setup(FCS handle, fcs_float cutoff)
{
  static const char func_name[] = "fcs_direct_setup";

  FCSResult result;

  DIRECT_HANDLE_CHECK(handle, func_name);

  result = fcs_direct_set_cutoff(handle, cutoff);
  if (result != FCS_RESULT_SUCCESS) return result;

  return FCS_RESULT_SUCCESS;
}


void fcs_direct_setup_f(void *handle, fcs_float cutoff, fcs_int *return_value)
{
  FCSResult result;

  result = fcs_direct_setup((FCS)handle, cutoff);

  if (result == FCS_RESULT_SUCCESS)
    *return_value = 0;
  else
    *return_value = fcsResult_getReturnCode(result);
}
