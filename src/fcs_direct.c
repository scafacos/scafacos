/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann, Rene Halver
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "direct/directc.h"
#include "fcs_direct.h"


#define DIRECT_CHECK_RETURN_RESULT(_h_, _f_)  do { \
  CHECK_HANDLE_RETURN_RESULT(_h_, _f_); \
  CHECK_METHOD_RETURN_RESULT(_h_, _f_, FCS_METHOD_DIRECT, "direct"); \
  } while (0)

#define DIRECT_CHECK_RETURN_VAL(_h_, _f_, _v_)  do { \
  CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_); \
  CHECK_METHOD_RETURN_VAL(_h_, _f_, FCS_METHOD_DIRECT, "direct", _v_); \
  } while (0)


FCSResult fcs_direct_init(FCS handle)
{
  fcs_direct_context_t *ctx;
  fcs_float default_cutoff = 0.0;
  fcs_int default_periodic_images[3] = { 1, 1, 1 };

  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  handle->direct_param = malloc(sizeof(*handle->direct_param));

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

  handle->shift_positions = 0;

  handle->destroy = fcs_direct_destroy;
  handle->set_parameter = fcs_direct_set_parameter;
  handle->print_parameters = fcs_direct_print_parameters;
  handle->tune = fcs_direct_tune;
  handle->run = fcs_direct_run;
  handle->set_compute_virial = fcs_direct_require_virial;
  handle->get_virial = fcs_direct_get_virial;

  handle->set_max_particle_move = fcs_direct_set_max_particle_move;
  handle->set_resort = fcs_direct_set_resort;
  handle->get_resort = fcs_direct_get_resort;
  handle->get_resort_availability = fcs_direct_get_resort_availability;
  handle->get_resort_particles = fcs_direct_get_resort_particles;
  handle->resort_ints = fcs_direct_resort_ints;
  handle->resort_floats = fcs_direct_resort_floats;
  handle->resort_bytes = fcs_direct_resort_bytes;

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_destroy(FCS handle)
{
  fcs_direct_context_t *ctx;

  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_destroy(&handle->direct_param->directc);

  ctx = fcs_get_method_context(handle);

  free(ctx);

  fcs_set_method_context(handle, NULL);

  free(handle->direct_param);
  handle->direct_param = NULL;

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_check(FCS handle)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

	return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_tune(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_run(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  fcs_direct_context_t *ctx;
  fcs_int i;
  fcs_float field_correction[3];
  const fcs_float *box_base, *box_a, *box_b, *box_c;
  const fcs_int *periodicity;
  fcs_int max_local_particles;
  fcs_float cutoff = 0;

  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  ctx = fcs_get_method_context(handle);

  box_base = fcs_get_box_origin(handle);
  box_a = fcs_get_box_a(handle);
  box_b = fcs_get_box_b(handle);
  box_c = fcs_get_box_c(handle);
  periodicity = fcs_get_periodicity(handle);

  max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  fcs_directc_set_system(&handle->direct_param->directc, box_base, box_a, box_b, box_c, periodicity);

  fcs_directc_set_particles(&handle->direct_param->directc, local_particles, max_local_particles, positions, charges, field, potentials);

  fcs_directc_run(&handle->direct_param->directc, ctx->comm);

  fcs_directc_get_cutoff(&handle->direct_param->directc, &cutoff);

  if (handle->direct_param->metallic_boundary_conditions && cutoff == 0 && (periodicity[0] || periodicity[1] || periodicity[2]))
  {
    fcs_compute_dipole_correction(handle, local_particles, positions, charges, 0.0, field_correction, NULL);

    for (i = 0; i < local_particles; ++i)
    {
      field[i * 3 + 0] -= field_correction[0];
      field[i * 3 + 1] -= field_correction[1];
      field[i * 3 + 2] -= field_correction[2];
    }
  }

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_periodic_images",  direct_set_periodic_images,  FCS_PARSE_SEQ(fcs_int, 3));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_cutoff",           direct_set_cutoff,           FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_cutoff_with_near", direct_set_cutoff_with_near, FCS_PARSE_VAL(fcs_int));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_print_parameters(FCS handle)
{
  fcs_float cutoff;
  fcs_int images[3];
  FCSResult result;

  FCS_DEBUG_FUNC_INTRO(__func__);

  result = fcs_direct_get_cutoff(handle, &cutoff);
  if (result != FCS_RESULT_SUCCESS)
  {
    printf("direct cutoff: FAILED!");
    fcs_result_print_result(result);
    fcs_result_destroy(result);
  } else printf("direct cutoff: %" FCS_LMOD_FLOAT "e\n", cutoff);

  result = fcs_direct_get_periodic_images(handle, images);
  if (result != FCS_RESULT_SUCCESS)
  {
    printf("direct cutoff: FAILED!");
    fcs_result_print_result(result);
    fcs_result_destroy(result);
  } else printf("direct periodic images: %5" FCS_LMOD_INT "d %5" FCS_LMOD_INT "d %5" FCS_LMOD_INT "d\n", images[0], images[1], images[2]);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_require_virial(FCS handle, fcs_int compute_virial)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_virial(FCS handle, fcs_float *virial)
{
  fcs_int i;

  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  for (i = 0; i < 9; ++i) virial[i] = handle->direct_param->directc.virial[i];

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_cutoff(FCS handle, fcs_float cutoff)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_set_cutoff(&handle->direct_param->directc, cutoff);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_cutoff(FCS handle, fcs_float *cutoff)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_get_cutoff(&handle->direct_param->directc, cutoff);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_cutoff_with_near(FCS handle, fcs_bool cutoff_with_near)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_set_cutoff_with_near(&handle->direct_param->directc, FCS_IS_TRUE(cutoff_with_near));

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_cutoff_with_near(FCS handle, fcs_bool *cutoff_with_near)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_int i;
  fcs_directc_get_cutoff_with_near(&handle->direct_param->directc, &i);

  *cutoff_with_near = (i)?FCS_TRUE:FCS_FALSE;

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_periodic_images(FCS handle, fcs_int *periodic_images)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_set_periodic_images(&handle->direct_param->directc, periodic_images);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_periodic_images(FCS handle, fcs_int *periodic_images)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_get_periodic_images(&handle->direct_param->directc, periodic_images);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_in_particles(FCS handle, fcs_int nin_particles, fcs_float *in_positions, fcs_float *in_charges)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_set_in_particles(&handle->direct_param->directc, nin_particles, in_positions, in_charges);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_set_max_particle_move(&handle->direct_param->directc, max_particle_move);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_set_resort(FCS handle, fcs_int resort)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_set_resort(&handle->direct_param->directc, resort);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_resort(FCS handle, fcs_int *resort)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_get_resort(&handle->direct_param->directc, resort);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_resort_availability(FCS handle, fcs_int *availability)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_get_resort_availability(&handle->direct_param->directc, availability);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_get_resort_particles(&handle->direct_param->directc, resort_particles);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_resort_ints(&handle->direct_param->directc, src, dst, n, comm);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_resort_floats(&handle->direct_param->directc, src, dst, n, comm);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_resort_bytes(FCS handle, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  fcs_directc_resort_bytes(&handle->direct_param->directc, src, dst, n, comm);

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);
  
  return FCS_RESULT_SUCCESS;
}


FCSResult fcs_direct_setup(FCS handle, fcs_float cutoff)
{
  FCSResult result;

  FCS_DEBUG_FUNC_INTRO(__func__);

  DIRECT_CHECK_RETURN_RESULT(handle, __func__);

  result = fcs_direct_set_cutoff(handle, cutoff);
  if (result != FCS_RESULT_SUCCESS)
  {
    FCS_DEBUG_FUNC_OUTRO(__func__, result);
    return result;
  }

  FCS_DEBUG_FUNC_OUTRO(__func__, FCS_RESULT_SUCCESS);

  return FCS_RESULT_SUCCESS;
}


void fcs_direct_setup_f(void *handle, fcs_float cutoff, fcs_int *return_value)
{
  FCSResult result;

  FCS_DEBUG_FUNC_INTRO(__func__);

  result = fcs_direct_setup((FCS) handle, cutoff);

  if (result == FCS_RESULT_SUCCESS)
    *return_value = 0;
  else
    *return_value = fcs_result_get_return_code(result);

  FCS_DEBUG_FUNC_OUTRO(__func__, result);
}
