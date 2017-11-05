/*
  Copyright (C) 2011-2012 Rene Halver, Lukas Arnold, Mathias Winkel
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

#include <string.h>
#include <mpi.h>

#include "fcs_pepc.h"
#include "FCSCommon.h"


#define PEPC_CHECK_RETURN_RESULT(_h_, _f_)  do { \
  CHECK_HANDLE_RETURN_RESULT(_h_, _f_); \
  CHECK_METHOD_RETURN_RESULT(_h_, _f_, FCS_METHOD_PEPC, "pepc"); \
  } while (0)

#define PEPC_CHECK_RETURN_VAL(_h_, _f_, _v_)  do { \
  CHECK_HANDLE_RETURN_VAL(_h_, _f_, _v_); \
  CHECK_METHOD_RETURN_VAL(_h_, _f_, FCS_METHOD_PEPC, "pepc", _v_); \
  } while (0)


/* combined setter function for all pepc parameters */
FCSResult fcs_pepc_setup(FCS handle, fcs_float epsilon, fcs_float theta)
{
  FCSResult result;

  result = fcs_pepc_set_epsilon(handle, epsilon);
  CHECK_RESULT_RETURN(result);

  result = fcs_pepc_set_theta(handle, theta);
  CHECK_RESULT_RETURN(result);

  return FCS_RESULT_SUCCESS;

}

/* method to check if pepc parameters are entered into checked FCS */
FCSResult fcs_pepc_check(FCS handle)
{
  const fcs_float *a,*b,*c;
  fcs_float eps, theta;
  FCSResult res;

  if ((res = fcs_pepc_get_epsilon(handle, &eps))) return res;
  if (eps == -1.0)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "pepc: epsilon not set");
  if ((res = fcs_pepc_get_theta(handle, &theta))) return res;
  if (theta == -1.0)
    return fcs_result_create(FCS_ERROR_MISSING_ELEMENT, __func__, "pepc: theta not set");

  a = fcs_get_box_a(handle);
  b = fcs_get_box_b(handle);
  c = fcs_get_box_c(handle);
  if (!fcs_uses_principal_axes(a,b,c))
    printf("%s\n", "WARNING: support of pepc for non-cubic simulation boxes currently is experimental.");

  if (!fcs_get_near_field_flag(handle))
    return fcs_result_create(FCS_ERROR_INCOMPATIBLE_METHOD, __func__, 
			    "pepc performs near field computations by itself!");

  return FCS_RESULT_SUCCESS;
}



/* initialization function for basic pepc parameters */
FCSResult fcs_pepc_init(FCS handle)
{
  handle->shift_positions = 1;

  handle->destroy = fcs_pepc_destroy;
  handle->set_parameter = fcs_pepc_set_parameter;
  handle->print_parameters = fcs_pepc_print_parameters;
  handle->tune = fcs_pepc_tune;
  handle->run = fcs_pepc_run;
  handle->set_compute_virial = fcs_pepc_require_virial;
  handle->get_virial = fcs_pepc_get_virial;

  handle->pepc_param = malloc(sizeof(*handle->pepc_param));
  handle->pepc_param->theta             = 0.6;
  handle->pepc_param->epsilon           = 0.0;
  handle->pepc_param->require_virial    = 0;
  handle->pepc_param->num_walk_threads  = 3;
  handle->pepc_param->load_balancing    = 0;
  handle->pepc_param->dipole_correction = 1;
  handle->pepc_param->npm               = -45.0;
  handle->pepc_param->debug_level       = 0;

  fcs_pepc_internal_t *pepc_internal;
  MPI_Comm comm  = fcs_get_communicator(handle);
  MPI_Fint fcomm = MPI_Comm_c2f(comm);

  pepc_scafacos_initialize(&fcomm);

  handle->method_context = malloc(sizeof(fcs_pepc_internal_t));
  pepc_internal = (fcs_pepc_internal_t*) handle->method_context;
  pepc_internal->work_length = -1;
  pepc_internal->work        = NULL;

  return FCS_RESULT_SUCCESS;
}

/* internal pepc-specific tuning function */
FCSResult fcs_pepc_tune(FCS handle, fcs_int local_particles, fcs_float *positions, fcs_float *charges)
{
  FCSResult result;

  /* currently, there is no tune functionality for pepc */
  result = fcs_pepc_check(handle);
  return result;
}

/* internal pepc-specific run function */
FCSResult fcs_pepc_run(FCS handle, fcs_int local_particles,
		       fcs_float *positions, fcs_float *charges, 
		       fcs_float *field, fcs_float *potentials)
{
  FCSResult result;
  fcs_int pcnt;
  fcs_int nparts_tot;
  fcs_pepc_internal_t *pepc_internal;

  nparts_tot = fcs_get_total_particles(handle);

  pepc_internal = (fcs_pepc_internal_t*) fcs_get_method_context(handle);
  
  if((local_particles != pepc_internal->work_length) || (handle->pepc_param->load_balancing==0)){
    
    if(pepc_internal->work != NULL){
      free(pepc_internal->work);
    }
    
    pepc_internal->work = (fcs_float*)malloc(sizeof(fcs_float)*local_particles);
    pepc_internal->work_length = local_particles;

    for(int pcnt=0; pcnt<pepc_internal->work_length; pcnt++)
      pepc_internal->work[pcnt] = 1.0;
  }
  
  result = fcs_pepc_check(handle);
  CHECK_RESULT_RETURN(result);

  fcs_int max_local_particles = fcs_get_max_local_particles(handle);
  if (local_particles > max_local_particles) max_local_particles = local_particles;

  if (handle->pepc_param->debug_level > 3) {
    printf("*** run pepc kernel\n");
    printf("** local particles:        %" FCS_LMOD_INT "d\n", local_particles);
    printf("** max local particles:    %" FCS_LMOD_INT "d\n", max_local_particles);
    printf("** total particles:        %" FCS_LMOD_INT "d\n", nparts_tot);
    printf("** epsilon:                %" FCS_LMOD_FLOAT "f\n", handle->pepc_param->epsilon);
    printf("** theta:                  %" FCS_LMOD_FLOAT "f\n", handle->pepc_param->theta);
    printf("** npm:                    %" FCS_LMOD_FLOAT "f\n", handle->pepc_param->npm);
    printf("** debug_level:            %" FCS_LMOD_INT "d\n", handle->pepc_param->debug_level);
    printf("** require virial:         %" FCS_LMOD_INT "d\n", handle->pepc_param->require_virial);
    printf("** num walk threads:       %" FCS_LMOD_INT "d\n", handle->pepc_param->num_walk_threads);
    printf("** dipole correction:      %" FCS_LMOD_INT "d\n", handle->pepc_param->dipole_correction);
    printf("** use load balancing:     %" FCS_LMOD_INT "d\n", handle->pepc_param->load_balancing);
    printf("** size int:               %d\n", (int)sizeof(fcs_int));
    printf("** size float:             %d\n", (int)sizeof(fcs_float));
    printf("** debug lattice pointers: %p\n", fcs_get_box_a(handle));
    printf("** debug lattice pointers: %p\n", fcs_get_box_b(handle));
    printf("** debug lattice pointers: %p\n", fcs_get_box_c(handle));

    printf("** debug output lattice in x %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f\n", 
	   fcs_get_box_a(handle)[0], fcs_get_box_a(handle)[1], fcs_get_box_a(handle)[2]);
    printf("** debug output lattice in y %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f\n", 
	   fcs_get_box_b(handle)[0], fcs_get_box_b(handle)[1], fcs_get_box_b(handle)[2]);
    printf("** debug output lattice in z %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f %" FCS_LMOD_FLOAT "f\n", 
	   fcs_get_box_c(handle)[0], fcs_get_box_c(handle)[1], fcs_get_box_c(handle)[2]);
  }

  pepc_scafacos_run(&local_particles, &nparts_tot, positions, charges,
		    field, potentials, pepc_internal->work,
		    ((fcs_pepc_internal_t*)(handle->method_context))->virial,
		    fcs_get_box_a(handle), fcs_get_box_b(handle), fcs_get_box_c(handle),
		    fcs_get_periodicity(handle), &handle->pepc_param->dipole_correction,
		    &handle->pepc_param->epsilon, &handle->pepc_param->theta, &handle->pepc_param->debug_level, &handle->pepc_param->num_walk_threads, &handle->pepc_param->npm);

  if (handle->pepc_param->debug_level > 3)
  {
    printf("virial(0,0) %12.4e\n", ((fcs_pepc_internal_t*)(handle->method_context))->virial[0]);

    for(pcnt=0; pcnt<local_particles; pcnt++)
      printf("** e-field dump, particle %" FCS_LMOD_INT "d, (ex,ey,ez): %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f\n", pcnt, field[pcnt+0], field[pcnt+1], field[pcnt+2]);
  }

  return FCS_RESULT_SUCCESS;
}

/* clean-up function for pepc */
FCSResult fcs_pepc_destroy(FCS handle)
{
  MPI_Comm comm  = fcs_get_communicator(handle);
  MPI_Fint fcomm = MPI_Comm_c2f(comm);
  pepc_scafacos_finalize(&fcomm);

  if(((fcs_pepc_internal_t*)fcs_get_method_context(handle))->work != NULL){
    free(((fcs_pepc_internal_t*)fcs_get_method_context(handle))->work);   
  }

  free(handle->method_context);

  free(handle->pepc_param);

  return FCS_RESULT_SUCCESS;
}

/******************************************************************************************************
 *
 *						Setter and Getter functions for pepc parameters
 *
 ******************************************************************************************************/


/* setter function for pepc parameter epsilon */
FCSResult fcs_pepc_set_epsilon(FCS handle, fcs_float epsilon)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  if (epsilon < 0 )
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
			    "epsilon > 0. has been violated");

  handle->pepc_param->epsilon = epsilon;

  return FCS_RESULT_SUCCESS;
}

/* getter function for pepc parameter epsilon */
FCSResult fcs_pepc_get_epsilon(FCS handle, fcs_float* epsilon)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  *epsilon = handle->pepc_param->epsilon;

  return FCS_RESULT_SUCCESS;
}


/* setter function for pepc parameter theta */
FCSResult fcs_pepc_set_theta(FCS handle, fcs_float theta)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  if (theta < 0)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
			    "0 <= theta has been violated");

  handle->pepc_param->theta = theta;

  return FCS_RESULT_SUCCESS;
}

/* getter function for pepc parameter theta */
FCSResult fcs_pepc_get_theta(FCS handle, fcs_float* theta)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  *theta = handle->pepc_param->theta;

  return FCS_RESULT_SUCCESS;
}

/* getter function for pepc parameter num_walk_threads */
FCSResult fcs_pepc_get_num_walk_threads(FCS handle, fcs_int *num_walk_threads)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  *num_walk_threads = handle->pepc_param->num_walk_threads;

  return FCS_RESULT_SUCCESS;
}

/* setter function for pepc parameter num_walk_threads */
FCSResult fcs_pepc_set_num_walk_threads(FCS handle, fcs_int num_walk_threads)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  if (num_walk_threads < 1 )
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, __func__,
			    "0 < num_walk_threads has been violated");

  handle->pepc_param->num_walk_threads = num_walk_threads;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_pepc_set_load_balancing(FCS handle, fcs_int load_balancing)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  handle->pepc_param->load_balancing = load_balancing;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_pepc_get_load_balancing(FCS handle, fcs_int* load_balancing)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  *load_balancing = handle->pepc_param->load_balancing;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_pepc_set_dipole_correction(FCS handle, fcs_int dipole_correction)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  handle->pepc_param->dipole_correction = dipole_correction;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_pepc_get_dipole_correction(FCS handle, fcs_int* dipole_correction)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  *dipole_correction = handle->pepc_param->dipole_correction;

  return FCS_RESULT_SUCCESS;
}

/* setter function for pepc parameter npm */
FCSResult fcs_pepc_set_npm(FCS handle, fcs_float npm)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  handle->pepc_param->npm = npm;

  return FCS_RESULT_SUCCESS;
}

/* getter function for pepc parameter npm */
FCSResult fcs_pepc_get_npm(FCS handle, fcs_float* npm)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  *npm = handle->pepc_param->npm;

  return FCS_RESULT_SUCCESS;
}

/* setter function for pepc parameter debug_level */
FCSResult fcs_pepc_set_debug_level(FCS handle, fcs_int level)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  handle->pepc_param->debug_level = level;

  return FCS_RESULT_SUCCESS;
}

/* getter function for pepc parameter debug_level */
FCSResult fcs_pepc_get_debug_level(FCS handle, fcs_int* level)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  *level = handle->pepc_param->debug_level;

  return FCS_RESULT_SUCCESS;
}

/******************************************************************************************************
 *
 *						additional pepc functions
 *
 ******************************************************************************************************/

/* setter function to (de)activate virial computation in pepc */
FCSResult fcs_pepc_require_virial(FCS handle, fcs_int choice)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  if ((choice < 0) || (choice > 1))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"require_virial must be set to 0 or 1");

  handle->pepc_param->require_virial = choice;

  return FCS_RESULT_SUCCESS;
}

/* function getting the virial result from pepc */
FCSResult fcs_pepc_get_virial(FCS handle, fcs_float virial[9])
{
  int i;

  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  if (handle->pepc_param->require_virial != 1)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT,__func__,"calculation of virial was not activated. Please call require_virial(fcs, 1) before calling run()");

  for (i=0;i<9;i++)
	  virial[i] = ((fcs_pepc_internal_t*)(handle->method_context))->virial[i];

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_pepc_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_epsilon",           pepc_set_epsilon,           FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_theta",             pepc_set_theta,             FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_num_walk_threads",  pepc_set_num_walk_threads,  FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_dipole_correction", pepc_set_dipole_correction, FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_load_balancing",    pepc_set_load_balancing,    FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_npm",               pepc_set_npm,               FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_debug_level",       pepc_set_debug_level,       FCS_PARSE_VAL(fcs_int));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}

FCSResult fcs_pepc_print_parameters(FCS handle)
{
  PEPC_CHECK_RETURN_RESULT(handle, __func__);

  printf("pepc epsilon: %" FCS_LMOD_FLOAT "f\n", handle->pepc_param->epsilon);
  printf("pepc theta: %" FCS_LMOD_FLOAT "f\n", handle->pepc_param->theta);
  printf("pepc num_walk_threads: %" FCS_LMOD_INT "d\n", handle->pepc_param->num_walk_threads);
  printf("pepc require virial: %" FCS_LMOD_INT "d\n", handle->pepc_param->require_virial);
  printf("pepc dipole correction: %" FCS_LMOD_INT "d\n", handle->pepc_param->dipole_correction);
  printf("pepc load balancing: %" FCS_LMOD_INT "d\n", handle->pepc_param->load_balancing);
  printf("pepc npm: %" FCS_LMOD_FLOAT "f\n", handle->pepc_param->npm);
  printf("pepc debug level: %" FCS_LMOD_INT "d\n", handle->pepc_param->debug_level);

  return FCS_RESULT_SUCCESS;  
}
