/*
  Copyright (C) 2011-2012 Rene Halver, Lukas Arnold, Mathias Winkel

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

#include "fcs_pepc.h"
#include "FCSCommon.h"
#include "mpi.h"
#include <string.h>


/* combined setter function for all pepc parameters */
extern FCSResult fcs_pepc_setup(FCS handle, fcs_float epsilon, fcs_float theta)
{
  /*	char* fnc_name = "fcs_pepc_setup"; */
  FCSResult result;

  result = fcs_pepc_set_epsilon(handle, epsilon);
  if (result != NULL)
    return result;


  result = fcs_pepc_set_theta(handle, theta);
  if (result != NULL)
    return result;

  return NULL;

}

/* method to check if pepc parameters are entered into checked FCS */
extern FCSResult fcs_pepc_check(FCS handle)
{
  char* fnc_name = "fcs_pepc_check";
  fcs_float *a,*b,*c;
  fcs_float eps, theta;
  FCSResult res;

  if ((res = fcs_pepc_get_epsilon(handle, &eps))) return res;
  if (eps == -1.0)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "pepc: epsilon not set");
  if ((res = fcs_pepc_get_theta(handle, &theta))) return res;
  if (theta == -1.0)
    return fcsResult_create(FCS_MISSING_ELEMENT, fnc_name, "pepc: theta not set");

  a = fcs_get_box_a(handle);
  b = fcs_get_box_b(handle);
  c = fcs_get_box_c(handle);
  if (!fcs_uses_principal_axes(a,b,c))
    printf("%s\n", "WARNING: support of pepc for non-cubic simulation boxes currently is experimental.");

  if (!fcs_get_near_field_flag(handle))
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, 
			    "pepc performs near field computations by itself!");

  return NULL;
}



/* initialization function for basic pepc parameters */
FCSResult fcs_pepc_init(FCS handle)
{
  /*	char* fnc_name = "fcs_pepc_init"; */
  FCSResult result;
  
  fcs_pepc_internal_t *pepc_internal;
  MPI_Comm comm  = fcs_get_communicator(handle);
  MPI_Fint fcomm = MPI_Comm_c2f(comm);

  pepc_scafacos_initialize(&fcomm);

  handle->method_context = malloc(sizeof(fcs_pepc_internal_t));
  pepc_internal = (fcs_pepc_internal_t*) handle->method_context;
  pepc_internal->work_length = -1;
  pepc_internal->work        = NULL;

  /* we only set some default values here */
  if ((result = fcs_pepc_set_theta(             handle,   0.6 ))) return result;
  if ((result = fcs_pepc_set_epsilon(           handle,   0.0 ))) return result;
  if ((result = fcs_pepc_require_virial(        handle,   0   ))) return result;
  if ((result = fcs_pepc_set_num_walk_threads(  handle,   3   ))) return result;
  if ((result = fcs_pepc_set_load_balancing(    handle,   0   ))) return result;
  if ((result = fcs_pepc_set_dipole_correction( handle,   1   ))) return result;
  if ((result = fcs_pepc_set_npm(               handle, -45.0 ))) return result;
  if ((result = fcs_pepc_set_debug_level(       handle,   0   ))) return result;

  return NULL;
}

/* internal pepc-specific tuning function */
FCSResult fcs_pepc_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles, fcs_float *positions, fcs_float *charges)
{
  /*	char* fnc_name = "fcs_pepc_tune"; */
  FCSResult result;

  /* currently, there is no tune functionality for pepc */
  result = fcs_pepc_check(handle);
  return result;
}

/* internal pepc-specific run function */
FCSResult fcs_pepc_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles, 
		       fcs_float *positions, fcs_float *charges, 
		       fcs_float *field, fcs_float *potentials)
{
  /*	char* fnc_name = "fcs_pepc_run"; */
  FCSResult result;
  fcs_int pcnt, requirevirial, db_level=0, num_walk_threads, load_balancing, lat_corr;
  fcs_int nparts_tot;
  fcs_float eps, theta, npm;
  fcs_pepc_internal_t *pepc_internal;

  nparts_tot = fcs_get_total_particles(handle);
  if ((result = fcs_pepc_get_epsilon(handle, &eps)))                        return result;
  if ((result = fcs_pepc_get_theta(handle, &theta)))                        return result;
  if ((result = fcs_pepc_get_num_walk_threads(handle, &num_walk_threads)))  return result;
  if ((result = fcs_pepc_get_load_balancing(handle, &load_balancing)))      return result;
  if ((result = fcs_pepc_get_dipole_correction(handle, &lat_corr)))         return result;
  if ((result = fcs_pepc_get_npm(handle, &npm)))                            return result;
  if ((result = fcs_pepc_get_debug_level(handle, &db_level)))               return result;
  requirevirial = handle->pepc_param->requirevirial;

  pepc_internal = (fcs_pepc_internal_t*) fcs_get_method_context(handle);
  
  if((local_particles != pepc_internal->work_length) || (load_balancing==0)){
    
    if(NULL != pepc_internal->work){
      free(pepc_internal->work);
    }
    
    pepc_internal->work = (fcs_float*)malloc(sizeof(fcs_float)*local_particles);
    pepc_internal->work_length = local_particles;

    for(int pcnt=0; pcnt<pepc_internal->work_length; pcnt++)
      pepc_internal->work[pcnt] = 1.0;
  }
  
  result = fcs_pepc_check(handle);
  if (result != NULL)
    return result;

  if (db_level > 3) {
    printf("*** run pepc kernel\n");
    printf("** local particles :       %" FCS_LMOD_INT "d\n", local_particles);
    printf("** local max particles :   %" FCS_LMOD_INT "d\n", local_max_particles);
    printf("** total particles :       %" FCS_LMOD_INT "d\n", nparts_tot);
    printf("** epsilon :               %f\n", eps);
    printf("** theta :                 %f\n", theta);
    printf("** npm :                   %f\n", npm);
    printf("** db_level :              %" FCS_LMOD_INT "d\n", db_level);
    printf("** require virial :        %" FCS_LMOD_INT "d\n", requirevirial);
    printf("** num walk threads :      %" FCS_LMOD_INT "d\n", num_walk_threads);
    printf("** dipole correction :     %" FCS_LMOD_INT "d\n", lat_corr);
    printf("** use load balancing :    %" FCS_LMOD_INT "d\n", load_balancing);
    printf("** size int   :            %d\n", (int)sizeof(fcs_int));
    printf("** size float :            %d\n", (int)sizeof(fcs_float));
    printf("** debug lattice pointers: %p\n", fcs_get_box_a(handle));
    printf("** debug lattice pointers: %p\n", fcs_get_box_b(handle));
    printf("** debug lattice pointers: %p\n", fcs_get_box_c(handle));

    printf("** debug output lattice in x %f %f %f\n", 
	   fcs_get_box_a(handle)[0], fcs_get_box_a(handle)[1], fcs_get_box_a(handle)[2]);
    printf("** debug output lattice in y %f %f %f\n", 
	   fcs_get_box_b(handle)[0], fcs_get_box_b(handle)[1], fcs_get_box_b(handle)[2]);
    printf("** debug output lattice in z %f %f %f\n", 
	   fcs_get_box_c(handle)[0], fcs_get_box_c(handle)[1], fcs_get_box_c(handle)[2]);
  }

  pepc_scafacos_run(&local_particles, &nparts_tot, positions, charges,
		    field, potentials, pepc_internal->work,
		    ((fcs_pepc_internal_t*)(handle->method_context))->virial,
		    fcs_get_box_a(handle), fcs_get_box_b(handle), fcs_get_box_c(handle),
		    fcs_get_periodicity(handle), &lat_corr,
		    &eps, &theta, &db_level, &num_walk_threads, &npm );

  if (db_level > 3)
    printf("virial(0,0) %12.4e\n", ((fcs_pepc_internal_t*)(handle->method_context))->virial[0]);

  if (db_level > 3)
    for(pcnt=0; pcnt<local_particles; pcnt++)
      printf("** e-field dump, particle %3" FCS_LMOD_INT "d, (ex,ey,ez): %f, %f, %f\n", 
	     pcnt, field[pcnt+0], field[pcnt+1], field[pcnt+2]);

  return NULL;
}

/* clean-up function for pepc */
extern FCSResult fcs_pepc_destroy(FCS handle)
{
  /*	char* fnc_name = "fcs_pepc_destroy"; */

  MPI_Comm comm  = fcs_get_communicator(handle);
  MPI_Fint fcomm = MPI_Comm_c2f(comm);
  pepc_scafacos_finalize(&fcomm);

  free(handle->method_context);

  if(NULL != ((fcs_pepc_internal_t*)fcs_get_method_context(handle))->work){
    free(((fcs_pepc_internal_t*)fcs_get_method_context(handle))->work);   
  }

  return NULL;
}

/******************************************************************************************************
 *
 *						Setter and Getter functions for pepc parameters
 *
 ******************************************************************************************************/


/* setter function for pepc parameter epsilon */
extern FCSResult fcs_pepc_set_epsilon(FCS handle, fcs_float epsilon)
{
  char* fnc_name = "fcs_pepc_set_epsilon";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name,
			    "null pointer supplied as handle");


  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name,
			    "wrong method chosen, please choose a method (method is not \"pepc\")");
  if (epsilon < 0 )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,
			    "epsilon > 0. has been violated");
  else
    {
      handle->pepc_param->epsilon = epsilon;
      return NULL;
    }
}

/* getter function for pepc parameter epsilon */
extern FCSResult fcs_pepc_get_epsilon(FCS handle, fcs_float* epsilon)
{
  char* fnc_name = "fcs_pepc_get_epsilon";

  if (!handle || !handle->pepc_param)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "received NULL pointer");

  *epsilon = handle->pepc_param->epsilon;
  return NULL;
}


/* setter function for pepc parameter theta */
extern FCSResult fcs_pepc_set_theta(FCS handle, fcs_float theta)
{
  char* fnc_name = "fcs_pepc_set_theta";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name,
			    "null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name,
			    "wrong method chosen, please choose a method (method is not \"pepc\")");
  if (theta < 0)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,
			    "0 <= theta has been violated");
  else
    {
      handle->pepc_param->theta = theta;
      return NULL;
    }
}

/* getter function for pepc parameter theta */
extern FCSResult fcs_pepc_get_theta(FCS handle, fcs_float* theta)
{
  char* fnc_name = "fcs_pepc_get_theta";

  if (!handle || !handle->pepc_param)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "received NULL pointer");

  *theta = handle->pepc_param->theta;
  return NULL;
}

/* getter function for pepc parameter num_walk_threads */
extern FCSResult fcs_pepc_get_num_walk_threads(FCS handle, fcs_int *num_walk_threads)
{
  char* fnc_name = "fcs_pepc_get_num_walk_threads";

  if (!handle || !handle->pepc_param)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "received NULL pointer");

  *num_walk_threads = handle->pepc_param->num_walk_threads;
  return NULL;
}

/* setter function for pepc parameter num_walk_threads */
extern FCSResult fcs_pepc_set_num_walk_threads(FCS handle, fcs_int num_walk_threads)
{
  char* fnc_name = "fcs_pepc_set_num_walk_threads";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name,
			    "null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name,
			    "wrong method chosen, please choose a method (method is not \"pepc\")");
  if (num_walk_threads < 1 )
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name,
			    "0 < num_walk_threads has been violated");
  else
    {
      handle->pepc_param->num_walk_threads = num_walk_threads;
      return NULL;
    }
}

FCSResult fcs_pepc_set_load_balancing(FCS handle, fcs_int load_balancing)
{
  char* fnc_name = "fcs_pepc_set_load_balancing";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name,
			    "null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name,
			    "wrong method chosen, please choose a method (method is not \"pepc\")");

  handle->pepc_param->load_balancing = load_balancing;
  return NULL;
}

FCSResult fcs_pepc_get_load_balancing(FCS handle, fcs_int* load_balancing)
{
  char* fnc_name = "fcs_pepc_get_load_balancing";

  if (!handle || !handle->pepc_param)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "received NULL pointer");

  *load_balancing = handle->pepc_param->load_balancing;
  return NULL;
}

FCSResult fcs_pepc_set_dipole_correction(FCS handle, fcs_int dipole_correction)
{
  char* fnc_name = "fcs_pepc_set_dipole_correction";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name,
			    "null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name,
			    "wrong method chosen, please choose a method (method is not \"pepc\")");

  handle->pepc_param->dipole_correction = dipole_correction;
  return NULL;
}

FCSResult fcs_pepc_get_dipole_correction(FCS handle, fcs_int* dipole_correction)
{
  char* fnc_name = "fcs_pepc_get_dipole_correction";

  if (!handle || !handle->pepc_param)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "received NULL pointer");

  *dipole_correction = handle->pepc_param->dipole_correction;
  return NULL;
}

/* setter function for pepc parameter npm */
extern FCSResult fcs_pepc_set_npm(FCS handle, fcs_float npm)
{
  char* fnc_name = "fcs_pepc_set_npm";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name,
			    "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name,
			    "wrong method chosen, please choose a method (method is not \"pepc\")");

  handle->pepc_param->npm = npm;
  return NULL;
}

/* getter function for pepc parameter npm */
extern FCSResult fcs_pepc_get_npm(FCS handle, fcs_float* npm)
{
  char* fnc_name = "fcs_pepc_get_npm";

  if (!handle || !handle->pepc_param)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "received NULL pointer");

  *npm = handle->pepc_param->npm;
  return NULL;
}

/* setter function for pepc parameter debug_level */
extern FCSResult fcs_pepc_set_debug_level(FCS handle, fcs_int level)
{
  char* fnc_name = "fcs_pepc_set_debuglevel";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name,
			    "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name,
			    "wrong method chosen, please choose a method (method is not \"pepc\")");

  handle->pepc_param->debug_level = level;
  return NULL;
}

/* getter function for pepc parameter debug_level */
extern FCSResult fcs_pepc_get_debug_level(FCS handle, fcs_int* level)
{
  char* fnc_name = "FCSResultfcs_pepc_get_debuglevel";

  if (!handle || !handle->pepc_param)
    return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "received NULL pointer");

  *level = handle->pepc_param->debug_level;
  return NULL;
}

/******************************************************************************************************
 *
 *						additional pepc functions
 *
 ******************************************************************************************************/

/* setter function to (de)activate virial computation in pepc */
extern FCSResult fcs_pepc_require_virial(FCS handle, fcs_int choice)
{
  char* fnc_name = "fcs_pepc_require_virial";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"pepc\")");
  if ((choice < 0) || (choice > 1))
    return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"require_virial must be set to 0 or 1");
  else
    {
      handle->pepc_param->requirevirial = choice;
      return NULL;
    }
}

/* function getting the virial result from pepc */
extern FCSResult fcs_pepc_get_virial(FCS handle, fcs_float virial[9])
{
  char* fnc_name = "fcs_pepc_get_virial";
  int i;

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (fcs_get_method(handle) != FCS_PEPC)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD,fnc_name,"wrong method chosen, please choose a method (method is not \"pepc\")");
  if (handle->pepc_param->requirevirial != 1)
    return fcsResult_create(FCS_LOGICAL_ERROR,fnc_name,"calculation of virial was not activated. Please call require_virial(fcs, 1) before calling run()");
  else
    {
      for (i=0;i<9;i++)
	virial[i] = ((fcs_pepc_internal_t*)(handle->method_context))->virial[i];

      return NULL;
    }
}

