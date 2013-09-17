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



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "FCSCommon.h"
#include "FCSInterface.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/******************************************************************************************************
 *
 *                     internal function definitions
 *
 ******************************************************************************************************/

/**
 * @brief function to check if initial information is set
 * @param handle FCS-object to be checked if obligatory
 * information is set, i.e. if fcs_init was run with this object
 * @return 0 if the information is not yet set, 1 if it is
 */
static fcs_int fcs_init_check(FCS handle)
{
    return (fcs_get_method(handle) != FCS_NO_METHOD_CHOSEN && fcs_get_communicator(handle) != MPI_COMM_NULL );
}

/**
 * @brief function to check if information inserted into fcs_tune
 * is set into the FCS-object (only non method related)
 * @param handle FCS-object to be checked if obligatory
 * information is set, i.e. if fcs_init was run with this object
 * @return 0 if the information is not yet set, 1 if it is
 */
static fcs_int fcs_tune_check(FCS handle)
{
    return   ( fcs_get_box_a(handle) && fcs_get_box_b(handle) && fcs_get_box_c(handle) &&
             (fcs_get_communicator(handle)!=MPI_COMM_NULL) && fcs_get_periodicity(handle) && (fcs_get_near_field_flag(handle)!=-1) &&
             (fcs_get_total_particles(handle)!=-1) );
}

/**
 * @brief function to check if the FCS-object is ready to
 * be passed to fcs_run
 * @param handle FCS-object to be checked if obligatory
 * information is set, i.e. if fcs_init was run with this object
 * @return 0 if the information is not yet set, 1 if it is
 */
static fcs_int fcs_run_check(FCS handle)
{
    fcs_int check = fcs_tune_check(handle);
    /*fcs_int method;*/
    if (check)
    {
        return check;
        /* commented out until the respective routines are ready */
        /*
        method = fcs_get_method(handle);
        if (strcmp_j(method, "pepc"))
        {
            return fcs_pepc_check(handle);
        }
        else if (strcmp_j(method, "fmm"))
        {
            return fcs_fmm_check(handle);
        }
        else if (strcmp_j(method, "p3m"))
        {
            return fcs_p3m_check(handle);
        }
        else if (strcmp_j(method, "p2nfft"))
        {
            return fcs_p2nfft_check(handle);
        }
        else if (strcmp_j(method, "pp3mg"))
        {
            return fcs_pp3mg_check(handle);
        }
        else if (strcmp_j(method, "direct"))
        {
            return fcs_direct_check(handle);
        }
        else
        {
            return 0;
        }
         */
    }
    else
        return check;
}


/******************************************************************************************************
 *
 *                         function definitions


 *
 ******************************************************************************************************/

/* init-method for ScaFaCoS library */
FCSResult fcs_init ( FCS *handle, const char* method, MPI_Comm communicator )
{
  char* fnc_name = "fcs_init";
  int i;

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  if (method == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as method name");

  *handle = (FCS)malloc(sizeof(FCS_t));
    
  if (*handle == NULL) 
    return fcsResult_create(FCS_ALLOC_FAILED,fnc_name,"allocation of FCS failed");



  /* Set all elements of the handle. The method is set below */
  (*handle)->communicator = communicator;
  (*handle)->near_field_flag = -1;
  for (i=0; i<3; ++i) {
    (*handle)->box_a[i] = 0.0;
    (*handle)->box_b[i] = 0.0;
    (*handle)->box_c[i] = 0.0;
    (*handle)->periodicity[i] = 0;
    (*handle)->offset[i] = 0.0;
  }
  (*handle)->total_particles = -1;
  (*handle)->total_charges = -1;
  (*handle)->method_context = NULL;

  /* Init all method parameter structs to NULL */


#ifdef FCS_ENABLE_DIRECT
  (*handle)->direct_param = NULL;
#endif
#ifdef FCS_ENABLE_FMM
  (*handle)->fmm_param = NULL;
#endif
#ifdef FCS_ENABLE_MEMD
  (*handle)->memd_param = NULL;
#endif
#ifdef FCS_ENABLE_P2NFFT
  (*handle)->p2nfft_param = NULL;
#endif
#ifdef FCS_ENABLE_PEPC
  (*handle)->pepc_param = NULL;
#endif
#ifdef FCS_ENABLE_PP3MG
  (*handle)->pp3mg_param = NULL;
#endif
#ifdef FCS_ENABLE_VMG
  (*handle)->vmg_param = NULL;
#endif
#ifdef FCS_ENABLE_WOLF
  (*handle)->wolf_param = NULL;
#endif

  (*handle)->set_max_particle_move = NULL;
  (*handle)->set_resort = NULL;
  (*handle)->get_resort = NULL;
  (*handle)->get_resort_availability = NULL;
  (*handle)->get_resort_particles = NULL;
  (*handle)->resort_ints = NULL;
  (*handle)->resort_floats = NULL;
  (*handle)->resort_bytes = NULL;

  /* Call the method-specific init functions */
#ifdef FCS_ENABLE_DIRECT
  if (strcmp(method, "direct") == 0) {
    (*handle)->method = FCS_DIRECT;
    (*handle)->direct_param = (fcs_direct_parameters) malloc(sizeof(fcs_direct_parameters_t));
    return fcs_direct_init(*handle);
  }
#endif

#ifdef FCS_ENABLE_EWALD
  if (strcmp(method, "ewald") == 0) {
    (*handle)->method = FCS_EWALD;
    return fcs_ewald_init(*handle, communicator);
  }
#endif

#ifdef FCS_ENABLE_FMM
  if (strcmp(method, "fmm") == 0) {
    (*handle)->method = FCS_FMM;
    (*handle)->fmm_param = (fcs_fmm_parameters)malloc(sizeof(fcs_fmm_parameters_t));
    /* setting fmm parameters to invalid values (or default values, if possible) */
    (*handle)->fmm_param->absrel = FCS_FMM_STANDARD_ERROR;
    (*handle)->fmm_param->tolerance_value = -1.0;
    (*handle)->fmm_param->dipole_correction = -1;
    (*handle)->fmm_param->potential = -1;
    (*handle)->fmm_param->cusp_radius = -1.0;
    return fcs_fmm_init(*handle);
  }
#endif

#ifdef FCS_ENABLE_MEMD
  if (strcmp(method, "memd") == 0) {
    (*handle)->method = FCS_MEMD;
    return fcs_memd_init(*handle, communicator);
  }
#endif

#ifdef FCS_ENABLE_MMM1D
  if (strcmp(method, "mmm1d") == 0) {
    (*handle)->method = FCS_MMM1D;
    return fcs_mmm1d_init(*handle, communicator);
  }
#endif

#ifdef FCS_ENABLE_MMM2D
  if (strcmp(method, "mmm2d") == 0) {
    (*handle)->method = FCS_MMM2D;
    return fcs_mmm2d_init(*handle, communicator);
  }
#endif

#ifdef FCS_ENABLE_P2NFFT
  if (strcmp(method, "p2nfft") == 0) {
    (*handle)->method = FCS_P2NFFT;
    return fcs_p2nfft_init(*handle);
  }
#endif

#ifdef FCS_ENABLE_P3M
  if (strcmp(method, "p3m") == 0) {
    (*handle)->method = FCS_P3M;
    return fcs_p3m_init(*handle, communicator);
  }
#endif

#ifdef FCS_ENABLE_PEPC
  if (strcmp(method, "pepc") == 0) {
    (*handle)->method = FCS_PEPC;
    (*handle)->pepc_param = (fcs_pepc_parameters)malloc(sizeof(fcs_pepc_parameters_t));
    (*handle)->pepc_param->theta = -1.0;
    (*handle)->pepc_param->epsilon = -1.0;
    (*handle)->pepc_param->num_walk_threads = 3;
    return fcs_pepc_init(*handle);
  }
#endif

#ifdef FCS_ENABLE_PP3MG
  if (strcmp(method, "pp3mg") == 0) {
    (*handle)->method = FCS_PP3MG;
    (*handle)->pp3mg_param = (fcs_pp3mg_parameters)malloc(sizeof(fcs_pp3mg_parameters_t));
    (*handle)->pp3mg_param->m = -1;
    (*handle)->pp3mg_param->n = -1;
    (*handle)->pp3mg_param->o = -1;
    (*handle)->pp3mg_param->ghosts = -1;
    (*handle)->pp3mg_param->max_particles = -1;
    (*handle)->pp3mg_param->degree = -1;
    (*handle)->pp3mg_param->maxiter = -1;
    (*handle)->pp3mg_param->tol = -1.0;
    (*handle)->pp3mg_param->distribution = 0;
    (*handle)->pp3mg_param->discretization = 0;
    return fcs_pp3mg_init(*handle);
  }
#endif

#ifdef FCS_ENABLE_VMG
  if (strcmp(method, "vmg") == 0) {
    (*handle)->method = FCS_VMG;
    (*handle)->vmg_param = (fcs_vmg_parameters)malloc(sizeof(fcs_vmg_parameters_t));
    (*handle)->vmg_param->max_level = -1;
    (*handle)->vmg_param->max_iterations = -1;
    (*handle)->vmg_param->smoothing_steps = -1;
    (*handle)->vmg_param->cycle_type = -1;
    (*handle)->vmg_param->precision = -1.0;
    (*handle)->vmg_param->near_field_cells = -1;
    (*handle)->vmg_param->interpolation_order = -1;
    (*handle)->vmg_param->discretization_order = -1;
    return fcs_vmg_init(*handle);
  }
#endif

#ifdef FCS_ENABLE_WOLF
  if (strcmp(method, "wolf") == 0) {
    (*handle)->method = FCS_WOLF;
    (*handle)->wolf_param = (fcs_wolf_parameters) malloc(sizeof(fcs_wolf_parameters_t));
    return fcs_wolf_init(*handle);
  }
#endif

  /* if we got here, method was chosen wrongly */
  return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
}



FCSResult fcs_destroy(FCS handle) {
  char* fnc_name = "fcs_destroy";
  FCSResult result;

  if (handle == NULL)
    return NULL;

  switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
    case FCS_DIRECT:
      result = fcs_direct_destroy(handle);
      free(handle->direct_param);
      if (result != NULL) return result;
      break;
#endif
#ifdef FCS_ENABLE_EWALD
    case FCS_EWALD:
      result = fcs_ewald_destroy(handle);
      if (result != NULL) return result;
      break;
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_FMM:
      result = fcs_fmm_destroy(handle);
      free(handle->fmm_param);
      if (result != NULL)
    return result;
      break;
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_MEMD:
      result = fcs_memd_destroy(handle);
      free(handle->memd_param);
      if (result != NULL)
    return result;
      break;
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_MMM1D:
      result = fcs_mmm1d_destroy(handle);
      if (result != NULL)
   return result;
      break;
#endif
#ifdef FCS_ENABLE_MMM2D
    case FCS_MMM2D:
      result = fcs_mmm2d_destroy(handle);
      if (result != NULL)
   return result;
      break;
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_PEPC:
      result = fcs_pepc_destroy(handle);
      free(handle->pepc_param);
      if (result != NULL)
    return result;
      break;
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_P2NFFT:
      result = fcs_p2nfft_destroy(handle);
      free(handle->p2nfft_param);
      if (result != NULL)
    return result;
      break;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_P3M:
      result = fcs_p3m_destroy(handle);
      if (result != NULL)
    return result;
      break;
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_PP3MG:
      result = fcs_pp3mg_destroy(handle);
      free(handle->pp3mg_param);
      if (result != NULL)
    return result;
      break;
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_VMG:
      result = fcs_vmg_destroy(handle);
      free(handle->vmg_param);
      if (result != NULL)
    return result;
      break;
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_WOLF:
      result = fcs_wolf_destroy(handle);
      free(handle->wolf_param);
      if (result != NULL) return result;
      break;
#endif
    default:
      return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
    }


  
  free(handle);

  return NULL;
}
  
/* tuning function to be called by the user */
FCSResult fcs_tune(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
           fcs_float *positions,  fcs_float *charges )
{
    char* fnc_name = "fcs_tune";
    FCSResult result;

    fcs_set_values_changed(handle,0);
    
    if (local_particles < 0)
      return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"number of local particles must be non negative");


    if (local_max_particles < 0)
      return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"maximum number of local particles must be non negative");
    if (local_particles > local_max_particles)
      return fcsResult_create(FCS_LOGICAL_ERROR,fnc_name,"maximum number of local particles should be greater or equal to current number of particles");

    
      if (fcs_init_check(handle) && fcs_tune_check(handle))
      {
        switch(fcs_get_method(handle))
          {
#ifdef FCS_ENABLE_DIRECT
            case FCS_DIRECT:
                result = fcs_direct_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_EWALD
            case FCS_EWALD:
                result = fcs_ewald_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_FMM
            case FCS_FMM:
                result = fcs_fmm_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_MEMD
            case FCS_MEMD:
                result = fcs_memd_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_MMM1D
            case FCS_MMM1D:
                result = fcs_mmm1d_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_MMM2D
            case FCS_MMM2D:
                result = fcs_mmm2d_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_P2NFFT
            case FCS_P2NFFT:
                result = fcs_p2nfft_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_P3M
            case FCS_P3M:
                result = fcs_p3m_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_PEPC
              case FCS_PEPC:
                result = fcs_pepc_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_PP3MG
            case FCS_PP3MG:
                result = fcs_pp3mg_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_VMG
            case FCS_VMG:
                result = fcs_vmg_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
#ifdef FCS_ENABLE_WOLF
            case FCS_WOLF:
                result = fcs_wolf_tune(handle, local_particles, local_max_particles, positions, charges);
                return result;
#endif
            default:
                return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
        }
    }

    return (fcsResult_create(FCS_MISSING_ELEMENT,fnc_name,"not all needed data has been inserted into the given handle"));
}



/* tuning function to be called by the user */
FCSResult fcs_run(FCS handle, fcs_int local_particles, fcs_int local_max_particles,
          fcs_float *positions, fcs_float *charges, 
          fcs_float *field, fcs_float *potentials )
{
    char* fnc_name = "fcs_run";
    FCSResult result;

    if (local_particles < 0)
        return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"number of local particles must be non negative");
    if (local_max_particles < 0)
        return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"maximum number of local particles must be non negative");
    if (local_particles > local_max_particles)
        return fcsResult_create(FCS_LOGICAL_ERROR,fnc_name,"maximum number of local particles should be greater or equal to current number of particles");

    if (fcs_get_values_changed(handle))
    {
        result = fcs_tune(handle, local_particles, local_max_particles, positions, charges);
        if (result != NULL) return result;
    }

    

    if (fcs_init_check(handle) && fcs_run_check(handle))
    {
        switch(fcs_get_method(handle))
        {
#ifdef FCS_ENABLE_DIRECT
            case FCS_DIRECT:
                result = fcs_direct_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_EWALD
            case FCS_EWALD:
                result = fcs_ewald_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_FMM
            case FCS_FMM:
                result = fcs_fmm_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_MEMD
            case FCS_MEMD:
                result = fcs_memd_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_MMM1D
            case FCS_MMM1D:
                result = fcs_mmm1d_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_MMM2D
            case FCS_MMM2D:
                result = fcs_mmm2d_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_PEPC
            case FCS_PEPC:
                result = fcs_pepc_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_P2NFFT
            case FCS_P2NFFT:
                result = fcs_p2nfft_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_P3M
            case FCS_P3M:
                result = fcs_p3m_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_PP3MG
            case FCS_PP3MG:
                result = fcs_pp3mg_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_VMG
            case FCS_VMG:
                result = fcs_vmg_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
#ifdef FCS_ENABLE_WOLF
            case FCS_WOLF:
                result = fcs_wolf_run(handle, local_particles, local_max_particles, positions, charges, field, potentials);
                return result;
#endif
            default:
                return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
        }
    }

    return (fcsResult_create(FCS_MISSING_ELEMENT,fnc_name,"not all needed data has been inserted into the given handle"));
}


FCSResult fcs_method_has_near(FCS handle, fcs_int *has_near)
{
    const char* fnc_name = "fcs_method_has_near";

    if (fcs_init_check(handle) && fcs_tune_check(handle))
    {
        switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
            case FCS_DIRECT:
                *has_near = 0;
                return NULL;
#endif
#ifdef FCS_ENABLE_EWALD
            case FCS_EWALD:
                *has_near = 0;
                return NULL;
#endif
#ifdef FCS_ENABLE_FMM
            case FCS_FMM:
                *has_near = 0;
                return NULL;
#endif
#ifdef FCS_ENABLE_MEMD
            case FCS_MEMD:
                *has_near = 0;
                return NULL;
#endif
#ifdef FCS_ENABLE_MMM1D
            case FCS_MMM1D:
                *has_near = 0; 
                return NULL;
#endif
#ifdef FCS_ENABLE_MMM2D
            case FCS_MMM2D:
                *has_near = 0; 
                return NULL;
#endif
#ifdef FCS_ENABLE_PEPC
            case FCS_PEPC:
                *has_near = 0;
                return NULL;
#endif
#ifdef FCS_ENABLE_P2NFFT
            case FCS_P2NFFT:
                *has_near = 1;
                return NULL;
#endif
#ifdef FCS_ENABLE_P3M
            case FCS_P3M:
                *has_near = 1; 
                return NULL;
#endif
#ifdef FCS_ENABLE_PP3MG
            case FCS_PP3MG:
                *has_near = 0;
                return NULL;
#endif
#ifdef FCS_ENABLE_VMG
            case FCS_VMG:
                *has_near = 0;
                return NULL;
#endif
            default:
                return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
        }
    }
    else
        return (fcsResult_create(FCS_MISSING_ELEMENT,fnc_name,"not all needed data has been inserted into the given handle"));
}

FCSResult fcs_compute_near_potential(FCS handle, fcs_float dist, fcs_float* potential) 
{
    const char* fnc_name = "fcs_compute_near_potential";

    if (fcs_init_check(handle) && fcs_tune_check(handle))
    {
        switch(fcs_get_method(handle))
        {
#ifdef FCS_ENABLE_P2NFFT
            case FCS_P2NFFT:
                *potential = fcs_p2nfft_compute_near_potential(handle, dist);
                return NULL;
#endif
#ifdef FCS_ENABLE_P3M
            case FCS_P3M:
            {
                fcs_p3m_near_parameters_t params;
                fcs_p3m_get_near_parameters(handle, &params);
                *potential = fcs_p3m_compute_near_potential(params, dist);
                return NULL;
            }
#endif
            default:
                return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
        }
    }
    else
        return (fcsResult_create(FCS_MISSING_ELEMENT,fnc_name,"not all needed data has been inserted into the given handle"));
}

FCSResult fcs_compute_near_field(FCS handle, fcs_float dist, fcs_float* field)
{
    const char* fnc_name = "fcs_compute_near_field";

    if (fcs_init_check(handle) && fcs_tune_check(handle))
    {
        switch(fcs_get_method(handle))
        {
#ifdef FCS_ENABLE_P2NFFT
            case FCS_P2NFFT:
                *field = fcs_p2nfft_compute_near_field(handle, dist);
                return NULL;
#endif
#ifdef FCS_ENABLE_P3M
            case FCS_P3M:
            {
                fcs_p3m_near_parameters_t params;
                fcs_p3m_get_near_parameters(handle, &params);
                *field = fcs_p3m_compute_near_field(params, dist);
                return NULL;
            }
#endif
            default:
                return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
        }
    }
    else
        return (fcsResult_create(FCS_MISSING_ELEMENT,fnc_name,"not all needed data has been inserted into the given handle"));
}

FCSResult fcs_compute_near(FCS handle, fcs_float dist, fcs_float* potential, fcs_float* field)
{
    const char* fnc_name = "fcs_compute_near";

    if (fcs_init_check(handle) && fcs_tune_check(handle))
    {
        switch(fcs_get_method(handle))
        {
#ifdef FCS_ENABLE_P2NFFT
            case FCS_P2NFFT:
                fcs_p2nfft_compute_near(handle, dist, potential, field);
                return NULL;
#endif
#ifdef FCS_ENABLE_P3M
            case FCS_P3M:
            {
                fcs_p3m_near_parameters_t params;
                fcs_p3m_get_near_parameters(handle, &params);
                fcs_p3m_compute_near(params, dist, potential, field);
                return NULL;
            }
#endif
            default:
                return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
        }
    }
    else
        return (fcsResult_create(FCS_MISSING_ELEMENT,fnc_name,"not all needed data has been inserted into the given handle"));
}

FCSResult fcs_compute_dipole_correction(FCS handle,
                    fcs_int local_particles,
                    fcs_float* positions, fcs_float *charges, 
                    fcs_float epsilon,
                    fcs_float *field_correction,
                    fcs_float *energy_correction)
{
  /* Local dipole moment */
  fcs_float local_dipole_moment[3] = {0.0, 0.0, 0.0};
  /* Global dipole moment */
  fcs_float dipole_moment[3];

  fcs_int pid;
  fcs_int dim;

  if (fcs_float_is_zero(epsilon) || epsilon > 0.0) {
    /* Compute the global dipole moment */
    for (pid = 0; pid < local_particles; pid++)
      for (dim = 0; dim < 3; dim++)
        local_dipole_moment[dim] += charges[pid]*positions[pid*3+dim];
    MPI_Allreduce(local_dipole_moment, dipole_moment, 3, FCS_MPI_FLOAT,
      MPI_SUM, handle->communicator);

    fcs_float *a = fcs_get_box_a(handle);
    fcs_float *b = fcs_get_box_b(handle);
    fcs_float *c = fcs_get_box_c(handle);
    
    /* Volume of the parallelepiped */
    fcs_float volume = 
        a[0] * (b[1]*c[2] - b[2]*c[1]) 
      + a[1] * (b[2]*c[0] - b[0]*c[2])
      + a[2] * (b[0]*c[1] - b[1]*c[0]);

    fcs_float pref = 4.0*3.14159265358979323846264338328 
      / (3.0*volume*(epsilon + 1.0));

    if (energy_correction)
      *energy_correction = 0.5*pref*(dipole_moment[0]*dipole_moment[0]
        + dipole_moment[1]*dipole_moment[1]
        + dipole_moment[2]*dipole_moment[2]);
    
    if (field_correction) {
      field_correction[0] = -pref*dipole_moment[0];
      field_correction[1] = -pref*dipole_moment[1];
      field_correction[2] = -pref*dipole_moment[2];
    }
  } else {
    /* metallic BC (epsilon=+infty) */
    if (energy_correction)
      *energy_correction = 0.0;
    if (field_correction) {
      field_correction[0] = 0.0;
      field_correction[1] = 0.0;
      field_correction[2] = 0.0;
    }
  }

  return NULL;
}

/* getter function for the method */
fcs_int fcs_get_method(FCS handle)
{
    char* fnc_name = "fcs_get_method";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return NULL", fnc_name);
        return FCS_NO_METHOD_CHOSEN;
    }
    return handle->method;
}

/* setter function for the dimension */
FCSResult fcs_set_dimension(FCS handle, fcs_int dim)
{
    char* fnc_name = "fcs_set_dimension";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (dim > 3 || dim < 1)
        return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"dimension must be between 1 and 3");
    else
    {
        handle->dimension = dim;
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the dimension */
fcs_int fcs_get_dimension(FCS handle)
{
    char* fnc_name = "fcs_get_dimension";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return -1", fnc_name);
        return -1;
    }
    return handle->dimension;
}

/* setter function for the near field flag */
FCSResult fcs_set_near_field_flag(FCS handle, fcs_int srf)
{
    char* fnc_name = "fcs_set_near_field_flag";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    else
    {
        handle->near_field_flag = srf;
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the near field flag */
fcs_int fcs_get_near_field_flag(FCS handle)
{
    char* fnc_name = "fcs_get_near_field_flag";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return -1", fnc_name);
        return -1;
    }
    return handle->near_field_flag;
}

/* getter function for the communicator */
MPI_Comm fcs_get_communicator(FCS handle)
{
    char* fnc_name = "fcs_get_communicator";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return MPI_COMM_NULL", fnc_name);
        return MPI_COMM_NULL;
    }
    return handle->communicator;
}

/* setter function for the first box base vector */
FCSResult fcs_set_box_a(FCS handle, fcs_float box_a[])
{
    char* fnc_name = "fcs_set_box_a";
    int i;

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (box_a == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as vector");
    else
    {
        for (i = 0; i < 3; ++i)
            handle->box_a[i] = box_a[i];
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the first box base vector */
fcs_float* fcs_get_box_a(FCS handle)
{
    char* fnc_name = "fcs_get_box_a";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return NULL", fnc_name);
        return NULL;
    }
    return handle->box_a;
}


/* setter function for the second box base vector */
FCSResult fcs_set_box_b(FCS handle, fcs_float box_b[])
{
    char* fnc_name = "fcs_set_box_b";
    int i;

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (box_b == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as vector");
    else
    {
        for (i = 0; i < 3; ++i)
            handle->box_b[i] = box_b[i];
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the second box base vector */
fcs_float* fcs_get_box_b(FCS handle)
{
    char* fnc_name = "fcs_get_box_b";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return NULL", fnc_name);
        return NULL;
    }
    return handle->box_b;
}



/* setter function for the third box base vector */
FCSResult fcs_set_box_c(FCS handle, fcs_float box_c[])
{
    char* fnc_name = "fcs_set_box_c";
    int i;

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (box_c == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as vector");
    else
    {
        for (i = 0; i < 3; ++i)
          handle->box_c[i] = box_c[i];
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the first box base vector */


fcs_float* fcs_get_box_c(FCS handle)
{
    char* fnc_name = "fcs_get_box_c";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return NULL", fnc_name);
        return NULL;
    }
    return handle->box_c;

}

/* setter function for periodicity of the box */
FCSResult fcs_set_periodicity(FCS handle, fcs_int* periodicity)
{
    char* fnc_name = "fcs_set_periodicity";
    int i;
    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (periodicity == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as vector");
    else
    {
        for (i = 0; i < 3; ++i)
            *(handle->periodicity+i) = *(periodicity+i);
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the periodicity of the box */
fcs_int* fcs_get_periodicity(FCS handle)
{
    char* fnc_name = "fcs_get_periodicity";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return NULL", fnc_name);
        return NULL;
    }
    return handle->periodicity;
}

/* setter function for periodicity of the box */
FCSResult fcs_set_offset(FCS handle, fcs_float* offset)
{
    char* fnc_name = "fcs_set_offset";
    int i;
    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (offset == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as vector");
    else
    {
        for (i = 0; i < 3; ++i)
            *(handle->offset+i) = *(offset+i);
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the offset of the box */
fcs_float* fcs_get_offset(FCS handle)
{
    char* fnc_name = "fcs_get_offset";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return NULL", fnc_name);
        return NULL;
    }
    return handle->offset;
}
/* setter function for the total particle count */
FCSResult fcs_set_total_particles(FCS handle, fcs_int total_particles)
{
    char* fnc_name = "fcs_set_total_particles";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (total_particles < 1)
        return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"total particles must be at least 1");
    else
    {
        handle->total_particles = total_particles;
        fcs_set_values_changed(handle,1);
        return NULL;
    }
}

/* getter function for the total particle count */
fcs_int fcs_get_total_particles(FCS handle)
{
    char* fnc_name = "fcs_get_total_particles";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return -1", fnc_name);
        return 0;
    }
    return handle->total_particles;
}

/* function to fill in all obligatory data at once */
FCSResult fcs_set_common(FCS handle, fcs_int near_field_flag, fcs_float* box_a, fcs_float* box_b, fcs_float* box_c, fcs_float* offset, fcs_int* periodicity, fcs_int total_particles)
{
    /*char* fnc_name = "fcs_set_common";*/
    FCSResult result;

    result = fcs_set_near_field_flag(handle,near_field_flag);
    if (result != NULL)
        return    result;
    result = fcs_set_box_a(handle,box_a);
    if (result != NULL)
        return    result;
    result = fcs_set_box_b(handle,box_b);
    if (result != NULL)
        return    result;
    result = fcs_set_box_c(handle,box_c);
    if (result != NULL)
        return    result;
    result = fcs_set_offset(handle,offset);
    if (result != NULL)
        return  result;
    result = fcs_set_periodicity(handle, periodicity);
    if (result != NULL)
        return    result;
    result = fcs_set_total_particles(handle,total_particles);
    if (result != NULL)
        return    result;
    return NULL;
}

/* setter function for the method context */


FCSResult fcs_set_method_context(FCS handle, void *method_context)
{
    char* fnc_name = "fcs_set_method_context";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");

    handle->method_context = method_context;

    return NULL;
}

/* getter function for the method context */
void *fcs_get_method_context(FCS handle)
{
    char* fnc_name = "fcs_get_method_context";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
  return handle->method_context;
}

/* output function */
void fcs_printHandle(FCS handle)
{
  if (handle)
    {
      switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
    case FCS_DIRECT:
      printf("chosen method: direct solver\n");
      break;
#endif
#ifdef FCS_ENABLE_EWALD
    case FCS_EWALD:
      printf("chosen method: ewald\n");
      break;
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_FMM:
      printf("chosen method: fmm\n");
      break;
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_MEMD:
      printf("chosen method: memd\n");
      break;
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_MMM1D:
      printf("chosen method: mmm1d\n");
      break;
#endif
#ifdef FCS_ENABLE_MMM2D
    case FCS_MMM2D:
      printf("chosen method: mmm2d\n");
      break;
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_PEPC:
      printf("chosen method: pepc\n");
      break;
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_P2NFFT:
      printf("chosen method: p2nfft\n");
      break;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_P3M:
      printf("chosen method: p3m\n");
      break;
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_PP3MG:
      printf("chosen method: pp3mg\n");
      break;
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_VMG:
      printf("chosen method: vmg\n");
      break;
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_WOLF:
      printf("chosen method: wolf solver\n");
      break;
#endif
    default:
      printf("chosen method: UNKNOWN\n");
      break;
    }
      printf("near field computations done by solver: %c\n", 
         (fcs_get_near_field_flag(handle)?'T':'F'));

      printf("box vectors: [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f], [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f]\n",
         (fcs_float)fcs_get_box_a(handle)[0],(fcs_float)fcs_get_box_a(handle)[1],(fcs_float)fcs_get_box_a(handle)[2],
         (fcs_float)fcs_get_box_b(handle)[0],(fcs_float)fcs_get_box_b(handle)[1],(fcs_float)fcs_get_box_b(handle)[2],
         (fcs_float)fcs_get_box_c(handle)[0],(fcs_float)fcs_get_box_c(handle)[1],(fcs_float)fcs_get_box_c(handle)[2]);
      printf("offset: [%10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f %10.4" FCS_LMOD_FLOAT "f]\n", (fcs_float)fcs_get_offset(handle)[0],
         (fcs_float)fcs_get_offset(handle)[1],(fcs_float)fcs_get_offset(handle)[2]);
      printf("periodicity: %c %c %c\n", ((fcs_get_periodicity(handle)[0] == 1)?'T':'F'), ((fcs_get_periodicity(handle)[1] == 1)?'T':'F'),((fcs_get_periodicity(handle)[2] == 1)?'T':'F') );
      printf("total particles: %" FCS_LMOD_INT "d\n", fcs_get_total_particles(handle));
      printf("------------------------");
      printf("solver specific data:\n");

      switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
    case FCS_DIRECT:
      {
        fcs_float cutoff;
        fcs_int images[3];
        FCSResult res;
        if ((res = fcs_direct_get_cutoff(handle, &cutoff))) fcsResult_printResult(res);
        printf("direct cutoff: %" FCS_LMOD_FLOAT "e\n", cutoff);
        if ((res = fcs_direct_get_periodic_images(handle, images))) fcsResult_printResult(res);
        printf("direct periodic images: %5" FCS_LMOD_INT "d %5" FCS_LMOD_INT "d %5" FCS_LMOD_INT "d\n", images[0], images[1], images[2]);
        break;
      }
#endif
#ifdef FCS_ENABLE_EWALD
    case FCS_EWALD: {
      FCSResult res = fcs_ewald_print_parameters(handle);
      if (res) 
        fcsResult_printResult(res);
      break;
    }
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_FMM: {
      fcs_int absrel;
      fcs_float tolerance_value;
      fcs_int dipole_correction;
      long long tuning;
      long long maxdepth;
      long long limit;
      long long load;
      fcs_fmm_get_absrel(handle, &absrel);
      fcs_fmm_get_tolerance_energy(handle, &tolerance_value);
      fcs_fmm_get_dipole_correction(handle, &dipole_correction);
      fcs_fmm_get_internal_tuning(handle, &tuning);
      fcs_fmm_get_balanceload(handle, &load);
      fcs_fmm_get_maxdepth(handle, &maxdepth);
      fcs_fmm_get_unroll_limit(handle, &limit);
      printf("fmm absrel: %" FCS_LMOD_INT "d\n", absrel);
      printf("fmm tolerance value: %e\n", tolerance_value);
      printf("fmm dipole correction: %" FCS_LMOD_INT "d\n", dipole_correction);
      printf("fmm internal tuning: %c\n", (tuning)?'T':'F');
      printf("fmm maxdepth: %lld\n", maxdepth);
      printf("fmm unroll limit: %lld\n", limit);
      printf("fmm internal balance load: %c\n", (load)?'T':'F');
      break;
    }
#endif
#ifdef FCS_ENABLE_MEMD
    case FCS_MEMD:
      break;
#endif
#ifdef FCS_ENABLE_MMM1D
    case FCS_MMM1D: {
      fcs_float radius;
      fcs_float PWerror;
      fcs_int cutoff;
      fcs_mmm1d_get_far_switch_radius(handle, &radius);
      fcs_mmm1d_get_bessel_cutoff(handle, &cutoff);
      fcs_mmm1d_get_maxPWerror(handle, &PWerror);
      printf("mmm1d bessel cutoff: %" FCS_LMOD_INT "d\n", cutoff);
      printf("mmm1d far switch radius: %e\n", radius);
      printf("mmm1d maximum PWerror: %e\n", PWerror);
      break;
    }
#endif
#ifdef FCS_ENABLE_MMM2D
  case FCS_MMM2D: {
      fcs_float contrasts_min, contrasts_max;
      fcs_float PWerror;
      fcs_float cutoff;
      fcs_int layers;
      fcs_float skin;
      fcs_mmm2d_get_dielectric_contrasts(handle, &contrasts_min, &contrasts_max);
      fcs_mmm2d_get_far_cutoff(handle, &cutoff);
      fcs_mmm2d_get_layers_per_node(handle, &layers);
      fcs_mmm2d_get_maxPWerror(handle, &PWerror);
      fcs_mmm2d_get_skin(handle, &skin);
      printf("mmm2d dielectric contrasts: %e %e\n", contrasts_min, contrasts_max);
      printf("mmm2d far cutoff: %e\n", cutoff);
      printf("mmm2d layer per node: %" FCS_LMOD_INT "d\n", layers);
      printf("mmm2d maximum PWerror: %e\n", PWerror);
      printf("mmm2d skin: %e\n", skin);
      break;
    }
#endif
#ifdef FCS_ENABLE_PEPC
    case FCS_PEPC:
      {
        fcs_float theta, eps;
        fcs_int num_walk_threads;
        FCSResult res;
        if ((res = fcs_pepc_get_theta(handle, &theta)))                       fcsResult_printResult(res);
        if ((res = fcs_pepc_get_epsilon(handle, &eps)))                       fcsResult_printResult(res);
        if ((res = fcs_pepc_get_num_walk_threads(handle, &num_walk_threads))) fcsResult_printResult(res);
        printf("pepc theta: %e\n", theta);
        printf("pepc epsilon: %e\n", eps);
        printf("pepc num_walk_threads: %" FCS_LMOD_INT "d\n", num_walk_threads);
        printf("pepc user requires virial: %4" FCS_LMOD_INT "d\n", handle->pepc_param->requirevirial);
        break;
      }
#endif
#ifdef FCS_ENABLE_P2NFFT
    case FCS_P2NFFT:
      {
        fcs_int tolerance_type;
        fcs_float tolerance;
        fcs_get_tolerance(handle, &tolerance_type, &tolerance);
        if(tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
          printf("p2nfft: tolerance_type = FCS_TOLERANCE_TYPE_FIELD, tolerance = %" FCS_LMOD_FLOAT "e\n", tolerance);
        else if(tolerance_type == FCS_TOLERANCE_TYPE_POTENTIAL)
          printf("p2nfft: tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL, tolerance = %" FCS_LMOD_FLOAT "e\n", tolerance);
        break;
      }      
      break;
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_P3M:
      {
        fcs_float tolerance;
        fcs_p3m_get_tolerance_field(handle, &tolerance);
        printf("p3m absolute field tolerance: %" FCS_LMOD_FLOAT "e\n", tolerance);
        break;
      }
#endif
#ifdef FCS_ENABLE_PP3MG
    case FCS_PP3MG:
      {
	fcs_int cells_x, cells_y, cells_z;
	fcs_int ghosts;
	fcs_int degree;
	fcs_int max_particles;
	fcs_int max_iterations;
	fcs_float tol;
	fcs_int distribution;
	fcs_int discretization;

	fcs_pp3mg_get_cells_x(handle, &cells_x);
	fcs_pp3mg_get_cells_y(handle, &cells_y);
	fcs_pp3mg_get_cells_z(handle, &cells_z);
	fcs_pp3mg_get_ghosts(handle, &ghosts);
	fcs_pp3mg_get_degree(handle, &degree);
	fcs_pp3mg_get_max_particles(handle, &max_particles);
	fcs_pp3mg_get_max_iterations(handle, &max_iterations);
	fcs_pp3mg_get_tol(handle, &tol);
	fcs_pp3mg_get_distribution(handle, &distribution);
	fcs_pp3mg_get_discretization(handle, &discretization);
 
	printf("pp3mg cells x: %" FCS_LMOD_INT "d\n",cells_x);
	printf("pp3mg cells y: %" FCS_LMOD_INT "d\n",cells_y);
	printf("pp3mg cells z: %" FCS_LMOD_INT "d\n",cells_z);
	printf("pp3mg ghosts: %" FCS_LMOD_INT "d\n",ghosts);
	printf("pp3mg degree: %" FCS_LMOD_INT "d\n",degree);
	printf("pp3mg max_particles: %" FCS_LMOD_INT "d\n",max_particles);
	printf("pp3mg max_iterations: %" FCS_LMOD_INT "d\n",max_iterations);
	printf("pp3mg tol: %e\n",tol);
	printf("pp3mg distribution: %" FCS_LMOD_INT "d\n",distribution);
	printf("pp3mg discretization: %" FCS_LMOD_INT "d\n",discretization);
      }
      break;
#endif
#ifdef FCS_ENABLE_VMG
    case FCS_VMG:
    {
      fcs_int level;
      fcs_int max_iter;
      fcs_int smoothing_steps;
      fcs_int cycle_type;
      fcs_float precision;
      fcs_int near_field_cells;
      fcs_int interpolation_order;
      fcs_int discretization_order;

      fcs_vmg_get_max_level(handle, &level);
      fcs_vmg_get_max_iterations(handle, &max_iter);
      fcs_vmg_get_smoothing_steps(handle, &smoothing_steps);
      fcs_vmg_get_cycle_type(handle, &cycle_type);
      fcs_vmg_get_precision(handle, &precision);
      fcs_vmg_get_near_field_cells(handle, &near_field_cells);
      fcs_vmg_get_interpolation_order(handle, &interpolation_order);
      fcs_vmg_get_discretization_order(handle, &discretization_order);

      printf("vmg max level:            %" FCS_LMOD_INT "d\n", level);
      printf("vmg max iterations:       %" FCS_LMOD_INT "d\n", max_iter);
      printf("vmg smoothing steps:      %" FCS_LMOD_INT "d\n", smoothing_steps);
      printf("vmg cycle_type:           %" FCS_LMOD_INT "d\n", cycle_type);
      printf("vmg precision:            %e\n", precision);
      printf("vmg near field cells:     %" FCS_LMOD_INT "d\n", near_field_cells);
      printf("vmg interpolation degree: %" FCS_LMOD_INT "d\n", interpolation_order);
      printf("vmg discretization order: %" FCS_LMOD_INT "d\n", discretization_order);
    }
      break;
#endif
#ifdef FCS_ENABLE_WOLF
    case FCS_WOLF:
      break;
#endif
    default:
      printf("no solver specific data found or unknown solver specified\n");
      break;
    }
    }
  else
    printf("passed handle is NULL pointer, no additional output!");
}


#if defined(FCS_ENABLE_DEBUG) || 0
# define PRINT_PARAM_BEGIN(_f_)     printf("%s: calling " #_f_ "(handle", func_name)
# define PRINT_PARAM_VAL(_f_, _v_)  printf( _f_, _v_)
# define PRINT_PARAM_STR(_str_)     printf("%s", (_str_))
# define PRINT_PARAM_END(_r_)       printf(") -> %p\n", _r_)
#else
# define PRINT_PARAM_BEGIN(_f_)     do {} while (0)
# define PRINT_PARAM_VAL(_f_, _v_)  do {} while (0)
# define PRINT_PARAM_STR(_str_)     do {} while (0)
# define PRINT_PARAM_END(_r_)       do {} while (0)
#endif

typedef long long fcs_long_long_t;
typedef char *fcs_p_char_t;

#define MAKE_TYPE_FUNC(_type_, _atox_, _format_) \
  static inline _type_ *parse_##_type_(char **s, _type_ *v) { \
    if (v == NULL) return NULL; \
    *v = (_type_) _atox_(*s); \
    *s = strchr(*s, ','); \
    if (*s) { **s = 0; *s += 1; } \
    PRINT_PARAM_VAL(_format_, *v); \
    return v; \
  } \
  static inline _type_ *const_##_type_(_type_ c, _type_ *v) { \
    if (v == NULL) return NULL; \
    *v = c; \
    PRINT_PARAM_VAL(_format_, *v); \
    return v; \
  }

static inline fcs_bool atob(const char *nptr)
{
  const char false_str[] = "false";
  if ((strlen(nptr) == 1 && strncmp(nptr, "0", 1) == 0) || (strlen(nptr) == strlen(false_str) && strncasecmp(nptr, false_str, strlen(false_str)))) return FCS_FALSE;
  return FCS_TRUE;
}

MAKE_TYPE_FUNC(fcs_int, atoll, "%" FCS_LMOD_INT "d")
MAKE_TYPE_FUNC(fcs_float, atof, "%" FCS_LMOD_FLOAT "f")
MAKE_TYPE_FUNC(fcs_bool, atob, "%" FCS_LMOD_INT "d")
MAKE_TYPE_FUNC(fcs_long_long_t, atoll, "%lld")
MAKE_TYPE_FUNC(fcs_p_char_t, , "%s")

#define PARSE_SEQ_MAX  3

#define PARAM_SELECTED(_str_, _param_) \
  (strcmp(param, #_param_) == 0 || strcmp(param, _str_) == 0)

#define IF_PARAM_INTRO(_str_, _param_) \
  if (PARAM_SELECTED(_str_, _param_)) { \
    FCSResult _r; \
    struct { \
      void *t; \
      fcs_int v_fcs_int[PARSE_SEQ_MAX]; \
      fcs_float v_fcs_float[PARSE_SEQ_MAX]; \
      fcs_bool v_fcs_bool[PARSE_SEQ_MAX]; \
      fcs_long_long_t v_fcs_long_long_t[PARSE_SEQ_MAX]; \
      fcs_p_char_t v_fcs_p_char_t[PARSE_SEQ_MAX]; \
    } _t; \
    char *_n=NULL; \
    PRINT_PARAM_BEGIN(_param_);

#define IF_PARAM_EXTRO() \
    PRINT_PARAM_END(_r); \
    if (_r != FCS_RESULT_SUCCESS && FCS_IS_FALSE(continue_on_errors)) return _r; \
    goto next_param; \
  }

#define PARSE_VAL(_type_) \
  _t.v_##_type_; \
  if (cur) { \
    PRINT_PARAM_STR(", "); \
    parse_##_type_(&cur, &_t.v_##_type_[0]); \
    _n = cur; cur = NULL; \
  } _type_

#define PARSE_SEQ(_type_, _n_) \
  &_t.v_##_type_; \
  if (cur) { \
    PRINT_PARAM_STR(", ["); \
    for (int _i = 0; _i < (_n_) && _i < PARSE_SEQ_MAX; ++_i) { \
      PRINT_PARAM_STR((_i == 0)?"":", "); \
      parse_##_type_(&cur, &_t.v_##_type_[_i]); \
    } \
    _n = cur; cur = NULL; \
    PRINT_PARAM_STR("]"); \
  } _type_ *

#define CONST_VAL(_type_, _c_) \
  _t.v_##_type_; \
  if (cur) { \
    PRINT_PARAM_STR(", "); \
    const_##_type_(_c_, &_t.v_##_type_[0]); \
    _n = cur; cur = NULL; \
  } _type_

#define IF_PARAM_THEN_FUNC0_GOTO_NEXT(_str_, _param_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _r = fcs_##_param_(handle); \
  IF_PARAM_EXTRO()

#define IF_PARAM_THEN_FUNC1_GOTO_NEXT(_str_, _param_, _p0_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _r = fcs_##_param_(handle, _vv0); \
  IF_PARAM_EXTRO()

#define IF_PARAM_THEN_FUNC2_GOTO_NEXT(_str_, _param_, _p0_, _p1_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _t.t = _p1_ _v1 = *_p1_ _vv1 = _v1; cur = _n; \
    _r = fcs_##_param_(handle, _vv0, _vv1); \
  IF_PARAM_EXTRO()

#define IF_PARAM_THEN_FUNC3_GOTO_NEXT(_str_, _param_, _p0_, _p1_, _p2_) \
  IF_PARAM_INTRO(_str_, _param_) \
    _t.t = _p0_ _v0 = *_p0_ _vv0 = _v0; cur = _n; \
    _t.t = _p1_ _v1 = *_p1_ _vv1 = _v1; cur = _n; \
    _t.t = _p2_ _v2 = *_p2_ _vv2 = _v2; cur = _n; \
    _r = fcs_##_param_(handle, _vv0, _vv1, _vv2); \
  IF_PARAM_EXTRO()

#define DUMMY_REFERENCE_TO_STATIC_FUNCTIONS(_str_) \
  if (PARAM_SELECTED(_str_, "EVEN_IF_IT_MATCHES_IT_DOES_NOTHING")) { \
    parse_fcs_int(NULL, NULL); \
    const_fcs_int(0, NULL); \
    parse_fcs_float(NULL, NULL); \
    const_fcs_float(0, NULL); \
    parse_fcs_bool(NULL, NULL); \
    const_fcs_bool(0, NULL); \
    parse_fcs_long_long_t(NULL, NULL); \
    const_fcs_long_long_t(0, NULL); \
    parse_fcs_p_char_t(NULL, NULL); \
    const_fcs_p_char_t(NULL, NULL); \
  }

/*
 * parser function for string based input
 */
FCSResult fcs_parser(FCS handle, const char *parameters, fcs_bool continue_on_errors)
{
  static const char func_name[] = "fcs_parser";

  FCSResult r = FCS_RESULT_SUCCESS;

  char *cur;
  char *params, *param;
  fcs_int params_strlen;

  params_strlen = strlen(parameters) + 1;
  params = malloc(params_strlen * sizeof(char));
  strncpy(params, parameters, params_strlen);

  cur = params;

  while (cur)
  {
    param = cur;

    cur = strchr(cur, ',');

    if (cur)
    {
      *cur = 0;
      ++cur;
    }

/*    printf("param: %s\n", param);
    printf("cur: %s\n", cur);*/

    IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_a",                   set_box_a,            PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_b",                   set_box_b,            PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("box_c",                   set_box_c,            PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("offset",                  set_offset,           PARSE_SEQ(fcs_float, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("periodicity",             set_periodicity,      PARSE_SEQ(fcs_int, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("near_field_flag",         set_near_field_flag,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("total_particles",         set_total_particles,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("",                        require_virial,       PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("",                        set_tolerance,        PARSE_VAL(fcs_int),                                   PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_energy",        set_tolerance,        CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_ENERGY),        PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_energy_rel",    set_tolerance,        CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_ENERGY_REL),    PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_potential",     set_tolerance,        CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_POTENTIAL),     PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_potential_rel", set_tolerance,        CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_POTENTIAL_REL), PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_field",         set_tolerance,        CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_FIELD),         PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("tolerance_field_rel",     set_tolerance,        CONST_VAL(fcs_int, FCS_TOLERANCE_TYPE_FIELD_REL),     PARSE_VAL(fcs_float));
#ifdef FCS_ENABLE_DIRECT
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_periodic_images",  direct_set_periodic_images,  PARSE_SEQ(fcs_int, 3));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_cutoff",           direct_set_cutoff,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("direct_cutoff_with_near", direct_set_cutoff_with_near, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_EWALD
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_maxkmax", ewald_set_maxkmax,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_kmax", ewald_set_kmax,                PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_r_cut", ewald_set_r_cut,               PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_alpha", ewald_set_alpha,               PARSE_VAL(fcs_float));
#endif
#ifdef FCS_ENABLE_FMM
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_absrel",            fmm_set_absrel,            PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_tolerance_energy",  fmm_set_tolerance_energy,  PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_dipole_correction", fmm_set_dipole_correction, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_potential",         fmm_set_potential,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_cusp_radius",       fmm_set_cusp_radius,       PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_internal_tuning",   fmm_set_internal_tuning,   PARSE_VAL(fcs_long_long_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_maxdepth",          fmm_set_maxdepth,          PARSE_VAL(fcs_long_long_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_unroll_limit",      fmm_set_unroll_limit,      PARSE_VAL(fcs_long_long_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("fmm_balanceload",       fmm_set_balanceload,       PARSE_VAL(fcs_long_long_t));
#endif
#ifdef FCS_ENABLE_MEMD
/*    IF_PARAM_THEN_FUNC3_GOTO_NEXT("", fcs_memd_set_box_size,                  PARSE_VAL(fcs_float), PARSE_VAL(fcs_float), PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_time_step,                 PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_total_number_of_particles, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_local_number_of_particles, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_init_flag,                 PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_mesh_size_1D,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_speed_of_light,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_permittivity,              PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_temperature,               PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("", fcs_memd_set_bjerrum_length,            PARSE_VAL(fcs_float));*/
#endif
#ifdef FCS_ENABLE_MMM1D
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_far_switch_radius", mmm1d_set_far_switch_radius, PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_bessel_cutoff",     mmm1d_set_bessel_cutoff,     PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm1d_maxPWerror",        mmm1d_set_maxPWerror,        PARSE_VAL(fcs_float));
#endif
#ifdef FCS_ENABLE_MMM2D
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_maxPWerror",           mmm2d_set_maxPWerror,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_far_cutoff",           mmm2d_set_far_cutoff,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC2_GOTO_NEXT("mmm2d_dielectric_contrasts", mmm2d_set_dielectric_contrasts, PARSE_VAL(fcs_float), PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_layers_per_node",      mmm2d_set_layers_per_node,      PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("mmm2d_skin",                 mmm2d_set_skin,                 PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("",                           mmm2d_require_total_energy,     PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_P2NFFT
    /* P2NFFT specific parameters */
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_r_cut",                p2nfft_set_r_cut,                     PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_epsI",                 p2nfft_set_epsI,                      PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_epsB",                 p2nfft_set_epsB,                      PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_c",                    p2nfft_set_c,                         PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_alpha",                p2nfft_set_alpha,                     PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_intpol_order",         p2nfft_set_interpolation_order,       PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_near",             p2nfft_set_reg_near,                  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_near_name",        p2nfft_set_reg_near_by_name,          PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_far",              p2nfft_set_reg_far,                   PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_reg_far_name",         p2nfft_set_reg_far_by_name,           PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_p",                    p2nfft_set_p,                         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_require_virial",       p2nfft_require_virial,                PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_ignore_tolerance",     p2nfft_set_ignore_tolerance,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("p2nfft_grid",                 p2nfft_set_grid,                      PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("p2nfft_oversampled_grid",     p2nfft_set_oversampled_grid,          PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p2nfft_cao",                  p2nfft_set_cao,                       PARSE_VAL(fcs_int));

    /* PNFFT specific parameters */
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("pnfft_N",                     p2nfft_set_pnfft_N,                   PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC3_GOTO_NEXT("pnfft_n",                     p2nfft_set_pnfft_n,                   PARSE_VAL(fcs_int), PARSE_VAL(fcs_int), PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_window_name",           p2nfft_set_pnfft_window_by_name,      PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_window",                p2nfft_set_pnfft_window,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_m",                     p2nfft_set_pnfft_m,                   PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_intpol_order",          p2nfft_set_pnfft_interpolation_order, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_phi_hat",           p2nfft_set_pnfft_pre_phi_hat,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_fg_psi",                p2nfft_set_pnfft_fg_psi,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_fft_in_place",          p2nfft_set_pnfft_fft_in_place,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_sort_nodes",            p2nfft_set_pnfft_sort_nodes,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_interlaced",            p2nfft_set_pnfft_interlaced,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_grad_ik",               p2nfft_set_pnfft_grad_ik,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_psi",               p2nfft_set_pnfft_pre_psi,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_fg_psi",            p2nfft_set_pnfft_pre_fg_psi,          PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_full_psi",          p2nfft_set_pnfft_pre_full_psi,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_pre_full_fg_psi",       p2nfft_set_pnfft_pre_full_fg_psi,     PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pnfft_real_f",                p2nfft_set_pnfft_real_f,              PARSE_VAL(fcs_int));

    /* PFFT specific parameters */
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_patience",               p2nfft_set_pfft_patience,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_patience_name",          p2nfft_set_pfft_patience_by_name,      PARSE_VAL(fcs_p_char_t));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_preserve_input",         p2nfft_set_pfft_preserve_input,       PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pfft_tune",                   p2nfft_set_pfft_tune,                 PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_P3M
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_r_cut",           p3m_set_r_cut,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_alpha",           p3m_set_alpha,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_grid",            p3m_set_grid,             PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_cao",             p3m_set_cao,              PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("p3m_require_total_energy",\
                                                         p3m_require_total_energy, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_PEPC
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_epsilon",           pepc_set_epsilon,           PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_theta",             pepc_set_theta,             PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_num_walk_threads",  pepc_set_num_walk_threads,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_dipole_correction", pepc_set_dipole_correction, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_load_balancing",    pepc_set_load_balancing,    PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_npm",               pepc_set_npm,               PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pepc_debug_level",       pepc_set_debug_level,       PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_PP3MG
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_cells_x",        pp3mg_set_cells_x,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_cells_y",        pp3mg_set_cells_y,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_cells_z",        pp3mg_set_cells_z,        PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_ghosts",         pp3mg_set_ghosts,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_degree",         pp3mg_set_degree,         PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_max_particles",  pp3mg_set_max_particles,  PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_max_iterations", pp3mg_set_max_iterations, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_tol",            pp3mg_set_tol,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_distribution",   pp3mg_set_distribution,   PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("pp3mg_discretization", pp3mg_set_discretization, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_VMG
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_max_level",            vmg_set_max_level,            PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_max_iterations",       vmg_set_max_iterations,       PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_smoothing_steps",      vmg_set_smoothing_steps,      PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_cycle_type",                vmg_set_cycle_type,                PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_precision",            vmg_set_precision,            PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_near_field_cells",     vmg_set_near_field_cells,     PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_interpolation_order", vmg_set_interpolation_order, PARSE_VAL(fcs_int));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("vmg_discretization_order", vmg_set_discretization_order, PARSE_VAL(fcs_int));
#endif
#ifdef FCS_ENABLE_WOLF
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("wolf_cutoff", wolf_set_cutoff, PARSE_VAL(fcs_float));
    IF_PARAM_THEN_FUNC1_GOTO_NEXT("wolf_alpha", wolf_set_alpha, PARSE_VAL(fcs_float));
#endif

    DUMMY_REFERENCE_TO_STATIC_FUNCTIONS(param);

    if (r == FCS_RESULT_SUCCESS)
      r = fcsResult_create(FCS_WRONG_ARGUMENT, func_name, "interface (parser): error in parameter string at '%s'!", param); 

    if (FCS_IS_FALSE(continue_on_errors)) break;

next_param:
    ;
  }

  free(params);

  return r;
}


FCSResult fcs_require_virial(FCS handle, fcs_int flag)
{
    char* fnc_name = "fcs_require_virial";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if ((flag<0) || (flag >1))
        return fcsResult_create(FCS_WRONG_ARGUMENT,fnc_name,"flag must be 0 or 1");

    switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
        case FCS_DIRECT:
            return fcs_direct_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_PEPC
        case FCS_PEPC:
            return fcs_pepc_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_FMM
        case FCS_FMM:
            return fcs_fmm_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_P3M
        case FCS_P3M:
            return fcs_p3m_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_MMM1D
        case FCS_MMM1D:
            return fcs_mmm1d_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_MMM2D
        case FCS_MMM2D:
            return fcs_mmm2d_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_MEMD
        case FCS_MEMD:
            return fcs_memd_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_P2NFFT
        case FCS_P2NFFT:
            return fcs_p2nfft_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_PP3MG
        case FCS_PP3MG:
            return fcs_pp3mg_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_VMG
        case FCS_VMG:
            return fcs_vmg_require_virial(handle, flag);
#endif
#ifdef FCS_ENABLE_WOLF
        case FCS_WOLF:
            return fcs_wolf_require_virial(handle, flag);
#endif
        default:
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
    }
}


FCSResult fcs_get_virial(FCS handle, fcs_float *virial)
{
    char* fnc_name = "fcs_get_virial";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");
    if (virial == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied for virial");

    switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
        case FCS_DIRECT:
            return fcs_direct_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_PEPC
        case FCS_PEPC:
            return fcs_pepc_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_FMM
        case FCS_FMM:
            return fcs_fmm_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_P3M
        case FCS_P3M:
            return fcs_p3m_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_MMM1D
        case FCS_MMM1D:
            return fcs_mmm1d_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_MMM2D
        case FCS_MMM2D:
            return fcs_mmm2d_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_MEMD
        case FCS_MEMD:
            return fcs_memd_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_P2NFFT
        case FCS_P2NFFT:
            return fcs_p2nfft_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_PP3MG
        case FCS_PP3MG:
            return fcs_pp3mg_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_VMG
        case FCS_VMG:
            return fcs_vmg_get_virial(handle, virial);
#endif
#ifdef FCS_ENABLE_WOLF
        case FCS_WOLF:
            return fcs_wolf_get_virial(handle, virial);
#endif
        default:
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
    }
}

FCSResult fcs_set_r_cut(FCS handle, fcs_float r_cut)
{
    char* fnc_name = "fcs_set_r_cut";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");

    switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_PEPC
        case FCS_PEPC:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_FMM
        case FCS_FMM:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_P3M
        case FCS_P3M:
            return fcs_p3m_set_r_cut(handle, r_cut);
#endif
#ifdef FCS_ENABLE_MMM1D
        case FCS_MMM1D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_MMM2D
        case FCS_MMM2D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_MEMD
        case FCS_MEMD:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_P2NFFT
        case FCS_P2NFFT:
            return fcs_p2nfft_set_r_cut(handle, r_cut);
#endif
#ifdef FCS_ENABLE_PP3MG
        case FCS_PP3MG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_VMG
        case FCS_VMG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_DIRECT
        case FCS_DIRECT:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
        default:
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
    }
}

FCSResult fcs_unset_r_cut(FCS handle)
{
    char* fnc_name = "fcs_set_r_cut";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");

    switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_PEPC
        case FCS_PEPC:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_FMM
        case FCS_FMM:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_P3M
        case FCS_P3M:
            return fcs_p3m_set_r_cut_tune(handle);
#endif
#ifdef FCS_ENABLE_MMM1D
        case FCS_MMM1D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_MMM2D
        case FCS_MMM2D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_MEMD
        case FCS_MEMD:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_P2NFFT
        case FCS_P2NFFT:
            return fcs_p2nfft_set_r_cut_tune(handle);
#endif
#ifdef FCS_ENABLE_PP3MG
        case FCS_PP3MG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_VMG
        case FCS_VMG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_DIRECT
        case FCS_DIRECT:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
        default:
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
    }
}

FCSResult fcs_get_r_cut(FCS handle, fcs_float *r_cut)
{
    char* fnc_name = "fcs_get_r_cut";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");

    switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_PEPC
        case FCS_PEPC:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_FMM
        case FCS_FMM:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_P3M
        case FCS_P3M:
            return fcs_p3m_get_r_cut(handle, r_cut);
#endif
#ifdef FCS_ENABLE_MMM1D
        case FCS_MMM1D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_MMM2D
        case FCS_MMM2D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_MEMD
        case FCS_MEMD:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_P2NFFT
        case FCS_P2NFFT:
            return fcs_p2nfft_get_r_cut(handle, r_cut);
#endif
#ifdef FCS_ENABLE_PP3MG
        case FCS_PP3MG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_VMG
        case FCS_VMG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
#ifdef FCS_ENABLE_DIRECT
        case FCS_DIRECT:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of r_cut not implemented");
#endif
        default:
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
    }
}

/* setter for tolerance */
FCSResult fcs_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance_value)
{
    char* fnc_name = "fcs_set_tolerance";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");

    switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
        case FCS_DIRECT:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of tolerance not implemented. DIRECT method computes with machine precision.");
#endif
#ifdef FCS_ENABLE_EWALD
        case FCS_EWALD:
            if(tolerance_type == FCS_TOLERANCE_TYPE_FIELD){
              fcs_ewald_set_tolerance_field(handle, tolerance_value);
              return NULL;
            } else
              return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Unsupported tolerance type. EWALD only supports FCS_TOLERANCE_TYPE_FIELD.");
#endif
#ifdef FCS_ENABLE_FMM
        case FCS_FMM:
            if(tolerance_type == FCS_TOLERANCE_TYPE_ENERGY){
              fcs_fmm_set_absrel(handle, FCS_FMM_CUSTOM_ABSOLUTE);
              fcs_fmm_set_tolerance_energy(handle, tolerance_value);
              return NULL;
            } else if(tolerance_type == FCS_TOLERANCE_TYPE_ENERGY_REL){
              fcs_fmm_set_absrel(handle, FCS_FMM_CUSTOM_RELATIVE);
              fcs_fmm_set_tolerance_energy(handle, tolerance_value);
              return NULL;
            }
            else
              return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Unsupported tolerance type. FMM only supports FCS_TOLERANCE_TYPE_ENERGY and FCS_TOLERANCE_TYPE_ENERGY_REL.");
#endif
#ifdef FCS_ENABLE_MEMD
        case FCS_MEMD:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of tolerance not implemented");
#endif
#ifdef FCS_ENABLE_MMM1D
        case FCS_MMM1D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of tolerance not implemented");
#endif
#ifdef FCS_ENABLE_MMM2D
        case FCS_MMM2D:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of tolerance not implemented");
#endif
#ifdef FCS_ENABLE_P2NFFT
        case FCS_P2NFFT:
            return fcs_p2nfft_set_tolerance(handle, tolerance_type, tolerance_value);
#endif
#ifdef FCS_ENABLE_P3M
        case FCS_P3M:
            if(tolerance_type == FCS_TOLERANCE_TYPE_FIELD){
              fcs_p3m_set_tolerance_field(handle, tolerance_value);
              return NULL;
            } else
              return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Unsupported tolerance type. P3M only supports FCS_TOLERANCE_TYPE_FIELD.");
#endif
#ifdef FCS_ENABLE_PEPC
        case FCS_PEPC:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of tolerance not implemented");
#endif
#ifdef FCS_ENABLE_PP3MG
        case FCS_PP3MG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of tolerance not implemented");
#endif
#ifdef FCS_ENABLE_VMG
        case FCS_VMG:
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Change of tolerance not implemented");
#endif
        default:
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
    }
}

FCSResult fcs_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance_value)
{
    char* fnc_name = "fcs_get_tolerance";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");

    switch(fcs_get_method(handle))
    {
#ifdef FCS_ENABLE_DIRECT
        case FCS_DIRECT:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented. DIRECT method computes with machine precision.");
          }
#endif
#ifdef FCS_ENABLE_EWALD
        case FCS_EWALD:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
            return fcs_ewald_get_tolerance_field(handle, tolerance_value);
          }
#endif
#ifdef FCS_ENABLE_FMM
        case FCS_FMM:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented");
          }
#endif
#ifdef FCS_ENABLE_MEMD
        case FCS_MEMD:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented");
          }
#endif
#ifdef FCS_ENABLE_MMM1D
        case FCS_MMM1D:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented");
          }
#endif
#ifdef FCS_ENABLE_MMM2D
        case FCS_MMM2D:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented");
          }
#endif
#ifdef FCS_ENABLE_P2NFFT
        case FCS_P2NFFT:
            return fcs_p2nfft_get_tolerance(handle, tolerance_type, tolerance_value);
#endif
#ifdef FCS_ENABLE_P3M
        case FCS_P3M:
            {
              *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
              fcs_p3m_get_tolerance_field(handle, tolerance_value);
              return NULL;
            }
#endif
#ifdef FCS_ENABLE_PEPC
        case FCS_PEPC:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented");
          }
#endif
#ifdef FCS_ENABLE_PP3MG
        case FCS_PP3MG:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented");
          }
#endif
#ifdef FCS_ENABLE_VMG
        case FCS_VMG:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"Getter for tolerance not implemented");
          }
#endif
        default:
          {
            *tolerance_type = FCS_TOLERANCE_TYPE_UNDEFINED;
            *tolerance_value = -1.0; 
            return fcsResult_create(FCS_WRONG_ARGUMENT, fnc_name, "unknown method chosen");
          }
    }
}

/* getter / setter for value_changed */
fcs_int fcs_get_values_changed(FCS handle)
{
    char* fnc_name = "fcs_get_box_values_changed";

    if (!handle)
    {
        fprintf(stderr, "Error in %s: NULL pointer received, return NULL", fnc_name);
        return -1;
    }
    return handle->values_changed;
}

FCSResult fcs_set_values_changed(FCS handle, fcs_int changed)
{
    char* fnc_name = "fcs_set_method_context";

    if (handle == NULL)
        return fcsResult_create(FCS_NULL_ARGUMENT,fnc_name,"null pointer supplied as handle");

    handle->values_changed = changed;

    return NULL;
}

FCSResult fcs_set_max_particle_move(FCS handle, fcs_float max_particle_move)
{
  char* fnc_name = "fcs_set_max_particle_move";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->set_max_particle_move == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "max. particle move not supported");

  return handle->set_max_particle_move(handle, max_particle_move);
}

FCSResult fcs_set_resort(FCS handle, fcs_int resort)
{
  char* fnc_name = "fcs_set_resort";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->set_resort == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->set_resort(handle, resort);
}


FCSResult fcs_get_resort(FCS handle, fcs_int *resort)
{
  char* fnc_name = "fcs_get_resort";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->get_resort == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->get_resort(handle, resort);
}


FCSResult fcs_get_resort_availability(FCS handle, fcs_int *availability)
{
  char* fnc_name = "fcs_get_resort_availability";

  *availability = 0;

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->get_resort_availability == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->get_resort_availability(handle, availability);
}


FCSResult fcs_get_resort_particles(FCS handle, fcs_int *resort_particles)
{
  char* fnc_name = "fcs_get_resort_particles";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->get_resort_particles == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->get_resort_particles(handle, resort_particles);
}


FCSResult fcs_resort_ints(FCS handle, fcs_int *src, fcs_int *dst, fcs_int n)
{
  char* fnc_name = "fcs_resort_ints";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->resort_ints == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->resort_ints(handle, src, dst, n, fcs_get_communicator(handle));
}


FCSResult fcs_resort_floats(FCS handle, fcs_float *src, fcs_float *dst, fcs_int n)
{
  char* fnc_name = "fcs_resort_floats";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->resort_floats == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->resort_floats(handle, src, dst, n, fcs_get_communicator(handle));
}


FCSResult fcs_resort_bytes(FCS handle, void *src, void *dst, fcs_int n)
{
  char* fnc_name = "fcs_resort_bytes";

  if (handle == NULL)
    return fcsResult_create(FCS_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (handle->resort_bytes == NULL)
    return fcsResult_create(FCS_INCOMPATIBLE_METHOD, fnc_name, "resorting not supported");

  return handle->resort_bytes(handle, src, dst, n, fcs_get_communicator(handle));
}


/*****************************/
/* FORTRAN wrapper functions */
/*****************************/

/* void fcs_init_f (void **handle, const char* method, MPI_Fint communicator, FCSResult return_value) */
FCSResult fcs_init_f (FCS *handle, const char* method, MPI_Fint communicator)
{
  MPI_Comm c_comm = MPI_Comm_f2c(communicator);
  return fcs_init(handle, method, c_comm);
//   FCSResult result = fcs_init((FCS*)handle, method, c_comm);
//   if (NULL == result)
//     return_value = 0;
//   else
//     return_value = fcsResult_getReturnCode(result);
}
