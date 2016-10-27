/*
  Copyright (C) 2011,2012 Olaf Lenz
  
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
/* Implementation of the Ewald summation

   This implementation uses the same notation and variable names as
   Deserno, Holm. "How to mesh up Ewald sums". J. Chem. Phys. 109(18), 1998.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "FCSCommon.h"
#include "ewald/ewald.h"
#include "fcs_ewald.h"


/* method to check whether a given FCS is EWALD */
static FCSResult fcs_ewald_check(FCS handle, const char* fnc_name) {
  if (handle == NULL)
    return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "null pointer supplied as handle");

  if (fcs_get_method(handle) != FCS_METHOD_EWALD)
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "Wrong method chosen. You should choose \"EWALD\".");

  return NULL;
}

FCSResult fcs_ewald_init(FCS handle) {
  const char* fnc_name = "ewald_init";
  FCSResult result;

  result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  handle->shift_positions = 1;

  handle->destroy = fcs_ewald_destroy;
  handle->set_tolerance = fcs_ewald_set_tolerance;
  handle->get_tolerance = fcs_ewald_get_tolerance;
  handle->set_parameter = fcs_ewald_set_parameter;
  handle->print_parameters = fcs_ewald_print_parameters;
  handle->tune = fcs_ewald_tune;
  handle->run = fcs_ewald_run;

  /* initialize ewald struct */
  ewald_data_struct *d;
  if (handle->method_context == NULL) {
    d = malloc(sizeof(ewald_data_struct));
    handle->method_context = d;
  } else {
    d = (ewald_data_struct*)handle->method_context;
  }

  /* MPI communicator */
  d->comm = handle->communicator;
  MPI_Comm_size(d->comm, &d->comm_size);
  MPI_Comm_rank(d->comm, &d->comm_rank);

  /* set up cartesian communicator */
  int node_grid[3] = {0, 0, 0};
  
  /* compute node grid */
  MPI_Dims_create(d->comm_size, 3, node_grid);

  /* create cartesian communicator */
  int periodicity[3] = {1, 1, 1};
  MPI_Cart_create(d->comm, 3, node_grid, periodicity, 
		  1, &d->comm_cart);
  
  d->box_l[0] = 0.0;
  d->box_l[1] = 0.0;
  d->box_l[2] = 0.0;
  d->tolerance_field = EWALD_DEFAULT_TOLERANCE_FIELD;
  d->alpha = 0.0;
  d->r_cut = 0.0;
  d->kmax = 0;
  d->maxkmax = MAXKMAX_DEFAULT;
  /* d->alpha = 1.0; */
  /* d->r_cut = 3.0; */
  /* d->kmax = 40; */

  d->needs_retune = 1;
  d->tune_alpha = 1;
  d->tune_r_cut = 1;
  d->tune_kmax = 1;
  d->G = NULL;

  d->far_fields = NULL;
  d->near_fields = NULL;
  d->far_potentials = NULL;
  d->near_potentials = NULL;

  return NULL;
}

FCSResult fcs_ewald_set_maxkmax(FCS handle, fcs_int maxkmax) {
  const char *fnc_name = "fcs_ewald_set_maxkmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->maxkmax = maxkmax;

  return NULL;
}

FCSResult fcs_ewald_set_maxkmax_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_maxkmax_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->maxkmax = MAXKMAX_DEFAULT;
  
  return NULL;
}

FCSResult fcs_ewald_get_maxkmax(FCS handle, fcs_int *maxkmax) {
  const char *fnc_name = "fcs_ewald_get_kmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *maxkmax = d->maxkmax;

  return NULL;
}

FCSResult fcs_ewald_set_kmax(FCS handle, fcs_int kmax) {
  const char *fnc_name = "fcs_ewald_set_kmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  if (kmax != d->kmax)
    d->needs_retune = 1;
  d->kmax = kmax;
  d->tune_kmax = 0;

  return NULL;
}

FCSResult fcs_ewald_set_kmax_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_kmax_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->needs_retune = 1;
  d->tune_kmax = 1;

  return NULL;
}

FCSResult fcs_ewald_get_kmax(FCS handle, fcs_int *kmax) {
  const char *fnc_name = "fcs_ewald_get_kmax";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *kmax = d->kmax;

  return NULL;
}

FCSResult fcs_ewald_set_r_cut(FCS handle, fcs_float r_cut) {
  const char *fnc_name = "fcs_ewald_set_r_cut";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  if (!fcs_float_is_equal(r_cut, d->r_cut))
    d->needs_retune = 1;
  d->r_cut = r_cut;
  d->tune_r_cut = 0;

  return NULL;
}

FCSResult fcs_ewald_set_r_cut_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_r_cut_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->needs_retune = 1;
  d->tune_r_cut = 1;

  return NULL;
}

FCSResult fcs_ewald_get_r_cut(FCS handle, fcs_float *r_cut) {
  const char *fnc_name = "fcs_ewald_get_r_cut";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *r_cut = d->r_cut;

  return NULL;
}

FCSResult fcs_ewald_set_alpha(FCS handle, fcs_float alpha) {
  const char *fnc_name = "fcs_ewald_set_alpha";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  if (!fcs_float_is_equal(alpha, d->alpha))
    d->needs_retune = 1;
  d->alpha = alpha;
  d->tune_alpha = 0;

  return NULL;
}

FCSResult fcs_ewald_set_alpha_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_alpha_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  d->needs_retune = 1;
  d->tune_alpha = 1;
  return NULL;
}

FCSResult fcs_ewald_get_alpha(FCS handle, fcs_float *alpha) {
  const char *fnc_name = "fcs_ewald_get_alpha";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  *alpha = d->alpha;
  return NULL;
}

FCSResult fcs_ewald_set_tolerance_field(FCS handle, fcs_float tolerance_field) {
  const char *fnc_name = "fcs_ewald_set_tolerance_field";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  if (!fcs_float_is_equal(tolerance_field, d->tolerance_field))
    d->needs_retune = 1;
  d->tolerance_field = tolerance_field;
  return NULL;
}

FCSResult fcs_ewald_set_tolerance_field_tune(FCS handle) {
  const char *fnc_name = "fcs_ewald_set_tolerance_field_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  d->needs_retune = 1;
  d->tolerance_field = EWALD_DEFAULT_TOLERANCE_FIELD;
  return NULL;
}  

FCSResult fcs_ewald_get_tolerance_field(FCS handle, fcs_float* tolerance_field) {
  const char *fnc_name = "fcs_ewald_get_tolerance_field";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
  *tolerance_field = d->tolerance_field;
  return NULL;
}

FCSResult fcs_ewald_set_tolerance(FCS handle, fcs_int tolerance_type, fcs_float tolerance)
{
  const char *fnc_name = "fcs_ewald_set_tolerance";

  if (tolerance_type == FCS_TOLERANCE_TYPE_FIELD)
  {
    fcs_ewald_set_tolerance_field(handle, tolerance);
    return FCS_RESULT_SUCCESS;

  } else return fcs_result_create(FCS_ERROR_NULL_ARGUMENT, fnc_name, "Unsupported tolerance type. EWALD only supports FCS_TOLERANCE_TYPE_FIELD.");
}


FCSResult fcs_ewald_get_tolerance(FCS handle, fcs_int *tolerance_type, fcs_float *tolerance)
{
  *tolerance_type = FCS_TOLERANCE_TYPE_FIELD;
  return fcs_ewald_get_tolerance_field(handle, tolerance);
}


FCSResult fcs_ewald_get_components(FCS handle, 
			 fcs_float** far_fields, 
			 fcs_float** near_fields,
			 fcs_float** far_potentials, 
			 fcs_float** near_potentials) {
  const char* fnc_name =  "fcs_P3M_get_components";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  if (far_fields != NULL) {
    if (d->far_fields == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"far_fields wanted but not computed."); 
    else
      *far_fields = d->far_fields;
  }

  if (near_fields != NULL) {
    if (d->near_fields == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"near_fields wanted but not computed."); 
    else
      *near_fields = d->near_fields;
  }

  if (far_potentials != NULL) {
    if (d->far_potentials == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"far_potentials wanted but not computed."); 
    else
      *far_potentials = d->far_potentials;
  }

  if (near_potentials != NULL) {
    if (d->near_potentials == NULL)
      return fcs_result_create(FCS_ERROR_NULL_ARGUMENT,fnc_name,"near_potentials wanted but not computed."); 
    else
      *near_potentials = d->near_potentials;
  }
  return NULL;
}


FCSResult fcs_ewald_tune(FCS handle,
			 fcs_int num_particles,
			 fcs_float *positions, 
			 fcs_float *charges) {
  const char *fnc_name = "fcs_ewald_tune";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  FCS_INFO(fprintf(stderr, "fcs_ewald_tune() started...\n"));

  /* Handle periodicity */
  const fcs_int *periodicity = fcs_get_periodicity(handle);
  if (! (periodicity[0] && periodicity[1] && periodicity[2]))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
  			    "EWALD requires periodic boundary conditions.");
  
  /* Handle box size */
  const fcs_float *a = fcs_get_box_a(handle);
  const fcs_float *b = fcs_get_box_b(handle);
  const fcs_float *c = fcs_get_box_c(handle);
  if (!fcs_is_orthogonal(a, b, c))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
  			    "EWALD requires the box to be orthorhombic.");
  
  if (!fcs_uses_principal_axes(a, b, c))
    return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name,
  			    "EWALD requires the box vectors to be parallel to the principal axes.");
  
  if (!fcs_float_is_equal(a[0], d->box_l[0])) {
    d->box_l[0] = a[0];
    d->needs_retune = 1;
  }
  if (!fcs_float_is_equal(b[1], d->box_l[1])) {
    d->box_l[1] = b[1];
    d->needs_retune = 1;
  }
  if (!fcs_float_is_equal(c[2], d->box_l[2])) {
    d->box_l[2] = c[2];
    d->needs_retune = 1;
  }
  
  /* Check whether retuning is required */
  if (!d->needs_retune) {
    FCS_INFO(fprintf(stderr, "  Retuning is not required.\n"));
    FCS_INFO(fprintf(stderr, "fcs_ewald_tune() finished.\n"));
    return NULL;
  }
  
  FCS_INFO(fprintf(stderr, "  Retuning is required.\n"));
  
  /* Check whether the input parameters are sane */
  if (!d->tune_r_cut) {
    if (d->r_cut < 0.0) {
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "r_cut is negative!");
    }
    
    if (fcs_float_is_zero(d->r_cut)) {
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "r_cut is too small!");
    }
    
    /* check whether cutoff is larger than half a box length */
    if ((d->r_cut > 0.5*d->box_l[0]) ||
  	(d->r_cut > 0.5*d->box_l[1]) ||
  	(d->r_cut > 0.5*d->box_l[2]))
      return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, "r_cut is larger than half a system box length.");
  }

  /* Init tuning */
  /* count the number of charges */
  fcs_int local_num_charges = 0;
  fcs_int num_charges;
  fcs_float local_q2 = 0;
  fcs_float q2;
  for (fcs_int i = 0; i < num_particles; i++) {
    if (!fcs_float_is_zero(charges[i]))
      local_num_charges++;
    local_q2 += charges[i]*charges[i];
  }
  MPI_Reduce(&local_num_charges, &num_charges, 1, FCS_MPI_INT,
  	     MPI_SUM, 0, d->comm);
  MPI_Reduce(&local_q2, &q2, 1, FCS_MPI_FLOAT,
  	     MPI_SUM, 0, d->comm);
  
  if (d->comm_rank == 0) {
    /* Tune r_cut */
    if (d->tune_r_cut) {
      /* compute the average distance between two charges */
      fcs_float avg_dist = pow((d->box_l[0]*d->box_l[1]*d->box_l[2])
  			       / num_charges, 0.33333);
      /* set r_cut to 3 times the average distance between charges */
      d->r_cut = 3.0 * avg_dist;
      /* tune r_cut to half the box length */
      if (0.5*d->box_l[0] < d->r_cut) d->r_cut = 0.5*d->box_l[0];
      if (0.5*d->box_l[1] < d->r_cut) d->r_cut = 0.5*d->box_l[1];
      if (0.5*d->box_l[2] < d->r_cut) d->r_cut = 0.5*d->box_l[2];
      
      FCS_INFO(fprintf(stderr, "    tuned r_cut to %" FCS_LMOD_FLOAT "f\n", d->r_cut));
    }

    fcs_float estimated_accuracy;
    /* Tune kmax and alpha */
    if (d->tune_kmax) {
      estimated_accuracy = 1.0e10;
      d->kmax = 0;
      while (estimated_accuracy > d->tolerance_field) {
  	d->kmax++;
  	if (d->kmax > d->maxkmax) {
  	  static char msg[255];
  	  sprintf(msg, "EWALD cannot achieve required accuracy with given r_cut=%" FCS_LMOD_FLOAT "f with a kmax<=%" FCS_LMOD_INT "d",
  		  d->r_cut, d->maxkmax);
  	  return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, msg);
  	}

	ewald_tune_alpha(num_charges, q2, d->box_l, d->r_cut, d->kmax, 
			 &d->alpha, &estimated_accuracy, d->tune_alpha);
      }
      FCS_INFO(if (d->tune_alpha) fprintf(stderr, "    tuned alpha to %" FCS_LMOD_FLOAT "f\n", d->alpha));
      FCS_INFO(fprintf(stderr, "    tuned kmax to %" FCS_LMOD_INT "d\n", d->kmax));
      FCS_INFO(fprintf(stderr, "    estimated_accuracy=%" FCS_LMOD_FLOAT "e\n", estimated_accuracy));
    } else {
      /* Tune only alpha */
      ewald_tune_alpha(num_charges, q2, d->box_l, d->r_cut, d->kmax, 
		       &d->alpha, &estimated_accuracy, d->tune_alpha);
      FCS_INFO(fprintf(stderr, "    tuned alpha to %" FCS_LMOD_FLOAT "f\n", d->alpha));
      FCS_INFO(fprintf(stderr, "    estimated_accuracy=%" FCS_LMOD_FLOAT "e\n", estimated_accuracy));
      if (estimated_accuracy > d->tolerance_field) {
  	char msg[255];
	if (d->tune_alpha)
	  sprintf(msg, "EWALD cannot achieve required accuracy with given kmax=%" FCS_LMOD_INT "d and r_cut=%" FCS_LMOD_FLOAT "f", d->kmax, d->r_cut);
	else
	  sprintf(msg, "EWALD cannot achieve required accuracy with given kmax=%" FCS_LMOD_INT "d, r_cut=%" FCS_LMOD_FLOAT "f and alpha=%" FCS_LMOD_FLOAT "f", d->kmax, d->r_cut, d->alpha);
  	return fcs_result_create(FCS_ERROR_WRONG_ARGUMENT, fnc_name, msg);
      }
    }
  }

  /* broadcast the tuned parameters */
  if (d->tune_r_cut)
    MPI_Bcast(&d->r_cut, 1, FCS_MPI_FLOAT, 0, d->comm);
  if (d->tune_alpha)
    MPI_Bcast(&d->alpha, 1, FCS_MPI_FLOAT, 0, d->comm);
  if (d->tune_kmax)
    MPI_Bcast(&d->kmax, 1, FCS_MPI_INT, 0, d->comm);

  /* Compute influence function */
  /* The influence function is 1/V*g_tilde(k_vec)*gamma_tilde(k_vec) */
  fcs_int size = d->kmax+1;
  size = size*size*size;
  if (d->G == NULL)
    d->G = (fcs_float*)malloc(sizeof(fcs_float)*size);
  else
    d->G = (fcs_float*)realloc(d->G, sizeof(fcs_float)*size);

  for (fcs_int nx=0; nx <= d->kmax; nx++)
    for (fcs_int ny=0; ny <= d->kmax; ny++)
      for (fcs_int nz=0; nz <= d->kmax; nz++) {
	// system length vector L
	const fcs_float Lx = d->box_l[0];
	const fcs_float Ly = d->box_l[1];
	const fcs_float Lz = d->box_l[2];
	// volume V
	const fcs_float V = Lx*Ly*Lz;
	// reciprocal vector k
	const fcs_float kx = 2.0*M_PI*nx / Lx;
	const fcs_float ky = 2.0*M_PI*ny / Ly;
	const fcs_float kz = 2.0*M_PI*nz / Lz;
	// k^2
	const fcs_float ksqr = kx*kx + ky*ky + kz*kz;
	// alpha
	const fcs_float alpha = d->alpha;

  	if ((nx==0 && ny==0 && nz==0) || 
	    nx*nx + ny*ny + nz*nz > d->kmax*d->kmax)
	  // set influence function to 0 for k=0 and outside the sphere
  	  d->G[linindex(nx, ny, nz, d->kmax)] = 0.0;
  	else 
	  // compute the influence function inside the sphere
	  /* Compare to Deserno, Holm (1998) eqs. (5) and (15) */
	  d->G[linindex(nx, ny, nz, d->kmax)]
	    = 1.0/V * 4.0*M_PI/ksqr * exp(-ksqr/(4*SQR(alpha)));
      }

  d->needs_retune = 0;

  FCS_INFO(fprintf(stderr, "fcs_ewald_tune() finished.\n"));
  return NULL;
}


FCSResult fcs_ewald_run(FCS handle,
			fcs_int num_particles,
			fcs_float *positions, 
			fcs_float *charges,
			fcs_float *fields,
			fcs_float *potentials) {
  const char *fnc_name = "fcs_ewald_run";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  /* First run tune */
  fcs_ewald_tune(handle, num_particles, positions, charges);

  FCS_INFO(fprintf(stderr, "fcs_ewald_run() started...\n"));
  FCS_INFO(fprintf(stderr,						\
		   "    system parameters: box_l=(%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f)\n",	\
		   d->box_l[0], d->box_l[1], d->box_l[2]));
  FCS_INFO(fprintf(stderr,						\
		   "    ewald params: r_cut=%" FCS_LMOD_FLOAT "f, alpha=%" FCS_LMOD_FLOAT "f kmax=%" FCS_LMOD_INT "d\n",		\
		   d->r_cut, d->alpha, d->kmax));

  if (fields != NULL) {
    if (d->far_fields == NULL)
      d->far_fields = malloc(num_particles*sizeof(fcs_float)*3);
    else
      d->far_fields = realloc(d->far_fields, 
			      num_particles*sizeof(fcs_float)*3);
  }
  
  if (potentials != NULL) {
    if (d->far_potentials == NULL)
      d->far_potentials = malloc(num_particles*sizeof(fcs_float));
      else
  	d->far_potentials = realloc(d->far_potentials,
				    num_particles*sizeof(fcs_float));
  }

  /* Compute far field component */
  ewald_compute_kspace(d, num_particles, positions, charges,
  		       fields==NULL ? NULL: d->far_fields,
		       potentials==NULL ? NULL : d->far_potentials
		       );


  if (fields != NULL) {
    if (d->near_fields == NULL)
      d->near_fields = malloc(num_particles*sizeof(fcs_float)*3);
    else
      d->near_fields = realloc(d->near_fields,
			       num_particles*sizeof(fcs_float)*3);
  }
  if (potentials != NULL) {
    if (d->near_potentials == NULL)
      d->near_potentials = malloc(num_particles*sizeof(fcs_float));
    else
      d->near_potentials = realloc(d->near_potentials,
				   num_particles*sizeof(fcs_float));
  }

  fcs_int max_local_particles = fcs_get_max_local_particles(handle);
  if (num_particles > max_local_particles) max_local_particles = num_particles;

  /* Compute near field component */
  ewald_compute_rspace(d, num_particles, max_local_particles, positions, charges,
		       fields==NULL ? NULL : d->near_fields, 
		       potentials==NULL ? NULL : d->near_potentials);

  /* Add up components */
  if (fields != NULL) 
    for (fcs_int pid=0; pid < num_particles; pid++) {
      fields[3*pid] = d->far_fields[3*pid] + d->near_fields[3*pid];
      fields[3*pid+1] = d->far_fields[3*pid+1] + d->near_fields[3*pid+1];
      fields[3*pid+2] = d->far_fields[3*pid+2] + d->near_fields[3*pid+2];

      /* printf("%d: far_field=(%e, %e, %e)\n", pid,  */
      /* 	     d->far_fields[3*pid], */
      /* 	     d->far_fields[3*pid+1], */
      /* 	     d->far_fields[3*pid+2]); */
      /* printf("%d: near_field=(%e, %e, %e)\n", pid,  */
      /* 	     d->near_fields[3*pid], */
      /* 	     d->near_fields[3*pid+1], */
      /* 	     d->near_fields[3*pid+2]); */
    }

  if (potentials != NULL)
    for (fcs_int pid=0; pid < num_particles; pid++) {
      potentials[pid] = d->far_potentials[pid] + d->near_potentials[pid];
      /* printf("%d: potential=%e (near=%e, far=%e)\n", pid,  */
      /* 	     potentials[pid], d->near_potentials[pid], d->far_potentials[pid]); */
    }

  FCS_INFO(fprintf(stderr, "fcs_ewald_run() finished.\n"));
  return NULL;
}


FCSResult fcs_ewald_set_parameter(FCS handle, fcs_bool continue_on_errors, char **current, char **next, fcs_int *matched)
{
//   const char *fnc_name = "fcs_ewald_set_parameter";

  char *param = *current;
  char *cur = *next;

  *matched = 0;

  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_maxkmax", ewald_set_maxkmax, FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_kmax",    ewald_set_kmax,    FCS_PARSE_VAL(fcs_int));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_r_cut",   ewald_set_r_cut,   FCS_PARSE_VAL(fcs_float));
  FCS_PARSE_IF_PARAM_THEN_FUNC1_GOTO_NEXT("ewald_alpha",   ewald_set_alpha,   FCS_PARSE_VAL(fcs_float));

  return FCS_RESULT_SUCCESS;

next_param:
  *current = param;
  *next = cur;

  *matched = 1;

  return FCS_RESULT_SUCCESS;
}


FCSResult
fcs_ewald_print_parameters(FCS handle) {
  const char *fnc_name = "fcs_ewald_print_parameters";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;
  
  ewald_data_struct *d = (ewald_data_struct*)handle->method_context;

  printf("ewald tolerance_field=%" FCS_LMOD_FLOAT "f\n", d->tolerance_field);

  if (d->needs_retune && d->tune_kmax)
    printf("ewald kmax=(to be tuned)\n");
  else if (d->tune_kmax)
    printf("ewald kmax=%" FCS_LMOD_INT "d (tuned)\n", d->kmax);
  else
    printf("ewald kmax=%" FCS_LMOD_INT "d\n", d->kmax);

  printf("ewald maxkmax=%" FCS_LMOD_INT "d\n", d->maxkmax);

  if (d->needs_retune && d->tune_r_cut)
    printf("ewald r_cut=(to be tuned)\n");
  else if (d->tune_r_cut)
    printf("ewald r_cut=%" FCS_LMOD_FLOAT "f (tuned)\n", d->r_cut);
  else
    printf("ewald r_cut=%" FCS_LMOD_FLOAT "f\n", d->r_cut);

  if (d->needs_retune && d->tune_alpha)
    printf("ewald alpha=(to be tuned)\n");
  else if (d->tune_alpha)
    printf("ewald alpha=%" FCS_LMOD_FLOAT "f (tuned)\n", d->alpha);
  else
    printf("ewald alpha=%" FCS_LMOD_FLOAT "f\n", d->alpha);
  return NULL;
}

/* safe free */
static void sfree(void* ptr) {
  if (ptr != NULL)
    free(ptr);
}

FCSResult fcs_ewald_destroy(FCS handle) {
  const char *fnc_name = "fcs_ewald_destroy";
  FCSResult result = fcs_ewald_check(handle, fnc_name);
  if (result != NULL) return result;

  if (handle->method_context != NULL) {
    ewald_data_struct *d = (ewald_data_struct*)handle->method_context;
    sfree(d->G);
    sfree(d->far_fields);
    sfree(d->near_fields);
    sfree(d->far_potentials);
    sfree(d->near_potentials);
    sfree(d);
  }
  return NULL;
}
