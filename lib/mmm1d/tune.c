/*
  Copyright (C) 2011 Olaf Lenz
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "tune.h"
#include <stdlib.h>
#include <stdio.h> ///@TODO: remove this

#include "types.h"
//#include "../mmm-common/specfunc.h"

#include <math.h>

/***************************************************/
/* FORWARD DECLARATIONS OF INTERNAL FUNCTIONS */
/***************************************************/
static fcs_int mmm1d_determine_bessel_cutoff(mmm1d_data_struct *d, fcs_float switch_radius, fcs_int max_cutoff);
static FCSResult mmm1d_setup_constants(mmm1d_data_struct *d, fcs_float boxz);
static void mmm1d_recalcTables(mmm1d_data_struct *d);
static FCSResult mmm1d_check_system_charges(mmm1d_data_struct *d,
                       fcs_int num_particles, fcs_float *charges);

/***************************************************/
/* IMPLEMENTATION */
/***************************************************/
FCSResult mmm1d_tune(void* rd,
        fcs_int num_particles,
        fcs_float *positions,
        fcs_float *charges) {
  mmm1d_data_struct *d = (mmm1d_data_struct*)rd;
  const char* fnc_name = "mmm1d_tune";
  ///@TODO: check the conditions to make or not the tuning. Probably the better approach: 3 flags in mmm1d_data_struct, one per each settable parameter, to control which one has been modified
  
  /* Require charge neutrality */
  mmm1d_check_system_charges(d, num_particles, charges);
  if (!fcs_float_is_zero(d->total_charge)) {
    return fcs_result_create(FCS_ERROR_LOGICAL_ERROR, fnc_name, "MMM1D requires a zero net charge.");
  }
  
  /* Exit if retuning is unnecessary */
  if (!d->needs_retune) return NULL;
  
  /* check whether the input parameters are sane */
  /* * check whether cutoff too large */
  fcs_float rboxz=d->box_l[2]*d->box_l[2];
  
  if (d->far_switch_radius_2 >= rboxz)
    d->far_switch_radius_2 = 0.8*rboxz;
  
  fprintf(stderr,"tune\n");
  
  mmm1d_setup_constants(d, d->box_l[2]);
  
  fprintf(stderr,"contants\n");
  
  d->bessel_cutoff=mmm1d_determine_bessel_cutoff(rd, sqrt(d->far_switch_radius_2), MMM1D_MAXIMAL_B_CUT);
  
  fprintf(stderr,"bessel\n");
  
  mmm1d_recalcTables(d);
  
  fprintf(stderr,"tables\n");
  
  return NULL;

}

static FCSResult mmm1d_setup_constants(mmm1d_data_struct *d, fcs_float boxz) {
  d->uz  = 1/boxz;
  d->L2  = boxz*boxz;
  d->uz2 = d->uz*d->uz;
  
  d->prefuz2 = d->uz2;
  d->prefL3_i = d->prefuz2*d->uz;
  
  return NULL;
}

static fcs_int mmm1d_determine_bessel_cutoff(mmm1d_data_struct *d, fcs_float switch_radius, fcs_int max_cutoff) {
  /* this calculates an upper bound to all force components and the potential */
  fcs_float err;
  fcs_float rhores = 2*M_PI*d->uz*switch_radius;
  fcs_float pref = 4*d->uz*mmm_dmax(1, 2*M_PI*d->uz);
  fcs_int P = 1;
  do {
    err = pref*mmm_K1(rhores*P)*exp(rhores)/rhores*(P - 1 + 1/rhores);
    P++;
  } while (err > d->maxPWerror && P <= max_cutoff);
  P--;
  return P;
}

static void mmm1d_recalcTables(mmm1d_data_struct *d)
{
  /* polygamma, determine order */
  fcs_int n;
  fcs_float err;
  fcs_float rhomax2nm2, rhomax2 = d->uz2*d->far_switch_radius_2;
  n = 1;
  rhomax2nm2 = 1.0;
  do {
     //printf("recalc 1 %d\n", n);
     //if(!(d->modPsi)){printf("1. aki modPsi en NULL\n");}
     //if((d->modPsi)==NULL){printf("2. aki modPsi en NULL\n");}
    mmm_create_mod_psi_up_to(d->polTaylor, n+1);
    //printf("plofff %d\n", (&(d->modPsi)[0])->n);
    //printf("recalc 2 %d\n", n);
    err = 2*n*fabs(mmm_mod_psi_even(d->polTaylor, n, 0.5))*rhomax2nm2;
    //printf("rank %d, recalc 3 %d, err: %e\n", d->comm.rank, n, err);
    rhomax2nm2 *= rhomax2;
    n++;
  }
  while (err > 0.1*d->maxPWerror);
  //printf("rank %d, rhomax2: %e, recalc_tables: n_modPsi: %d\n", d->comm.rank, rhomax2, (d->polTaylor)->n_modPsi);
}

/* compute the net charge of the system */
FCSResult mmm1d_check_system_charges(mmm1d_data_struct *d,
                fcs_int num_particles, fcs_float *charges) {
  fcs_int i;
  fcs_float local_charge=0., total_charge=0.;
  for (i = 0; i < num_particles; i++) {
    if (!fcs_float_is_zero(charges[i])) {
      local_charge += charges[i];
    }
  }
  fprintf(stderr,"check charges: rank %d\n",d->comm.rank);
  fprintf(stderr, "local_charges %d, %e\n", num_particles, local_charge);
  MPI_Allreduce(&local_charge, &total_charge, 1, FCS_MPI_FLOAT, MPI_SUM, d->comm.mpicomm);
  d->total_charge=total_charge;
  fprintf(stderr, "total charge %e\n", d->total_charge);
  return NULL;
}
