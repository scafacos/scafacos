/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
#ifndef _MMM1D_TYPES_H
#define _MMM1D_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "communication.h"
#include "common/mmm-common/specfunc.h"

/* DEFAULTS */
/** Largest cutoff for Bessel function. Shouldn't be much larger than 30,
    otherwise we get numerical instabilities. */
#define MMM1D_DEFAULT_MAXIMAL_B_CUT 30

/** Default for the switching radius, in multiples of box-z. Should be optimal in most cases. */
#define MMM1D_DEFAULT_FAR_SWITCH_RADIUS 0.33

/** Default for the accuracy. */
#define MMM1D_DEFAULT_REQUIRED_ACCURACY 1e-5

/** if you define this, the Besselfunctions are calculated up
    to machine precision, otherwise 10^-14, which should be
    definitely enough for daily life. */
#undef MMM1D_BESSEL_MACHINE_PREC

#ifndef MMM1D_BESSEL_MACHINE_PREC
#define MMM1D_K0 LPK0
#define MMM1D_K1 LPK1
#endif

/** Structure that holds all data of the mmm1d algorithm */
typedef struct {
  
  mmm1d_comm_struct comm;
  
  /****************************************************
   * SYSTEM PARAMETERS
   ****************************************************/
  /* System size in x,y,z */
  fcs_float box_l[3];
  fcs_float total_charge;

  /****************************************************
   * PARAMETERS OF THE METHOD
     Most of the parameters can also be tuned automatically. Unlike
     P3M, this tuning is redone automatically whenever parameters change,
     but not immediately if you set this parameters.
   ****************************************************/
  /* square of the switching radius */
  fcs_float far_switch_radius_2;
  /* required accuracy */
  fcs_float maxPWerror;
  /* maximal possible cutoff of the Bessel sum */
  fcs_int   bessel_cutoff;

  /****************************************************
   * Derived parameters
   ****************************************************/
  fcs_float uz;   //  = 1/box_l[2]
  fcs_float L2;   //  = box_l[2]*box_l[2]
  fcs_float uz2;  // = uz*uz;
  fcs_float L3_i; // = uz2*uz
  
  //* table of the Taylor expansions of the modified polygamma functions */
  mmm_data_struct *polTaylor;

  double *bessel_radii;

  /****************************************************
   * Method data
   ****************************************************/
   /* containers for the local particles */
  fcs_int n_localpart;
  fcs_float *local_charges;
  fcs_float *local_positions;
  
} mmm1d_data_struct;

#endif
