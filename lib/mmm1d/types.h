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
/** How many trial calculations */
#define MMM1D_TEST_INTEGRATIONS 1000

/** Largest reasonable cutoff for Bessel function */
#define MMM1D_MAXIMAL_B_CUT 30

/** Granularity of the radius scan in multiples of box_l[2] */
#define MMM1D_RAD_STEPPING 0.1

/** Default for the switching radius. */
#define MMM1D_DEFAULT_FAR_SWITCH_RADIUS 0.05

/** Default for the cutoff of the bessel sum. */
#define MMM1D_DEFAULT_BESSEL_CUTOFF 5

/** Default for the flag to wether recalculate the Bessel cutoff automatically. */
#define MMM1D_DEFAULT_BESSEL_CALCULATED 1

/** Default for the accuracy. */
#define MMM1D_DEFAULT_REQUIRED_ACCURACY 1e-5

/** Default for the coulomb prefactor, corresponding to vacuum in gaussian units */
#define MMM1D_COULOMB_PREFACTOR 1.0

/** if you define this, the Besselfunctions are calculated up
    to machine precision, otherwise 10^-14, which should be
    definitely enough for daily life. */
#undef MMM1D_BESSEL_MACHINE_PREC

#ifndef MMM1D_BESSEL_MACHINE_PREC
#define MMM1D_K0 LPK0
#define MMM1D_K1 LPK1
#endif

/** Structure that holds all data of the P3M algorithm */
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

  /* TUNABLE PARAMETERS */
  /* cutoff of the bessel sum */
  fcs_int   bessel_cutoff;
  

  /****************************************************
   * FLAGS TO TURN ON/OFF COMPUTATION OF DIFFERENT COMPONENTS
   ****************************************************/
  /* Wether to recalculate the Bessel cutoff automatically.
     If some parameters like the box dimensions change, the current
     Bessel cutoff may not be suitable to achieve the required accuracy
     anymore. If the user did not specify the Bessel cutoff explicitly,
     this flag is 1, and whenever necessary, the Bessel cutoff is
     recalculated */
  fcs_int    bessel_calculated;
  fcs_int    needs_retune;
  fcs_int    radius_reset;
  fcs_int    cutoff_reset;
  fcs_int    error_reset;

  /****************************************************
   * Derived parameters
   ****************************************************/
  fcs_float uz; //  = 1/box_l[2]
  fcs_float L2; //  = box_l[2]*box_l[2]
  fcs_float uz2; // = uz*uz;
  fcs_float prefuz2; // = coulomb.prefactor*uz2
  fcs_float prefL3_i; // = prefuz2*uz
  
  //* table of the Taylor expansions of the modified polygamma functions */
  mmm_data_struct *polTaylor;

  /****************************************************
   * Method data
   ****************************************************/
   /* containers for the local particles */
  fcs_int n_localpart;
  fcs_float *local_charges;
  fcs_float *local_positions;
  
} mmm1d_data_struct;

#endif
