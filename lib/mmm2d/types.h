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

#ifndef _MMM2D_TYPES_H
#define _MMM2D_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "communication.h"
#include "common/mmm-common/specfunc.h"

/* DEFAULTS */
/** Default for the accuracy. */
#define MMM2D_DEFAULT_REQUIRED_ACCURACY 1e-5
/** Default for the far cutoff. */
#define MMM2D_DEFAULT_FAR_CUTOFF 0.0
/** Default for the far cutoff calculated flag */
#define MMM2D_DEFAULT_FAR_CUTOFF_CALCULATED 1
/** Default for the mid-top dielectric contrast. */
#define MMM2D_DEFAULT_DCONTRAST_TOP 0.0
/** Default for the mid-bottom dielectric contrast. */
#define MMM2D_DEFAULT_DCONTRAST_BOT 0.0
/** Default for the skin. */
#define MMM2D_DEFAULT_SKIN 0.0
/** Default for n_layers. */
#define MMM2D_DEFAULT_N_LAYERS 1

/** Largest reasonable cutoff for far formula. A double cannot overflow
    with this value. */
#define MMM2D_MAXIMAL_FAR_CUT 100

/** Largest reasonable cutoff for Bessel function. The Bessel functions
    are quite slow, so do not make too large. */
#define MMM2D_MAXIMAL_B_CUT 50

/** Largest reasonable order of polygamma series. These are pretty fast,
    so use more of them. Also, the real cutoff is determined at run time,
    so normally we are faster */
#define MMM2D_MAXIMAL_POLYGAMMA 100

/** internal relative precision of far formula. This controls how many
    p,q vectors are done at once. This has nothing to do with the effective
    precision, but rather controls how different values can be we add up without
    loosing the smallest values. In principle one could choose smaller values, but
    that would not make things faster */
#define MMM2D_FARRELPREC 1e-6

/** number of steps in the complex cutoff table */
#define MMM2D_COMPLEX_STEP 16
/** map numbers from 0 to 1/2 onto the complex cutoff table
    (with security margin) */
#define MMM2D_COMPLEX_FAC (MMM2D_COMPLEX_STEP/(.5 + 0.01))

/** Product decomposition data organization
    For the cell blocks
    it is assumed that the lower blocks part is in the lower half.
    This has to have positive sign, so that has to be first. */
#define MMM2D_POQESP 0
#define MMM2D_POQECP 1
#define MMM2D_POQESM 2
#define MMM2D_POQECM 3

#define MMM2D_PQESSP 0
#define MMM2D_PQESCP 1
#define MMM2D_PQECSP 2
#define MMM2D_PQECCP 3
#define MMM2D_PQESSM 4
#define MMM2D_PQESCM 5
#define MMM2D_PQECSM 6
#define MMM2D_PQECCM 7

#define MMM2D_QQEQQP 0
#define MMM2D_QQEQQM 1

#define MMM2D_ABEQQP 0
#define MMM2D_ABEQZP 1
#define MMM2D_ABEQQM 2
#define MMM2D_ABEQZM 3

typedef struct {
  fcs_float s, c;
} mmm2d_SCCache;

/** Structure that holds all data of the MMM2D algorithm */
typedef struct {
  /****************************************************
   * SYSTEM PARAMETERS
   ****************************************************/
  /* System size in x,y,z and its inverse */
  fcs_float box_l[3], box_l_i[3];
  
  fcs_float total_charge;
  
  /* lower z coordinate of the local box */
  fcs_float my_bottom;
  /****************************************************
   * PARAMETERS OF THE METHOD
   ****************************************************/
  /** maximal error of a pairwise interaction. Used at least by the
      near formula, since this does the error control at run time */
  fcs_float maxPWerror;
  /** far formula cutoff and its square */
  fcs_float far_cut, far_cut2;
  /** if nonzero, far_cut has been calculated by MMM2D_tune_far, and will be
      recalculated automatically, if important parameters, such as the number of cells, change. If this
      is zero, the far cutoff has been set explicitly by the user and will not be touched by Espresso. */
  fcs_int calculate_far;
  
  /** if zero, n_layers has been set explicitly. If nonzero, n_layers is set accordingly to the number of nodes available */
  fcs_int calculate_n_layers;
  
  /** Whether or not the total energy is to be computed. */
  fcs_int require_total_energy;
  /** The total energy. */
  fcs_float total_energy;
  
  /// whether there is dielectric contrast
  fcs_int dielectric_contrast_on;
  /** dielectric contrasts at the bottom and top of the simulation cell */
  fcs_float delta_mid_top, delta_mid_bot, delta_mult;
  
  
  /****************************************
   * AUXILIARY VARIABLES
   ****************************************/

  /** up to that error the sums in the NF are evaluated */
  fcs_float part_error;

  /** cutoffs for the bessel sum */
  IntList besselCutoff; // = {NULL, 0, 0};

  /** cutoffs for the complex sum */
  fcs_int  complexCutoff[MMM2D_COMPLEX_STEP + 1];
  /** bernoulli numbers divided by n */
  SizedList bon; // = {NULL, 0, 0};

  /** inverse box dimensions */
  fcs_float ux, ux2, uy, uy2, uz;
  
  /** sin/cos caching */ 
  mmm2d_SCCache *scxcache;
  fcs_int    n_scxcache;
  mmm2d_SCCache *scycache;
  fcs_int    n_scycache;

  /** height of the layers and minimal z for far formula */
  fcs_float layer_h;
  
  fcs_int layers_per_node, n_total_layers;
  //fcs_int n_layers;
  
  /* Skin related */
  fcs_float skin, min_far, max_near;
  
  fcs_int needs_tuning;
  
  mmm2d_comm_struct comm;
  
  //* table of the Taylor expansions of the modified polygamma functions */
  mmm_data_struct *polTaylor;

  fcs_float self_energy;
  
  /* containers for the local particles */
  fcs_int n_localpart;
  fcs_float *local_charges;
  fcs_float *local_positions;
  fcs_int *zslices_nparticles;
  
  /* temporary buffers for product decomposition */
  fcs_float *partblk;
  /* for all local cells including ghosts */
  fcs_float *lclcblk;
  /* collected data from the cells above the top neighbor
    of a cell rsp. below the bottom neighbor
    (P=below, M=above, as the signs in the exp). */
  fcs_float *gblcblk;

  /* contribution from the image charges */
  fcs_float lclimge[8];
  
} mmm2d_data_struct;

#endif
