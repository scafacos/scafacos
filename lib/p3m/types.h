/*
  Copyright (C) 2011,2012 Olaf Lenz
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
#ifndef _P3M_TYPES_H
#define _P3M_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "communication.h"
#include "fft.h"

/* DEFAULTS */
/** Default for number of interpolation points of the charge
    assignment function. */
#define P3M_DEFAULT_N_INTERPOL 32768
/** Default for dielectric constant of the surrounding medium:
    metallic. */
#define P3M_DEFAULT_EPSILON 0.0
/** Default for offset of first grid point from the origin (lower left
    corner of the simulation box. */
#define P3M_DEFAULT_GRIDOFF 0.5
/** Default for the accuracy. */
#define P3M_DEFAULT_TOLERANCE_FIELD 1.0e-3

/* COMPILE TIME SWITCHES */
/* Differentiation method */
/** ik-Differentiation */
//#define P3M_IK
/** analytical differentiation */
#define P3M_AD

/** Whether to use interlaced version of P3M alogorithm. */
#define P3M_INTERLACE

/* Sanity checks */
#if defined(P3M_AD) && defined(P3M_IK)
#error Cannot use P3M_AD and P3M_IK at the same time
#endif

#if defined(P3M_AD) && !defined(P3M_INTERLACE)
#error Cannot use P3M_AD without P3M_INTERLACE
#endif

/* CONSTANTS */
/** maximal charge assignment order available */
#define P3M_MAX_CAO 7
/** This value for epsilon indicates metallic boundary conditions. */
#define P3M_EPSILON_METALLIC 0.0
/** precision limit for the r_cut zero */
#define P3M_RCUT_PREC 1e-3
/** granularity of the time measurement */
#define P3M_TIME_GRAN 2
/** Whether to use the approximation of Abramowitz/Stegun
    AS_erfc_part() for \f$\exp(d^2) erfc(d)\f$, or the C function erfc
    in P3M and Ewald summation. */
#define P3M_USE_ERFC_APPROXIMATION 1

/** Number of Brillouin zones taken into account in the calculation of
    the optimal influence function (aliasing sums). */
#ifdef P3M_INTERLACE
static const fcs_int P3M_BRILLOUIN = 1;
#else 
static const fcs_int P3M_BRILLOUIN = 0;
#endif

/* Index helpers for direct and reciprocal space
 * After the FFT the data is in order YZX, which
 * means that Y is the slowest changing index.
 * The defines are here to not get confused and
 * be able to easily change the order.
 */
#define RX 0
#define RY 1
#define RZ 2
#define KY 0
#define KZ 1
#define KX 2 

#define REQ_P3M_GATHER 100
#define REQ_P3M_SPREAD 101

/***************************************************/
/* DATA TYPES */
/***************************************************/
/** Structure for local grid parameters. */
typedef struct {
  /* local grid characterization. */
  /** dimension (size) of local grid. */
  fcs_int dim[3];
  /** number of local grid points. */
  fcs_int size;
  /** index of lower left corner of the 
      local grid in the global grid. */
  fcs_int ld_ind[3];
  /** position of the first local grid point. */
  fcs_float ld_pos[3];
  /** dimension of grid inside node domain. */
  fcs_int inner[3];
  /** inner left down grid point */
  fcs_int in_ld[3];
  /** inner up right grid point + (1,1,1)*/
  fcs_int in_ur[3];
  /** number of margin grid points. */
  fcs_int margin[6];
  /** number of margin grid points from neighbour nodes */
  fcs_int r_margin[6];
  /** offset between grid lines of the last dimension */
  fcs_int q_2_off;
  /** offset between grid lines of the two last dimensions */
  fcs_int q_21_off;
} ifcs_p3m_local_grid;

/** Structure for send/recv grids. */
typedef struct {
  /** dimension of sub grids to send. */
  fcs_int s_dim[6][3];
  /** left down corners of sub grids to send. */
  fcs_int s_ld[6][3];
  /** up right corners of sub grids to send. */
  fcs_int s_ur[6][3];
  /** sizes for send buffers. */
  fcs_int s_size[6];
  /** dimensionof sub grids to recv. */
  fcs_int r_dim[6][3];
  /** left down corners of sub grids to recv. */
  fcs_int r_ld[6][3];
  /** up right corners of sub grids to recv. */
  fcs_int r_ur[6][3];
  /** sizes for recv buffers. */
  fcs_int r_size[6];
  /** maximal size for send/recv buffers. */
  fcs_int max;
} ifcs_p3m_send_grid;

/** Structure that holds all data of the P3M algorithm */
typedef struct {
  /****************************************************
   * SYSTEM PARAMETERS
   ****************************************************/
  /* System size in x,y,z */
  fcs_float box_l[3];
  /* Skin */
  fcs_float skin;

  /****************************************************
   * PARAMETERS OF THE METHOD
   ****************************************************/
  /** Tolerance in the field rms. */
  fcs_float tolerance_field;
  /** number of interpolation points for charge assignment function */
  fcs_int n_interpol;
  /** whether to compute the near field in the method */
  fcs_int near_field_flag;

  /* TUNABLE PARAMETERS */
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  fcs_float r_cut;
  /** unscaled \ref alpha_L for use with fast inline functions only */
  fcs_float alpha;
  /** number of grid points per coordinate direction (>0). */
  fcs_int    grid[3];
  /** charge assignment order ([0,P3M_MAX_CAO]). */
  fcs_int    cao;
  
  /* Whether or not it is necessary to retune the method before running it. */
  int needs_retune;
  /** Whether or not rcut is to be automatically tuned. */
  int tune_r_cut;
  /** Whether or not alpha is to be automatically tuned. */
  int tune_alpha;
  /** Whether or not the grid is to be automatically tuned. */
  int tune_grid;
  /** Whether or not the charge assignment order is to be automatically tuned. */
  int tune_cao;

  /****************************************************
   * FLAGS TO TURN ON/OFF COMPUTATION OF DIFFERENT COMPONENTS
   ****************************************************/
  /** Whether or not the total energy is to be computed. */
  fcs_int require_total_energy;
  /** The total energy. */
  fcs_float total_energy;

  /****************************************************
   * DERIVED PARAMETERS
   ****************************************************/
  /** The errors of the method according to the error formula. */
  fcs_float error, ks_error, rs_error;

  /** offset of the first grid point (lower left corner) from the
      coordinate origin ([0,1[). */
  fcs_float grid_off[3];

  /** Cutoff for charge assignment. */
  fcs_float cao_cut[3];
  /** grid constant. */
  fcs_float a[3];
  /** inverse grid constant. */
  fcs_float ai[3];
  /** additional points around the charge assignment grid, for method like dielectric ELC
      creating virtual charges. */
  fcs_float additional_grid[3];

  /****************************************************
   * METHOD DATA
   ****************************************************/
  ifcs_fft_data_struct fft;
  ifcs_p3m_comm_struct comm;

  /** local grid. */
  ifcs_p3m_local_grid local_grid;
  /** real space grid (local) for CA/FFT.*/
  fcs_float *rs_grid;
  /** k space grid (local) for k space calculation and FFT.*/
  fcs_float *ks_grid;
  
  /** number of charged particles */
  fcs_int sum_qpart;
  /** Sum of square of charges */
  fcs_float sum_q2;
  /** square of sum of charges */
  fcs_float square_sum_q;

  /** interpolation of the charge assignment function. */
  fcs_float *int_caf;
  /** interpolation of the gradient of charge assignment function */
  fcs_float *int_caf_d;

  /** position shift for calc. of first assignment grid point. */
  fcs_float pos_shift;
  /** helper variable for calculation of aliasing sums */
  fcs_int *meshift_x;
  fcs_int *meshift_y;
  fcs_int *meshift_z;

  /** Spatial differential operator in k-space. We use an i*k differentiation. */
  fcs_int *d_op[3];
  /** Force optimized influence function (k-space) */
  fcs_float *g_force;
  /** Energy optimized influence function (k-space) */
  fcs_float *g_energy;

  /** number of permutations in k_space */
  fcs_int ks_pnum;

  /** send/recv grid sizes */
  ifcs_p3m_send_grid sm;

  /** Field to store grid points to send. */
  fcs_float *send_grid; 
  /** Field to store grid points to recv */
  fcs_float *recv_grid;
} ifcs_p3m_data_struct;

#endif
