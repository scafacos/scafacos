/*
  Copyright (C) 2011,2012,2013 Olaf Lenz
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

#include "p3mconfig.hpp"
#include "communication.hpp"
#include "fft.hpp"

#include "caf.hpp"

namespace P3M {
  /* DEFAULTS */
  /** Default for number of interpolation points of the charge
      assignment function. */
  const p3m_int P3M_DEFAULT_N_INTERPOL = 32768;
  /** Default for dielectric constant of the surrounding medium:
      metallic. */
  const p3m_float P3M_DEFAULT_EPSILON = 0.0;
  /** Default for offset of first grid point from the origin (lower left
      corner of the simulation box. */
  const p3m_float P3M_DEFAULT_GRIDOFF = 0.5;
  /** Default for the accuracy. */
  const p3m_float P3M_DEFAULT_TOLERANCE_FIELD = 1.0e-3;

  /* DEBUG SWITCHES */
  /* Define to turn on additional sanity checks. */
  /* #define ADDITIONAL_CHECKS */
  /* Define to print out timings at the end of run */
#define P3M_PRINT_TIMINGS
  /*enumeration to specify the type of timings*/
  enum timingEnum{NONE,ESTIMATE_ALL, ESTIMATE_FFT, ESTIMATE_ASSIGNMENT,FULL};

  /* COMPILE TIME SWITCHES */
  /* Differentiation method */
  /** ik-Differentiation */
#define P3M_IK
  /** analytical differentiation */
  //#define P3M_AD
  
  /** Whether to use interlaced version of P3M algorithm. */
  //#define P3M_INTERLACE

  /* Sanity checks */
#if defined(P3M_AD) && defined(P3M_IK)
#error Cannot use P3M_AD and P3M_IK at the same time
#endif

#if defined(P3M_AD) && !defined(P3M_INTERLACE)
#error Cannot use P3M_AD without P3M_INTERLACE
#endif

  /* CONSTANTS */
  /** Search horizon for maximal grid size. */
  const p3m_int P3M_MAX_GRID_DIFF = 10;
  /** This value for epsilon indicates metallic boundary conditions. */
  const p3m_float P3M_EPSILON_METALLIC = 0.0;
  /** precision limit for the r_cut zero */
  const p3m_float P3M_RCUT_PREC = 1.0e-3;
  /** Whether to use the approximation of Abramowitz/Stegun
      AS_erfc_part() for \f$\exp(d^2) erfc(d)\f$, or the C function erfc
      in P3M and Ewald summation. */
#define P3M_USE_ERFC_APPROXIMATION 1

  /** Number of Brillouin zones taken into account in the calculation of
      the optimal influence function (aliasing sums). */
#ifdef P3M_INTERLACE
  const p3m_int P3M_BRILLOUIN = 1;
#else 
  const p3m_int P3M_BRILLOUIN = 0;
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
    p3m_int dim[3];
    /** number of local grid points. */
    p3m_int size;
    /** index of lower left corner of the 
        local grid in the global grid. */
    p3m_int ld_ind[3];
    /** position of the first local grid point. */
    p3m_float ld_pos[3];
    /** dimension of grid inside node domain. */
    p3m_int inner[3];
    /** inner left down grid point */
    p3m_int in_ld[3];
    /** inner up right grid point + (1,1,1)*/
    p3m_int in_ur[3];
    /** number of margin grid points. */
    p3m_int margin[6];
    /** number of margin grid points from neighbour nodes */
    p3m_int r_margin[6];
    /** offset between grid lines of the last dimension */
    p3m_int q_2_off;
    /** offset between grid lines of the two last dimensions */
    p3m_int q_21_off;
  } local_grid_t;

  /** Structure for send/recv grids. */
  typedef struct {
    /** dimension of sub grids to send. */
    p3m_int s_dim[6][3];
    /** left down corners of sub grids to send. */
    p3m_int s_ld[6][3];
    /** up right corners of sub grids to send. */
    p3m_int s_ur[6][3];
    /** sizes for send buffers. */
    p3m_int s_size[6];
    /** dimensionof sub grids to recv. */
    p3m_int r_dim[6][3];
    /** left down corners of sub grids to recv. */
    p3m_int r_ld[6][3];
    /** up right corners of sub grids to recv. */
    p3m_int r_ur[6][3];
    /** sizes for recv buffers. */
    p3m_int r_size[6];
    /** maximal size for send/recv buffers. */
    p3m_int max;
  } send_grid_t;

  /** Structure that holds all data of the P3M algorithm */
  typedef struct {
    /****************************************************
     * SYSTEM PARAMETERS
     ****************************************************/
    /* System size in x,y,z */
    p3m_float box_l[3];
    /* Skin */
    p3m_float skin;

    /****************************************************
     * PARAMETERS OF THE METHOD
     ****************************************************/
    /** Tolerance in the field rms. */
    p3m_float tolerance_field;
    /** number of interpolation points for charge assignment function */
    p3m_int n_interpol;
    /** whether to compute the near field in the method */
    bool near_field_flag;

    /* TUNABLE PARAMETERS */
    /** cutoff radius */
    p3m_float r_cut;
    /** Ewald splitting parameter */
    p3m_float alpha;
    /** number of grid points per coordinate direction (>0). */
    p3m_int grid[3];
    /** charge assignment order ([0,P3M_MAX_CAO]). */
    p3m_int cao;
  
    /* Whether or not it is necessary to retune the method before running it. */
    bool needs_retune;
    /** Whether or not rcut is to be automatically tuned. */
    bool tune_r_cut;
    /** Whether or not alpha is to be automatically tuned. */
    bool tune_alpha;
    /** Whether or not the grid is to be automatically tuned. */
    bool tune_grid;
    /** Whether or not the charge assignment order is to be automatically tuned. */
    bool tune_cao;

    /****************************************************
     * FLAGS TO TURN ON/OFF COMPUTATION OF DIFFERENT COMPONENTS
     ****************************************************/
    /** Whether or not the total energy is to be computed. */
    p3m_int require_total_energy;
    /** The total energy. */
    p3m_float total_energy;
    /** Whether and how timings are to be taken.
     * 0: run without timings
     * 1: partial timing to estimate without results
     * 2: full means all timings and correct results*/
    timingEnum require_timings;
#define TIMING 0
#define TIMING_NEAR 1
#define TIMING_FAR 2
#define TIMING_CA 3
#define TIMING_GATHER 4
#define TIMING_FORWARD 5
#define TIMING_BACK 6
#define TIMING_INFLUENCE 7
#define TIMING_SPREAD 8
#define TIMING_POTENTIALS 9
#define TIMING_FIELDS 10
#define TIMING_DECOMP 11
#define TIMING_COMP 12
#define NUM_TIMINGS 13
    double timings[NUM_TIMINGS];

    /****************************************************
     * DERIVED PARAMETERS
     ****************************************************/
    /** The errors of the method according to the error formula. */
    p3m_float error, ks_error, rs_error;

    /** offset of the first grid point (lower left corner) from the
        coordinate origin ([0,1[). */
    p3m_float grid_off[3];

    /** Cutoff for charge assignment. */
    p3m_float cao_cut[3];
    /** grid constant. */
    p3m_float a[3];
    /** inverse grid constant. */
    p3m_float ai[3];
    /** additional points around the charge assignment grid, for method like dielectric ELC
        creating virtual charges. */
    p3m_float additional_grid[3];

    /****************************************************
     * METHOD DATA
     ****************************************************/
    fft_data_struct fft;
    comm_struct comm;

    /** local grid. */
    local_grid_t local_grid;
    /** real space grid (local) for CA/FFT.*/
    p3m_float *rs_grid;
    /** k space grid (local) for k space calculation and FFT.*/
    p3m_float *ks_grid;
  
    /** number of charged particles */
    p3m_int sum_qpart;
    /** Sum of square of charges */
    p3m_float sum_q2;
    /** square of sum of charges */
    p3m_float square_sum_q;

    /** charge assignment function. */
    CAF *caf;
    CAF::Cache *cafx;
    CAF::Cache *cafy;
    CAF::Cache *cafz;
    /** gradient of charge assignment function */
    CAF *caf_d;
    CAF::Cache *cafx_d;
    CAF::Cache *cafy_d;
    CAF::Cache *cafz_d;

    /** position shift for calc. of first assignment grid point. */
    p3m_float pos_shift;
    /** helper variable for calculation of aliasing sums */
    p3m_int *meshift_x;
    p3m_int *meshift_y;
    p3m_int *meshift_z;

    /** Spatial differential operator in k-space. We use an i*k differentiation. */
    p3m_int *d_op[3];
    /** Force optimized influence function (k-space) */
    p3m_float *g_force;
    /** Energy optimized influence function (k-space) */
    p3m_float *g_energy;

    /** number of permutations in k_space */
    p3m_int ks_pnum;

    /** send/recv grid sizes */
    send_grid_t sm;

    /** Field to store grid points to send. */
    p3m_float *send_grid; 
    /** Field to store grid points to recv */
    p3m_float *recv_grid;
  } data_struct;
}
#endif
