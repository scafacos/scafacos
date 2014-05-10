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
#ifndef _P3M_TYPES_HPP
#define _P3M_TYPES_HPP

#include "p3mconfig.hpp"

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
//#define P3M_PRINT_TIMINGS
/*enumeration to specify the type of timings*/
enum timingEnum {
	NONE, ESTIMATE_ALL, ESTIMATE_FFT, ESTIMATE_ASSIGNMENT, FULL
};

/* COMPILE TIME SWITCHES */
/* Differentiation method */
/** ik-Differentiation */
//#define P3M_IK
/** analytical differentiation */
#define P3M_AD
/** Whether to use interlaced version of P3M algorithm. */
#define P3M_INTERLACE
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
    /** structure for near computation parameters: alpha and the potential shift. */
    typedef struct {p3m_float alpha; p3m_float potentialOffset;} near_params_t;
    
/** Structure for local grid parameters. */
struct local_grid_t {
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
};

/** Structure for send/recv grids. */
struct send_grid_t {
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
};

struct Parameters {
	/** cutoff radius */
	p3m_float r_cut;
	/** Ewald splitting parameter */
	p3m_float alpha;
	/** number of grid points per coordinate direction (>0). */
	p3m_int grid[3];
	/** charge assignment order ([0,P3M_MAX_CAO]). */
	p3m_int cao;
};

}
#endif
