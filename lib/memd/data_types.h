/*
 Copyright (C) 2010/2011/2012 Florian Fahrenberger
 
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

#ifndef _MEMD_DATA_TYPES_H
#define _MEMD_DATA_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>
#include "FCSCommon.h"
#include "common/gridsort/gridsort.h"

/* number of dimensions */
#define SPACE_DIM 3
/* number of directions */
#define NDIRS 6



/** \name Data structures */
/*@{*/

typedef fcs_int     t_ivector[SPACE_DIM]; /* integer vector for position in grid */
typedef fcs_float   t_dvector[SPACE_DIM]; /* double vector for fields etc. */
typedef fcs_int     t_dirs[NDIRS];     /* integer vector for directions */


typedef struct {
    /** particle charge */
    fcs_float  q;
    /** particle position */
    fcs_float*  r;
    /** particle force */
    fcs_float*  f;
    /** particle velocity */
    fcs_float*  v;
    /** particle ID (not yet used) */
    fcs_int*    identity;
} memd_particle;

typedef struct {
    /** particle data */
    memd_particle* part;
    /** number of particles in cell */
    fcs_int n;
} memd_cell;

typedef struct {
    /** cell data array */
    memd_cell** cell;
    /** number of cells in list */
    fcs_int     n;
} memd_cell_list;


/** Structure of local lattice parameters. */
typedef struct {
    t_dvector left_down_position;        /* spatial position of left down grid point */  
    t_dvector upper_right_position;      /* spatial positon of upper right grid point */ 
    fcs_int inner_left_down[SPACE_DIM];              /* inner left down grid point    */
    fcs_int inner_up_right[SPACE_DIM];               /* inner up right grid point + (1,1,1) */
    fcs_int halo_left_down[SPACE_DIM];               /* halo-region left down grid point  */
    fcs_int halo_upper_right[SPACE_DIM];             /* halo-region up right global grid point  */
    fcs_int margin[SPACE_DIM*2];             /* number of margin mesh points (even index - left, odd - right). */
    t_ivector dim;                       /* grid dimension (size + glue_patch region) of local mesh.  */
    t_ivector size;                      /* dimension of mesh inside node domain. */
    fcs_int volume;                          /* number of lattice sites in local domain */
} lattice_parameters;

/** surface_patch structure for communication. */
typedef struct {
    fcs_int offset;        /* source offset for the site index */
    fcs_int doffset;       /* destination offset for the site index */
    fcs_int stride;        /* minimal contiguous block */  
    fcs_int skip;          /* gap between two strides (from the first element of one stride to the first elem. of next stride */
    fcs_int nblocks;       /* number of strides */
    fcs_int coord[2];      /* coordinates of the vector fields which has to be exchanged */
    fcs_int volume;        /* number of lattice sites in surface patch */
} t_surf_patch;

/** structure for properties of each lattice site */
typedef struct {
    fcs_float    charge;                     /* charge on site */
    fcs_float    permittivity[SPACE_DIM];    /* dielectric properties on adjoined bonds. */
    fcs_int       r[SPACE_DIM];               /* position of site in space */
} t_site;


/** \struct parameter_struct 
 Parameter structure. Contains global system information for MEMD algorithm
 */
typedef struct {
    /** = 1/c^2    speed of light parameter. */
    fcs_float f_mass;
    /** inverse of square root of f_mass. */
    fcs_float invsqrt_f_mass;
    /** prefactor to convert field to force. */
    fcs_float prefactor;
    /** prefactor / (4*pi) */
    fcs_float pref2;
    /** global back ground permittivity of the system */
    fcs_float permittivity;
    /** bjerrum length of the system. */
    fcs_float bjerrum;
    /** temperature of the system */
    fcs_float temperature;
    /** mesh size in one dimension */
    fcs_int    mesh;
    /** global box size */
    fcs_float box_length[SPACE_DIM];
    /** = 1/a = mesh / box_length */
    fcs_float inva;
    /** size of mesh cube */
    fcs_float a;
    /** self energy coefficients of the system */
    fcs_float alpha[8][8];
    /** time step of the simulation */
    fcs_float time_step;
    /** local number of particles on processor */
    fcs_int n_part;
    /** total number of particles (not used yet) */
    fcs_int n_part_total;
} parameter_struct;

typedef struct {
    MPI_Comm communicator;
    MPI_Comm original_comm;
    fcs_int size;
    fcs_int node_grid[SPACE_DIM];
    fcs_int this_node;
    fcs_int node_pos[SPACE_DIM];
    fcs_int my_left[SPACE_DIM];
    fcs_int my_right[SPACE_DIM];
    fcs_int node_neighbors[NDIRS];
} memd_mpi_parameters;

typedef struct {
    parameter_struct parameters;
    lattice_parameters lparams;
    memd_mpi_parameters mpiparams;
    fcs_int init_flag;
    t_site* lattice;
    fcs_float* Dfield;
    fcs_float* Bfield;
    /* create cell lists */
    memd_cell_list  local_cells;
    memd_cell_list  ghost_cells;
    t_dirs* neighbor;
    fcs_int total_energy_flag;
} memd_struct;


/* create system structure handle. Filled in maggs_set_parameters(); */
memd_struct     memd;


/*@}*/






/* MPI tags for the maggs communications: */
/* Used in maggs_init() -> calc_glue_patch(). */
#define REQ_MAGGS_SPREAD 300
#define REQ_MAGGS_EQUIL  301

/* Factors for self-influence currection */
/* (stemming from epsilon_zero and 4*pi) */
#define SELF_FACTOR_1 1.57364595
#define SELF_FACTOR_2 1.5078141
#define ROUND_ERROR_PREC 1.0e-14
// from lattice Green's function: 0.5054620197
// and 0.126365505
// e_0 = 8.85418781762 * 10^-12

/* Define numbers for directions and dimensions: */
#define NOWHERE -1                  /* not a direction */
#define XPLUS 0                     /* add for specific direction */
#define YPLUS 1
#define ZPLUS 2
#define ZMINUS 3
#define YMINUS 4
#define XMINUS 5
#define OPP_DIR(dir)	(5-(dir))	/* Opposite direction */


/* Three often used macros for looping over 3D */
#define FOR3D(dir) for(dir=0; dir<SPACE_DIM; dir++)

#define FORALL_INNER_SITES(i,j,k)					\
for(i=memd->lparams.inner_left_down[0];i<memd->lparams.inner_up_right[0];i++)	\
for(j=memd->lparams.inner_left_down[1];j<memd->lparams.inner_up_right[1];j++)	\
for(k=memd->lparams.inner_left_down[2];k<memd->lparams.inner_up_right[2];k++) 

#define FORALL_SITES(i,j,k)			\
for(i=0;i<memd->lparams.dim[0];i++)			\
for(j=0;j<memd->lparams.dim[1];j++)		\
for(k=0;k<memd->lparams.dim[2];k++) 
/* from ifndef MAGGS_H */

#define SQR(x) (x*x)

/* Mathematical constants, from gcc's math.h */
#ifndef M_PI
#define M_E		2.7182818284590452353602874713526625L  /* e */
#define M_LOG2E		1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10E	0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2		0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10	2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PI		3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2	1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4	0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PI	0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PI	0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPI	1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2	       	1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2	0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif

#endif
