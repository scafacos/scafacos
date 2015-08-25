/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann
  
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
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "common/fcs-common/FCSCommon.h"

#include "common/gridsort/gridsort.h"
#include "common/near/near.h"

#include "z_tools.h"
#include "directc.h"


#if defined(FCS_ENABLE_DEBUG) || 0
# define DO_DEBUG
# define DEBUG_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DEBUG_CMD(_cmd_)  Z_NOP()
#endif
#define DEBUG_PRINT_PREFIX  "DIRECT_DEBUG: "

#if defined(FCS_ENABLE_INFO) || 0
# define DO_INFO
# define INFO_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define INFO_CMD(_cmd_)  Z_NOP()
#endif
#define INFO_PRINT_PREFIX  "DIRECT_INFO: "

#if defined(FCS_ENABLE_TIMING) || 0
# define DO_TIMING
# define TIMING_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define TIMING_CMD(_cmd_)  Z_NOP()
#endif
#define TIMING_PRINT_PREFIX  "DIRECT_TIMING: "

#define MASTER_RANK  0

/*#define PRINT_PARTICLES*/

#define DO_TIMING_SYNC

#ifdef DO_TIMING
# define TIMING_DECL(_decl_)       _decl_
# define TIMING_CMD(_cmd_)         Z_MOP(_cmd_)
#else
# define TIMING_DECL(_decl_)
# define TIMING_CMD(_cmd_)         Z_NOP()
#endif
#ifdef DO_TIMING_SYNC
# define TIMING_SYNC(_c_)          TIMING_CMD(MPI_Barrier(_c_);)
#else
# define TIMING_SYNC(_c_)          Z_NOP()
#endif
#define TIMING_START(_t_)          TIMING_CMD(((_t_) = MPI_Wtime());)
#define TIMING_STOP(_t_)           TIMING_CMD(((_t_) = MPI_Wtime() - (_t_));)
#define TIMING_STOP_ADD(_t_, _r_)  TIMING_CMD(((_r_) += MPI_Wtime() - (_t_));)


void fcs_directc_create(fcs_directc_t *directc)
{
  directc->box_base[0] = directc->box_base[1] = directc->box_base[2] = 0;
  directc->box_a[0] = directc->box_a[1] = directc->box_a[2] = 0;
  directc->box_b[0] = directc->box_b[1] = directc->box_b[2] = 0;
  directc->box_c[0] = directc->box_c[1] = directc->box_c[2] = 0;
  directc->periodicity[0] = directc->periodicity[1] = directc->periodicity[2] = -1;

  directc->nparticles = directc->max_nparticles = 0;
  directc->positions = NULL;
  directc->charges = NULL;
  directc->field = NULL;
  directc->potentials = NULL;

  directc->in_nparticles = 0;
  directc->in_positions = NULL;
  directc->in_charges = NULL;

/*  directc->out_nparticles = 0;
  directc->out_positions = NULL;
  directc->out_field = NULL;
  directc->out_potentials = NULL;*/

#if FCS_DIRECT_WITH_DIPOLES
  directc->dipole_nparticles = directc->max_dipole_nparticles = -1;
  directc->dipole_positions = NULL;
  directc->dipole_moments = NULL;
  directc->dipole_field = NULL;
  directc->dipole_potentials = NULL;
#endif

  directc->periodic_images[0] = directc->periodic_images[1] = directc->periodic_images[2] = 1;
  directc->cutoff = 0.0;
  directc->cutoff_with_near = 0;

  directc->max_particle_move = -1;

  directc->resort = 0;
  directc->near_resort = FCS_NEAR_RESORT_NULL;
}


void fcs_directc_destroy(fcs_directc_t *directc)
{
  fcs_near_resort_destroy(&directc->near_resort);
}


void fcs_directc_set_system(fcs_directc_t *directc, const fcs_float *box_base, const fcs_float *box_a, const fcs_float *box_b, const fcs_float *box_c, const fcs_int *periodicity)
{
  fcs_int i;

  for (i = 0; i < 3; ++i)
  {
    directc->box_base[i] = box_base[i];
    directc->box_a[i] = box_a[i];
    directc->box_b[i] = box_b[i];
    directc->box_c[i] = box_c[i];

    if (periodicity) directc->periodicity[i] = periodicity[i];
  }
}


void fcs_directc_set_particles(fcs_directc_t *directc, fcs_int nparticles, fcs_int max_nparticles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  directc->nparticles = nparticles;
  directc->max_nparticles = max_nparticles;
  directc->positions = positions;
  directc->charges = charges;
  directc->field = field;
  directc->potentials = potentials;
}


void fcs_directc_set_in_particles(fcs_directc_t *directc, fcs_int in_nparticles, fcs_float *in_positions, fcs_float *in_charges)
{
  directc->in_nparticles = in_nparticles;
  directc->in_positions = in_positions;
  directc->in_charges = in_charges;
}


/*void fcs_directc_set_out_particles(fcs_directc_t *directc, fcs_int out_nparticles, fcs_float *out_positions, fcs_float *out_field, fcs_float *out_potentials)
{
  directc->out_nparticles = out_nparticles;
  directc->out_positions = out_positions;
  directc->out_field = out_field;
  directc->out_potentials = out_potentials;
}*/


#if FCS_DIRECT_WITH_DIPOLES
void fcs_directc_set_dipole_particles(fcs_directc_t *directc, fcs_int dipole_nparticles, fcs_int max_dipole_nparticles, fcs_float *dipole_positions, fcs_float *dipole_moments, fcs_float *dipole_field, fcs_float *dipole_potentials)
{
  directc->dipole_nparticles = dipole_nparticles;
  directc->max_dipole_nparticles = max_dipole_nparticles;
  directc->dipole_positions = dipole_positions;
  directc->dipole_moments = dipole_moments;
  directc->dipole_field = dipole_field;
  directc->dipole_potentials = dipole_potentials;
}
#endif


void fcs_directc_set_periodic_images(fcs_directc_t *directc, fcs_int *periodic_images)
{
  fcs_int i;

  for (i = 0; i < 3; ++i)
    directc->periodic_images[i] = periodic_images[i];
}


void fcs_directc_get_periodic_images(fcs_directc_t *directc, fcs_int *periodic_images)
{
  fcs_int i;

  for (i = 0; i < 3; ++i)
    periodic_images[i] = directc->periodic_images[i];
}


void fcs_directc_set_cutoff(fcs_directc_t *directc, fcs_float cutoff)
{
  directc->cutoff = cutoff;
}


void fcs_directc_get_cutoff(fcs_directc_t *directc, fcs_float *cutoff)
{
  *cutoff = directc->cutoff;
}


void fcs_directc_set_cutoff_with_near(fcs_directc_t *directc, fcs_int cutoff_with_near)
{
  directc->cutoff_with_near = cutoff_with_near;
}


void fcs_directc_get_cutoff_with_near(fcs_directc_t *directc, fcs_int *cutoff_with_near)
{
  *cutoff_with_near = directc->cutoff_with_near;
}


void fcs_directc_set_max_particle_move(fcs_directc_t *directc, fcs_float max_particle_move)
{
  directc->max_particle_move = max_particle_move;
}


void fcs_directc_set_resort(fcs_directc_t *directc, fcs_int resort)
{
  directc->resort = resort;
}


void fcs_directc_get_resort(fcs_directc_t *directc, fcs_int *resort)
{
  *resort = directc->resort;
}


void fcs_directc_get_resort_availability(fcs_directc_t *directc, fcs_int *availability)
{
  *availability = fcs_near_resort_is_available(directc->near_resort);
}


void fcs_directc_get_resort_particles(fcs_directc_t *directc, fcs_int *resort_particles)
{
  if (directc->near_resort == FCS_NEAR_RESORT_NULL)
  {
    *resort_particles = directc->nparticles;
    return;
  }
  
  *resort_particles = fcs_near_resort_get_sorted_particles(directc->near_resort);
}


#ifdef PRINT_PARTICLES
static void directc_print_particles(fcs_int n, fcs_float *xyz, fcs_float *q, fcs_float *f, fcs_float *p)
{
  fcs_int i;

  for (i = 0; i < n; ++i)
  {
    printf("  %" FCS_LMOD_INT "d: [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f] [%" FCS_LMOD_FLOAT "f] -> [%.2" FCS_LMOD_FLOAT "f,%.2" FCS_LMOD_FLOAT "f,%.2" FCS_LMOD_FLOAT "f] [%.2" FCS_LMOD_FLOAT "f]\n",
      i, xyz[i * 3 + 0], xyz[i * 3 + 1], xyz[i * 3 + 2], q[i], f[i * 3 + 0], f[i * 3 + 1], f[i * 3 + 2], p[i]);
  }
}
#endif


#define DIRECT_COMPUTE_NEW  1

#if DIRECT_COMPUTE_NEW


typedef struct
{
  fcs_int nin;
  fcs_float *positions, *charges;

  fcs_int nout;
  fcs_float *field, *potentials;

} charges_t;


#if FCS_DIRECT_WITH_DIPOLES
typedef struct
{
  fcs_int nin;
  fcs_float *positions, *moments;

  fcs_int nout;
  fcs_float *field, *potentials;

} dipoles_t;
#endif


typedef struct
{
  fcs_int *periodic;
  fcs_float *box_a, *box_b, *box_c;

} periodics_t;


#define CHARGE_VS_CHARGE(_p_, _f_, _q_, _ir_, _ir3_, _dx_, _dy_, _dz_) Z_MOP( \
  (_p_)[0] += (_q_) * (_ir_); \
  (_f_)[0] += (_q_) * (_ir3_) * (_dx_); \
  (_f_)[1] += (_q_) * (_ir3_) * (_dy_); \
  (_f_)[2] += (_q_) * (_ir3_) * (_dz_);)


#if FCS_DIRECT_WITH_DIPOLES

#define DIPOLE_VS_DIPOLE(_p_, _f_, _m_, _ir_, _dx_, _dy_, _dz_) Z_MOP( \
  );

#define CHARGE_VS_DIPOLE(_p_, _f_, _m_, _ir_, _dx_, _dy_, _dz_) Z_MOP( \
  );

#define DIPOLE_VS_CHARGE(_p_, _f_, _q_, _ir_, _dx_, _dy_, _dz_) Z_MOP( \
  );

#endif


static void compute_charge_charge_loop(charges_t *charges, fcs_float cutoff)
{
  fcs_int i, j;
  fcs_float dx, dy, dz, ir, ir3;


  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (i = 0; i < charges->nout; ++i)
  {
    for (j = i + 1; j < charges->nout; ++j)
    {
      dx = charges->positions[i * 3 + 0] - charges->positions[j * 3 + 0];
      dy = charges->positions[i * 3 + 1] - charges->positions[j * 3 + 1];
      dz = charges->positions[i * 3 + 2] - charges->positions[j * 3 + 2];

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      ir3 = ir * ir * ir;

      CHARGE_VS_CHARGE(charges->potentials + i, charges->field + (i * 3), charges->charges[j], ir, ir3, dx, dy, dz);

      CHARGE_VS_CHARGE(charges->potentials + j, charges->field + (j * 3), charges->charges[i], ir, ir3, -dx, -dy, -dz);
    }

    for (j = charges->nout; j < charges->nin; ++j)
    {
      dx = charges->positions[i * 3 + 0] - charges->positions[j * 3 + 0];
      dy = charges->positions[i * 3 + 1] - charges->positions[j * 3 + 1];
      dz = charges->positions[i * 3 + 2] - charges->positions[j * 3 + 2];

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      ir3 = ir * ir * ir;

      CHARGE_VS_CHARGE(charges->potentials + i, charges->field + (i * 3), charges->charges[j], ir, ir3, dx, dy, dz);
    }
  }
}


static void compute_charge_from_charge_loop(charges_t *charges, charges_t *from_charges, periodics_t *periodics, fcs_int skip_origin, fcs_float cutoff)
{
  fcs_int i, j, pd[3];
  fcs_float dx, dy, dz, ir, ir3;

  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (pd[0] = -periodics->periodic[0]; pd[0] <= periodics->periodic[0]; ++pd[0])
  for (pd[1] = -periodics->periodic[1]; pd[1] <= periodics->periodic[1]; ++pd[1])
  for (pd[2] = -periodics->periodic[2]; pd[2] <= periodics->periodic[2]; ++pd[2])
  {
    if (skip_origin && pd[0] == 0 && pd[1] == 0 && pd[2] == 0) continue;

    for (i = 0; i < charges->nout; ++i)
    for (j = 0; j < from_charges->nin; ++j)
    {
      dx = charges->positions[i * 3 + 0] - (from_charges->positions[j * 3 + 0] + (pd[0] * periodics->box_a[0]) + (pd[1] * periodics->box_b[0]) + (pd[2] * periodics->box_c[0]));
      dy = charges->positions[i * 3 + 1] - (from_charges->positions[j * 3 + 1] + (pd[0] * periodics->box_a[1]) + (pd[1] * periodics->box_b[1]) + (pd[2] * periodics->box_c[1]));
      dz = charges->positions[i * 3 + 2] - (from_charges->positions[j * 3 + 2] + (pd[0] * periodics->box_a[2]) + (pd[1] * periodics->box_b[2]) + (pd[2] * periodics->box_c[2]));

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      ir3 = ir * ir * ir;

      CHARGE_VS_CHARGE(charges->potentials + i, charges->field + (i * 3), from_charges->charges[j], ir, ir3, dx, dy, dz);
    }
  }
}


#if FCS_DIRECT_WITH_DIPOLES

static void compute_dipole_dipole_loop(dipoles_t *dipoles, fcs_float cutoff)
{
  fcs_int i, j;
  fcs_float dx, dy, dz, ir;


  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (i = 0; i < dipoles->nout; ++i)
  {
    for (j = i + 1; j < dipoles->nout; ++j)
    {
      dx = dipoles->positions[i * 3 + 0] - dipoles->positions[j * 3 + 0];
      dy = dipoles->positions[i * 3 + 1] - dipoles->positions[j * 3 + 1];
      dz = dipoles->positions[i * 3 + 2] - dipoles->positions[j * 3 + 2];

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      DIPOLE_VS_DIPOLE(dipoles->potentials + (i * 3), dipoles->field + (i * 6), dipoles->moments + (j * 3), ir, dx, dy, dz);

      DIPOLE_VS_DIPOLE(dipoles->potentials + (j * 3), dipoles->field + (j * 6), dipoles->moments + (i * 3), ir, -dx, -dy, -dz);
    }

    for (j = dipoles->nout; j < dipoles->nin; ++j)
    {
      dx = dipoles->positions[i * 3 + 0] - dipoles->positions[j * 3 + 0];
      dy = dipoles->positions[i * 3 + 1] - dipoles->positions[j * 3 + 1];
      dz = dipoles->positions[i * 3 + 2] - dipoles->positions[j * 3 + 2];

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      DIPOLE_VS_DIPOLE(dipoles->potentials + (i * 3), dipoles->field + (i * 6), dipoles->moments + (j * 3), ir, dx, dy, dz);
    }
  }
}


static void compute_dipole_from_dipole_loop(dipoles_t *dipoles, dipoles_t *from_dipoles, periodics_t *periodics, fcs_int skip_origin, fcs_float cutoff)
{
  fcs_int i, j, pd[3];
  fcs_float dx, dy, dz, ir;

  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (pd[0] = -periodics->periodic[0]; pd[0] <= periodics->periodic[0]; ++pd[0])
  for (pd[1] = -periodics->periodic[1]; pd[1] <= periodics->periodic[1]; ++pd[1])
  for (pd[2] = -periodics->periodic[2]; pd[2] <= periodics->periodic[2]; ++pd[2])
  {
    if (skip_origin && pd[0] == 0 && pd[1] == 0 && pd[2] == 0) continue;

    for (i = 0; i < dipoles->nout; ++i)
    for (j = 0; j < from_dipoles->nin; ++j)
    {
      dx = dipoles->positions[i * 3 + 0] - (from_dipoles->positions[j * 3 + 0] + (pd[0] * periodics->box_a[0]) + (pd[1] * periodics->box_b[0]) + (pd[2] * periodics->box_c[0]));
      dy = dipoles->positions[i * 3 + 1] - (from_dipoles->positions[j * 3 + 1] + (pd[0] * periodics->box_a[1]) + (pd[1] * periodics->box_b[1]) + (pd[2] * periodics->box_c[1]));
      dz = dipoles->positions[i * 3 + 2] - (from_dipoles->positions[j * 3 + 2] + (pd[0] * periodics->box_a[2]) + (pd[1] * periodics->box_b[2]) + (pd[2] * periodics->box_c[2]));

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      DIPOLE_VS_DIPOLE(dipoles->potentials + (i * 3), dipoles->field + (i * 6), from_dipoles->moments + (j * 3), ir, dx, dy, dz);
    }
  }
}


static void compute_charge_from_dipole_loop(charges_t *charges, dipoles_t *dipoles, periodics_t *periodics, fcs_float cutoff)
{
  fcs_int i, j, pd[3];
  fcs_float dx, dy, dz, ir;

  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (pd[0] = -periodics->periodic[0]; pd[0] <= periodics->periodic[0]; ++pd[0])
  for (pd[1] = -periodics->periodic[1]; pd[1] <= periodics->periodic[1]; ++pd[1])
  for (pd[2] = -periodics->periodic[2]; pd[2] <= periodics->periodic[2]; ++pd[2])
  {
    for (i = 0; i < charges->nout; ++i)
    for (j = 0; j < dipoles->nin; ++j)
    {
      dx = charges->positions[i * 3 + 0] - (dipoles->positions[j * 3 + 0] + (pd[0] * periodics->box_a[0]) + (pd[1] * periodics->box_b[0]) + (pd[2] * periodics->box_c[0]));
      dy = charges->positions[i * 3 + 1] - (dipoles->positions[j * 3 + 1] + (pd[0] * periodics->box_a[1]) + (pd[1] * periodics->box_b[1]) + (pd[2] * periodics->box_c[1]));
      dz = charges->positions[i * 3 + 2] - (dipoles->positions[j * 3 + 2] + (pd[0] * periodics->box_a[2]) + (pd[1] * periodics->box_b[2]) + (pd[2] * periodics->box_c[2]));

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      CHARGE_VS_DIPOLE(charges->potentials + i, charges->field + (i * 3), dipoles->moments + (j * 3), ir, dx, dy, dz);
    }
  }
}


static void compute_dipole_from_charge_loop(dipoles_t *dipoles, charges_t *charges, periodics_t *periodics, fcs_float cutoff)
{
  fcs_int i, j, pd[3];
  fcs_float dx, dy, dz, ir;

  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (pd[0] = -periodics->periodic[0]; pd[0] <= periodics->periodic[0]; ++pd[0])
  for (pd[1] = -periodics->periodic[1]; pd[1] <= periodics->periodic[1]; ++pd[1])
  for (pd[2] = -periodics->periodic[2]; pd[2] <= periodics->periodic[2]; ++pd[2])
  {
    for (i = 0; i < dipoles->nout; ++i)
    for (j = 0; j < charges->nin; ++j)
    {
      dx = dipoles->positions[i * 3 + 0] - (charges->positions[j * 3 + 0] + (pd[0] * periodics->box_a[0]) + (pd[1] * periodics->box_b[0]) + (pd[2] * periodics->box_c[0]));
      dy = dipoles->positions[i * 3 + 1] - (charges->positions[j * 3 + 1] + (pd[0] * periodics->box_a[1]) + (pd[1] * periodics->box_b[1]) + (pd[2] * periodics->box_c[1]));
      dz = dipoles->positions[i * 3 + 2] - (charges->positions[j * 3 + 2] + (pd[0] * periodics->box_a[2]) + (pd[1] * periodics->box_b[2]) + (pd[2] * periodics->box_c[2]));

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      DIPOLE_VS_CHARGE(dipoles->potentials + (i * 3), dipoles->field + (i * 6), charges->charges[j], ir, dx, dy, dz);
    }
  }
}

#endif


static void compute_charge_charge(charges_t *charges, periodics_t *periodics, fcs_float cutoff)
{
  compute_charge_charge_loop(charges, cutoff);

  compute_charge_from_charge_loop(charges, charges, periodics, 1, cutoff);
}


#if FCS_DIRECT_WITH_DIPOLES

static void compute_dipole_dipole(dipoles_t *dipoles, periodics_t *periodics, fcs_float cutoff)
{
  compute_dipole_dipole_loop(dipoles, cutoff);

  compute_dipole_from_dipole_loop(dipoles, dipoles, periodics, 1, cutoff);
}


static void compute_charge_dipole(charges_t *charges, dipoles_t *dipoles, periodics_t *periodics, fcs_float cutoff)
{
  compute_charge_from_dipole_loop(charges, dipoles, periodics, cutoff);
  compute_dipole_from_charge_loop(dipoles, charges, periodics, cutoff);
}

#endif


static void compute_charge_from_charge(charges_t *charges, charges_t *from_charges, periodics_t *periodics, fcs_float cutoff)
{
  compute_charge_from_charge_loop(charges, from_charges, periodics, 0, cutoff);
}


#if FCS_DIRECT_WITH_DIPOLES

static void compute_dipole_from_dipole(dipoles_t *dipoles, dipoles_t *from_dipoles, periodics_t *periodics, fcs_float cutoff)
{
  compute_dipole_from_dipole_loop(dipoles, from_dipoles, periodics, 0, cutoff);
}


static void compute_charge_from_dipole(charges_t *charges, dipoles_t *dipoles, periodics_t *periodics, fcs_float cutoff)
{
  compute_charge_from_dipole_loop(charges, dipoles, periodics, cutoff);
}


static void compute_dipole_from_charge(dipoles_t *dipoles, charges_t *charges, periodics_t *periodics, fcs_float cutoff)
{
  compute_dipole_from_charge_loop(dipoles, charges, periodics, cutoff);
}

#endif


static void directc_compute(fcs_directc_t *directc, fcs_int *periodic, int size, int rank, MPI_Comm comm)
{
  fcs_int l;

  periodics_t periodics;

  charges_t my_charges, other_charges;
#if FCS_DIRECT_WITH_DIPOLES
  dipoles_t my_dipoles, other_dipoles;
#endif

  periodics.periodic = periodic;
  periodics.box_a = directc->box_a;
  periodics.box_b = directc->box_b;
  periodics.box_c = directc->box_c;

  my_charges.nin = directc->nparticles;
  my_charges.positions = directc->positions;
  my_charges.charges = directc->charges;
  my_charges.nout = directc->nparticles;
  my_charges.field = directc->field;
  my_charges.potentials = directc->potentials;

  other_charges.nin = directc->in_nparticles;
  other_charges.positions = directc->in_positions;
  other_charges.charges = directc->in_charges;
  other_charges.nout = 0;
  other_charges.field = NULL;
  other_charges.potentials = NULL;

  /* vs. myself */
  compute_charge_charge(&my_charges, &periodics, directc->cutoff);
  compute_charge_from_charge(&my_charges, &other_charges, &periodics, directc->cutoff);

#if FCS_DIRECT_WITH_DIPOLES
  my_dipoles.nin = directc->dipole_nparticles;
  my_dipoles.positions = directc->dipole_positions;
  my_dipoles.moments = directc->dipole_moments;
  my_dipoles.nout = directc->dipole_nparticles;
  my_dipoles.field = directc->dipole_field;
  my_dipoles.potentials = directc->dipole_potentials;

  other_dipoles.nin = 0;
  other_dipoles.positions = NULL;
  other_dipoles.moments = NULL;
  other_dipoles.nout = 0;
  other_dipoles.field = NULL;
  other_dipoles.potentials = NULL;

  compute_dipole_dipole(&my_dipoles, &periodics, directc->cutoff);
  compute_charge_dipole(&my_charges, &my_dipoles, &periodics, directc->cutoff);
  compute_dipole_from_charge(&my_dipoles, &other_charges, &periodics, directc->cutoff);
#endif

  /* vs. others */
  if (size > 1)
  {
    fcs_int my_n[3], all_n[3 * size], other_max_nfloats;

    my_n[0] = directc->nparticles;
    my_n[1] = directc->in_nparticles;
    my_n[2] =
#if FCS_DIRECT_WITH_DIPOLES
      directc->dipole_nparticles;
#else
      0;
#endif
    MPI_Allgather(my_n, 3, FCS_MPI_INT, all_n, 3, FCS_MPI_INT, comm);

    other_max_nfloats = 0;
    for (l = 0; l < size; ++l)
    {
      if (l == rank) continue;
      other_max_nfloats = z_max(other_max_nfloats, (all_n[3 * l + 0] + all_n[3 * l + 1]) * 4 + all_n[3 * l + 2] * 6);
    }

    fcs_float *other_floats = malloc(other_max_nfloats * sizeof(fcs_float));

    int prev_rank = (rank - 1 + size) % size;
    int next_rank = (rank + 1) % size;
    fcs_int *prev_n = NULL;

    for (l = 1; l < size; ++l)
    {
      int other_rank = (rank - l + size) % size;
      fcs_int *other_n = all_n + 3 * other_rank;

      fcs_float *f = other_floats;

      other_charges.nin = other_n[0] + other_n[1];
      other_charges.positions = f;
      f += 3 * other_charges.nin;
      other_charges.charges = f;
      f += 1 * other_charges.nin;

#if FCS_DIRECT_WITH_DIPOLES
      other_dipoles.nin = other_n[2];
      other_dipoles.positions = f;
      f += 3 * other_dipoles.nin;
      other_dipoles.moments = f;
      f += 3 * other_dipoles.nin;
#endif

      if (l == 1)
      {
        /* first round */
        MPI_Sendrecv(directc->positions, 3 * my_n[0], FCS_MPI_FLOAT, next_rank, 0, other_charges.positions, 3 * other_n[0], FCS_MPI_FLOAT, prev_rank, 0, comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(directc->in_positions, 3 * my_n[1], FCS_MPI_FLOAT, next_rank, 0, other_charges.positions + 3 * other_n[0], 3 * other_n[1], FCS_MPI_FLOAT, prev_rank, 0, comm, MPI_STATUS_IGNORE);

        MPI_Sendrecv(directc->charges, my_n[0], FCS_MPI_FLOAT, next_rank, 0, other_charges.charges, other_n[0], FCS_MPI_FLOAT, prev_rank, 0, comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(directc->in_charges, my_n[1], FCS_MPI_FLOAT, next_rank, 0, other_charges.charges + other_n[0], other_n[1], FCS_MPI_FLOAT, prev_rank, 0, comm, MPI_STATUS_IGNORE);

#if FCS_DIRECT_WITH_DIPOLES
        MPI_Sendrecv(directc->dipole_positions, 3 * my_n[2], FCS_MPI_FLOAT, next_rank, 0, other_dipoles.positions, 3 * other_n[2], FCS_MPI_FLOAT, prev_rank, 0, comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(directc->dipole_moments, 3 * my_n[2], FCS_MPI_FLOAT, next_rank, 0, other_dipoles.moments, 3 * other_n[2], FCS_MPI_FLOAT, prev_rank, 0, comm, MPI_STATUS_IGNORE);
#endif
      } else
      {
        /* further rounds */
        fcs_int prev_nfloats = (prev_n[0] + prev_n[1]) * 4 + prev_n[2] * 6;
        fcs_int other_nfloats = (other_n[0] + other_n[1]) * 4 + other_n[2] * 6;
        fcs_int replace_nfloats = z_min(prev_nfloats, other_nfloats);

        MPI_Sendrecv_replace(other_floats, replace_nfloats, FCS_MPI_FLOAT, next_rank, 0, prev_rank, 0, comm, MPI_STATUS_IGNORE);

        if (replace_nfloats < prev_nfloats) MPI_Send(other_floats + replace_nfloats, prev_nfloats - replace_nfloats, FCS_MPI_FLOAT, next_rank, 0, comm);
        else if (replace_nfloats < other_nfloats) MPI_Recv(other_floats + replace_nfloats, other_nfloats - replace_nfloats, FCS_MPI_FLOAT, prev_rank, 0, comm, MPI_STATUS_IGNORE);
      }

      compute_charge_from_charge(&my_charges, &other_charges, &periodics, directc->cutoff);
#if FCS_DIRECT_WITH_DIPOLES
      compute_dipole_from_dipole(&my_dipoles, &other_dipoles, &periodics, directc->cutoff);
      compute_charge_from_dipole(&my_charges, &other_dipoles, &periodics, directc->cutoff);
      compute_dipole_from_charge(&my_dipoles, &other_charges, &periodics, directc->cutoff);
#endif

      prev_n = other_n;
    }

    free(other_floats);
  }
}


#else /* DIRECT_COMPUTE_NEW */


static void directc_local_one(fcs_int nout, fcs_int nin, fcs_float *xyz, fcs_float *q, fcs_float *f, fcs_float *p, fcs_float cutoff)
{
  fcs_int i, j;
  fcs_float dx, dy, dz, ir;


  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (i = 0; i < nout; ++i)
  {
    for (j = i + 1; j < nout; ++j)
    {
      dx = xyz[i*3+0] - xyz[j*3+0];
      dy = xyz[i*3+1] - xyz[j*3+1];
      dz = xyz[i*3+2] - xyz[j*3+2];

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      p[i] += q[j] * ir;
      p[j] += q[i] * ir;

      f[i*3+0] += q[j] * dx * ir * ir * ir;
      f[i*3+1] += q[j] * dy * ir * ir * ir;
      f[i*3+2] += q[j] * dz * ir * ir * ir;

      f[j*3+0] -= q[i] * dx * ir * ir * ir;
      f[j*3+1] -= q[i] * dy * ir * ir * ir;
      f[j*3+2] -= q[i] * dz * ir * ir * ir;
    }
    
    for (j = nout; j < nin; ++j)
    {
      dx = xyz[i*3+0] - xyz[j*3+0];
      dy = xyz[i*3+1] - xyz[j*3+1];
      dz = xyz[i*3+2] - xyz[j*3+2];

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      p[i] += q[j] * ir;

      f[i*3+0] += q[j] * dx * ir * ir * ir;
      f[i*3+1] += q[j] * dy * ir * ir * ir;
      f[i*3+2] += q[j] * dz * ir * ir * ir;
    }
  }
}


static void directc_local_two(fcs_int n0, fcs_float *xyz0, fcs_float *q0, fcs_int n1, fcs_float *xyz1, fcs_float *q1, fcs_float *f, fcs_float *p, fcs_float cutoff)
{
  fcs_int i, j;
  fcs_float dx, dy, dz, ir;


  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (i = 0; i < n0; ++i)
  for (j = 0; j < n1; ++j)
  {
    dx = xyz0[i*3+0] - xyz1[j*3+0];
    dy = xyz0[i*3+1] - xyz1[j*3+1];
    dz = xyz0[i*3+2] - xyz1[j*3+2];
  
    ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

    if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

    p[i] += q1[j] * ir;
    
    f[i*3+0] += q1[j] * dx * ir * ir * ir;
    f[i*3+1] += q1[j] * dy * ir * ir * ir;
    f[i*3+2] += q1[j] * dz * ir * ir * ir;
  }
}


static void directc_local_periodic(fcs_int n0, fcs_float *xyz0, fcs_float *q0, fcs_int n1, fcs_float *xyz1, fcs_float *q1, fcs_float *f, fcs_float *p, fcs_int *periodic, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_float cutoff)
{
  fcs_int i, j, pd[3];
  fcs_float dx, dy, dz, ir;


  if (fcs_fabs(cutoff) > 0) cutoff = 1.0 / cutoff;

  for (pd[0] = -periodic[0]; pd[0] <= periodic[0]; ++pd[0])
  for (pd[1] = -periodic[1]; pd[1] <= periodic[1]; ++pd[1])
  for (pd[2] = -periodic[2]; pd[2] <= periodic[2]; ++pd[2])
  {
    if (pd[0] == 0 && pd[1] == 0 && pd[2] == 0) continue;

    for (i = 0; i < n0; ++i)
    for (j = 0; j < n1; ++j)
    {
      dx = xyz0[i*3+0] - (xyz1[j*3+0] + (pd[0] * box_a[0]) + (pd[1] * box_b[0]) + (pd[2] * box_c[0]));
      dy = xyz0[i*3+1] - (xyz1[j*3+1] + (pd[0] * box_a[1]) + (pd[1] * box_b[1]) + (pd[2] * box_c[1]));
      dz = xyz0[i*3+2] - (xyz1[j*3+2] + (pd[0] * box_a[2]) + (pd[1] * box_b[2]) + (pd[2] * box_c[2]));

      ir = 1.0 / fcs_sqrt(z_sqr(dx) + z_sqr(dy) + z_sqr(dz));

      if ((cutoff > 0 && cutoff > ir) || (cutoff < 0 && -cutoff < ir)) continue;

      p[i] += q1[j] * ir;

      f[i*3+0] += q1[j] * dx * ir * ir * ir;
      f[i*3+1] += q1[j] * dy * ir * ir * ir;
      f[i*3+2] += q1[j] * dz * ir * ir * ir;
    }
  }
}


static void directc_global(fcs_directc_t *directc, fcs_int *periodic, int size, int rank, MPI_Comm comm)
{
  fcs_int l;

  fcs_int my_n, max_n, all_n[size], other_n, other_n_next;

  fcs_float *other_xyzq, *other_xyz, *other_q;

  MPI_Status status;


  my_n = directc->nparticles + directc->in_nparticles;
  MPI_Allreduce(&my_n, &max_n, 1, FCS_MPI_INT, MPI_MAX, comm);
  MPI_Allgather(&my_n, 1, FCS_MPI_INT, all_n, 1, FCS_MPI_INT, comm);

  other_xyzq = calloc(max_n, 4*sizeof(fcs_float));

  other_n_next = all_n[rank];

  other_n = other_n_next;
  other_xyz = other_xyzq;
  other_q = other_xyzq + 3 * other_n;

  memcpy(other_xyz, directc->positions, directc->nparticles * 3 * sizeof(fcs_float));
  memcpy(other_q, directc->charges, directc->nparticles * sizeof(fcs_float));

  if (directc->in_nparticles > 0 && directc->in_positions && directc->in_charges)
  {
    memcpy(other_xyz + directc->nparticles * 3, directc->in_positions, directc->in_nparticles * 3 * sizeof(fcs_float));
    memcpy(other_q + directc->nparticles, directc->in_charges, directc->in_nparticles * sizeof(fcs_float));
  }

  directc_local_one(directc->nparticles, other_n, other_xyz, other_q, directc->field, directc->potentials, directc->cutoff);
  directc_local_periodic(directc->nparticles, directc->positions, directc->charges, other_n, other_xyz, other_q, directc->field, directc->potentials, periodic, directc->box_a, directc->box_b, directc->box_c, directc->cutoff);

  for (l = 1; l < size; ++l)
  {
    other_n_next = all_n[(rank - l + size) % size];
    
    MPI_Sendrecv_replace(other_xyzq, max_n * (3 + 1), FCS_MPI_FLOAT, (rank + 1) % size, 0, (rank - 1 + size) % size, 0, comm, &status);

    other_n = other_n_next;
    other_xyz = other_xyzq;
    other_q = other_xyzq + 3 * other_n;

    directc_local_two(directc->nparticles, directc->positions, directc->charges, other_n, other_xyz, other_q, directc->field, directc->potentials, directc->cutoff);
    directc_local_periodic(directc->nparticles, directc->positions, directc->charges, other_n, other_xyz, other_q, directc->field, directc->potentials, periodic, directc->box_a, directc->box_b, directc->box_c, directc->cutoff);
  }

  free(other_xyzq);
}


#endif /* DIRECT_COMPUTE_NEW */


static void directc_virial(fcs_int n, fcs_float *xyz, fcs_float *q, fcs_float *f, fcs_float *v, int size, int rank, MPI_Comm comm)
{
  fcs_int i;
  fcs_float my_v[9];


  for (i = 0; i < 9; ++i) my_v[i] = 0;

  for (i = 0; i < n; ++i)
  {
    my_v[0] += f[i*3+0] * q[i] * xyz[i*3+0];
    my_v[1] += f[i*3+0] * q[i] * xyz[i*3+1];
    my_v[2] += f[i*3+0] * q[i] * xyz[i*3+2];
    my_v[3] += f[i*3+1] * q[i] * xyz[i*3+0];
    my_v[4] += f[i*3+1] * q[i] * xyz[i*3+1];
    my_v[5] += f[i*3+1] * q[i] * xyz[i*3+2];
    my_v[6] += f[i*3+2] * q[i] * xyz[i*3+0];
    my_v[7] += f[i*3+2] * q[i] * xyz[i*3+1];
    my_v[8] += f[i*3+2] * q[i] * xyz[i*3+2];
  }

  MPI_Allreduce(my_v, v, 9, FCS_MPI_FLOAT, MPI_SUM, comm);
}


#define VSIZE(_v_)  fcs_sqrt(z_sqr((_v_)[0]) + z_sqr((_v_)[1]) + z_sqr((_v_)[2]))

static fcs_float get_periodic_factor(fcs_float *v0, fcs_float *v1, fcs_float *v2, fcs_float cutoff)
{
  fcs_float n[3], f;


  n[0] = v1[1] * v2[2] - v1[2] * v2[1];
  n[1] = v1[2] * v2[0] - v1[0] * v2[2];
  n[2] = v1[0] * v2[1] - v1[1] * v2[0];

  f = VSIZE(n) * cutoff / (n[0] * v0[0] + n[1] * v0[1] + n[2] * v0[2]);

  if (f < 0) f *= -1.0;

  return f;
}


static void directc_coulomb_field_potential(const void *param, fcs_float dist, fcs_float *f, fcs_float *p)
{
  *p = 1.0 / dist;
  *f = -(*p) * (*p);
}

static FCS_NEAR_LOOP_FP(directc_coulomb_loop_fp, directc_coulomb_field_potential)


void fcs_directc_run(fcs_directc_t *directc, MPI_Comm comm)
{
  fcs_int i;

  int comm_rank, comm_size;

  fcs_near_t near;
  fcs_int periodic[3] = { 0, 0, 0 };

#ifdef DO_TIMING
  double t;
#endif


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  INFO_CMD(
    if (comm_rank == MASTER_RANK)
    {
      printf(INFO_PRINT_PREFIX "periodicity: [%" FCS_LMOD_INT "d, %" FCS_LMOD_INT "d, %" FCS_LMOD_INT "d]\n", directc->periodicity[0], directc->periodicity[1], directc->periodicity[2]);
      printf(INFO_PRINT_PREFIX "box_base: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", directc->box_base[0], directc->box_base[1], directc->box_base[2]);
      printf(INFO_PRINT_PREFIX "box_a: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", directc->box_a[0], directc->box_a[1], directc->box_a[2]);
      printf(INFO_PRINT_PREFIX "box_b: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", directc->box_b[0], directc->box_b[1], directc->box_b[2]);
      printf(INFO_PRINT_PREFIX "box_c: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", directc->box_c[0], directc->box_c[1], directc->box_c[2]);
      printf(INFO_PRINT_PREFIX "cutoff: %" FCS_LMOD_FLOAT "f, with near: %" FCS_LMOD_INT "d\n", directc->cutoff, directc->cutoff_with_near);
      printf(INFO_PRINT_PREFIX "nparticles: %" FCS_LMOD_INT "d (max: %" FCS_LMOD_INT "d)\n", directc->nparticles, directc->max_nparticles);
      printf(INFO_PRINT_PREFIX "positions:  %p\n", directc->positions);
      printf(INFO_PRINT_PREFIX "charges:    %p\n", directc->charges);
      printf(INFO_PRINT_PREFIX "field:      %p\n", directc->field);
      printf(INFO_PRINT_PREFIX "potentials: %p\n", directc->potentials);
      printf(INFO_PRINT_PREFIX "dipole_nparticles: %" FCS_LMOD_INT "d (max: %" FCS_LMOD_INT "d)\n", directc->dipole_nparticles, directc->max_dipole_nparticles);
      printf(INFO_PRINT_PREFIX "dipole_positions:  %p\n", directc->dipole_positions);
      printf(INFO_PRINT_PREFIX "dipole_moments:    %p\n", directc->dipole_moments);
      printf(INFO_PRINT_PREFIX "dipole_field:      %p\n", directc->dipole_field);
      printf(INFO_PRINT_PREFIX "dipole_potentials: %p\n", directc->dipole_potentials);
    }
  );

  periodic[0] = directc->periodicity[0] * directc->periodic_images[0];
  periodic[1] = directc->periodicity[1] * directc->periodic_images[1];
  periodic[2] = directc->periodicity[2] * directc->periodic_images[2];

  if (directc->cutoff > 0.0)
  {
    if (directc->periodicity[0]) periodic[0] = (fcs_int) fcs_ceil(get_periodic_factor(directc->box_a, directc->box_b, directc->box_c, directc->cutoff));
    if (directc->periodicity[1]) periodic[1] = (fcs_int) fcs_ceil(get_periodic_factor(directc->box_b, directc->box_c, directc->box_a, directc->cutoff));
    if (directc->periodicity[2]) periodic[2] = (fcs_int) fcs_ceil(get_periodic_factor(directc->box_c, directc->box_a, directc->box_b, directc->cutoff));
  }

  INFO_CMD(
    if (comm_rank == MASTER_RANK) printf(INFO_PRINT_PREFIX "periodic: %" FCS_LMOD_INT "d / %" FCS_LMOD_INT "d / %" FCS_LMOD_INT "d\n", periodic[0], periodic[1], periodic[2]);
  );

  for (i = 0; i < directc->nparticles; ++i) directc->field[i * 3 + 0] = directc->field[i * 3 + 1] = directc->field[i * 3 + 2] = directc->potentials[i] = 0.0;

#ifdef PRINT_PARTICLES
  printf("%d:   particles IN:\n", comm_rank);
  directc_print_particles(directc->nparticles, directc->positions, directc->charges, directc->field, directc->potentials);
#endif

  TIMING_SYNC(comm); TIMING_START(t);

  if (directc->cutoff_with_near)
  {
    fcs_near_create(&near);

    fcs_near_set_loop(&near, directc_coulomb_loop_fp);
    fcs_near_set_system(&near, directc->box_base, directc->box_a, directc->box_b, directc->box_c, periodic);
    fcs_near_set_particles(&near, directc->nparticles, directc->max_nparticles, directc->positions, directc->charges, NULL, directc->field, directc->potentials);
    fcs_near_set_max_particle_move(&near, directc->max_particle_move);
    fcs_near_set_resort(&near, directc->resort);

    fcs_near_field_solver(&near, fabs(directc->cutoff), NULL, comm);

    if (directc->resort)
    {
      fcs_near_resort_destroy(&directc->near_resort);

      fcs_near_resort_create(&directc->near_resort, &near);
/*      fcs_near_resort_print(directc->near_resort, comm);*/
    }

    fcs_near_destroy(&near);

  } else
  {
#if DIRECT_COMPUTE_NEW
    directc_compute(directc, periodic, comm_size, comm_rank, comm);
#else
    directc_global(directc, periodic, comm_size, comm_rank, comm);
#endif
  }

  TIMING_SYNC(comm); TIMING_STOP(t);

  directc_virial(directc->nparticles, directc->positions, directc->charges, directc->field, directc->virial, comm_size, comm_rank, comm);

#ifdef PRINT_PARTICLES
  printf("%d:   particles OUT:\n", comm_rank);
  directc_print_particles(directc->nparticles, directc->positions, directc->charges, directc->field, directc->potentials);
#endif

  TIMING_CMD(
    if (comm_rank == MASTER_RANK)
      printf(TIMING_PRINT_PREFIX "directc: %f\n", t);
  );
}


void fcs_directc_resort_ints(fcs_directc_t *directc, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  if (directc->near_resort == FCS_NEAR_RESORT_NULL) return;
  
  fcs_near_resort_ints(directc->near_resort, src, dst, n, comm);
}


void fcs_directc_resort_floats(fcs_directc_t *directc, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  if (directc->near_resort == FCS_NEAR_RESORT_NULL) return;
  
  fcs_near_resort_floats(directc->near_resort, src, dst, n, comm);
}


void fcs_directc_resort_bytes(fcs_directc_t *directc, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  if (directc->near_resort == FCS_NEAR_RESORT_NULL) return;
  
  fcs_near_resort_bytes(directc->near_resort, src, dst, n, comm);
}
