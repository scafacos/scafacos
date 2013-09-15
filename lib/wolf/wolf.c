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
#include "wolf.h"


#if defined(FCS_ENABLE_DEBUG) || 0
# define DO_DEBUG
# define DEBUG_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define DEBUG_CMD(_cmd_)  Z_NOP()
#endif
#define DEBUG_PRINT_PREFIX  "WOLF_DEBUG: "

#if defined(FCS_ENABLE_INFO) || 0
# define DO_INFO
# define INFO_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define INFO_CMD(_cmd_)  Z_NOP()
#endif
#define INFO_PRINT_PREFIX  "WOLF_INFO: "

#if defined(FCS_ENABLE_TIMING) || 0
# define DO_TIMING
# define TIMING_CMD(_cmd_)  Z_MOP(_cmd_)
#else
# define TIMING_CMD(_cmd_)  Z_NOP()
#endif
#define TIMING_PRINT_PREFIX  "WOLF_TIMING: "

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


void ifcs_wolf_create(ifcs_wolf_t *wolf)
{
  wolf->box_base[0] = wolf->box_base[1] = wolf->box_base[2] = 0;
  wolf->box_a[0] = wolf->box_a[1] = wolf->box_a[2] = 0;
  wolf->box_b[0] = wolf->box_b[1] = wolf->box_b[2] = 0;
  wolf->box_c[0] = wolf->box_c[1] = wolf->box_c[2] = 0;
  wolf->periodicity[0] = wolf->periodicity[1] = wolf->periodicity[2] = -1;

  wolf->nparticles = wolf->max_nparticles = 0;
  wolf->positions = NULL;
  wolf->charges = NULL;
  wolf->field = NULL;
  wolf->potentials = NULL;

  wolf->cutoff = 0.0;
  wolf->alpha = 0.0;

  wolf->max_particle_move = -1;

  wolf->resort = 0;
  wolf->near_resort = FCS_NEAR_RESORT_NULL;
}


void ifcs_wolf_destroy(ifcs_wolf_t *wolf)
{
  fcs_near_resort_destroy(&wolf->near_resort);
}


void ifcs_wolf_set_system(ifcs_wolf_t *wolf, fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c, fcs_int *periodicity)
{
  fcs_int i;

  for (i = 0; i < 3; ++i)
  {
    wolf->box_base[i] = box_base[i];
    wolf->box_a[i] = box_a[i];
    wolf->box_b[i] = box_b[i];
    wolf->box_c[i] = box_c[i];

    if (periodicity) wolf->periodicity[i] = periodicity[i];
  }
}


void ifcs_wolf_set_particles(ifcs_wolf_t *wolf, fcs_int nparticles, fcs_int max_nparticles, fcs_float *positions, fcs_float *charges, fcs_float *field, fcs_float *potentials)
{
  wolf->nparticles = nparticles;
  wolf->max_nparticles = max_nparticles;
  wolf->positions = positions;
  wolf->charges = charges;
  wolf->field = field;
  wolf->potentials = potentials;
}


void ifcs_wolf_set_cutoff(ifcs_wolf_t *wolf, fcs_float cutoff)
{
  wolf->cutoff = cutoff;
}


void ifcs_wolf_get_cutoff(ifcs_wolf_t *wolf, fcs_float *cutoff)
{
  *cutoff = wolf->cutoff;
}


void ifcs_wolf_set_alpha(ifcs_wolf_t *wolf, fcs_float alpha)
{
  wolf->alpha = alpha;
}


void ifcs_wolf_get_alpha(ifcs_wolf_t *wolf, fcs_float *alpha)
{
  *alpha = wolf->alpha;
}


void ifcs_wolf_set_max_particle_move(ifcs_wolf_t *wolf, fcs_float max_particle_move)
{
  wolf->max_particle_move = max_particle_move;
}


void ifcs_wolf_set_resort(ifcs_wolf_t *wolf, fcs_int resort)
{
  wolf->resort = resort;
}


void ifcs_wolf_get_resort(ifcs_wolf_t *wolf, fcs_int *resort)
{
  *resort = wolf->resort;
}


void ifcs_wolf_get_resort_availability(ifcs_wolf_t *wolf, fcs_int *availability)
{
  *availability = fcs_near_resort_is_available(wolf->near_resort);
}


void ifcs_wolf_get_resort_particles(ifcs_wolf_t *wolf, fcs_int *resort_particles)
{
  if (wolf->near_resort == FCS_NEAR_RESORT_NULL)
  {
    *resort_particles = wolf->nparticles;
    return;
  }
  
  *resort_particles = fcs_near_resort_get_sorted_particles(wolf->near_resort);
}


#ifdef PRINT_PARTICLES
static void wolf_print_particles(fcs_int n, fcs_float *xyz, fcs_float *q, fcs_float *f, fcs_float *p)
{
  fcs_int i;

  for (i = 0; i < n; ++i)
  {
    printf("  %" FCS_LMOD_INT "d: [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f] [%" FCS_LMOD_FLOAT "f] -> [%.2" FCS_LMOD_FLOAT "f,%.2" FCS_LMOD_FLOAT "f,%.2" FCS_LMOD_FLOAT "f] [%.2" FCS_LMOD_FLOAT "f]\n",
      i, xyz[i * 3 + 0], xyz[i * 3 + 1], xyz[i * 3 + 2], q[i], f[i * 3 + 0], f[i * 3 + 1], f[i * 3 + 2], p[i]);
  }
}
#endif


/*static void wolf_virial(fcs_int n, fcs_float *xyz, fcs_float *q, fcs_float *f, fcs_float *v, int size, int rank, MPI_Comm comm)
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
}*/


typedef struct {
  fcs_float alpha, p_shift, f_shift;

} wolf_coulomb_field_potential_t;


static void wolf_coulomb_field_potential(const void *param, fcs_float dist, fcs_float *f, fcs_float *p)
{
  wolf_coulomb_field_potential_t *wcfp = (wolf_coulomb_field_potential_t *) param;

/*#if 1
  printf("alpha = %" FCS_LMOD_FLOAT "e\n", wcfp->alpha);
  printf("p_shift = %" FCS_LMOD_FLOAT "e\n", wcfp->p_shift);
  printf("f_shift = %" FCS_LMOD_FLOAT "e\n", wcfp->f_shift);
#endif*/

  fcs_float adist = wcfp->alpha * dist;
  fcs_float erfc_part_ri = erfc(adist) / dist;

  *p = erfc_part_ri - wcfp->p_shift;
  *f = -(erfc_part_ri + 2.0 * wcfp->alpha * 0.56418958354775627928034964498 * exp(-adist * adist)) / dist - wcfp->f_shift; /* FIXME: use fcs-type constant for 1/sqrt(pi) */
}

static FCS_NEAR_LOOP_FP(wolf_coulomb_loop_fp, wolf_coulomb_field_potential)


void ifcs_wolf_run(ifcs_wolf_t *wolf, MPI_Comm comm)
{
  fcs_int i;

  int comm_rank, comm_size;

  fcs_near_t near;
  wolf_coulomb_field_potential_t wcfp;

#ifdef DO_TIMING
  double t;
#endif


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  INFO_CMD(
    if (comm_rank == MASTER_RANK)
    {
      printf(INFO_PRINT_PREFIX "nparticles: %" FCS_LMOD_INT "d (max: %" FCS_LMOD_INT "d)\n", wolf->nparticles, wolf->max_nparticles);
      printf(INFO_PRINT_PREFIX "positions:  %p\n", wolf->positions);
      printf(INFO_PRINT_PREFIX "charges:    %p\n", wolf->charges);
      printf(INFO_PRINT_PREFIX "field:      %p\n", wolf->field);
      printf(INFO_PRINT_PREFIX "potentials: %p\n", wolf->potentials);
      printf(INFO_PRINT_PREFIX "periodicity: [%" FCS_LMOD_INT "d, %" FCS_LMOD_INT "d, %" FCS_LMOD_INT "d]\n", wolf->periodicity[0], wolf->periodicity[1], wolf->periodicity[2]);
      printf(INFO_PRINT_PREFIX "box_base: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", wolf->box_base[0], wolf->box_base[1], wolf->box_base[2]);
      printf(INFO_PRINT_PREFIX "box_a: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", wolf->box_a[0], wolf->box_a[1], wolf->box_a[2]);
      printf(INFO_PRINT_PREFIX "box_b: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", wolf->box_b[0], wolf->box_b[1], wolf->box_b[2]);
      printf(INFO_PRINT_PREFIX "box_c: [%" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f, %" FCS_LMOD_FLOAT "f]\n", wolf->box_c[0], wolf->box_c[1], wolf->box_c[2]);
      printf(INFO_PRINT_PREFIX "cutoff: %" FCS_LMOD_FLOAT "f\n", wolf->cutoff);
      printf(INFO_PRINT_PREFIX "alpha: %" FCS_LMOD_FLOAT "f\n", wolf->alpha);
    }
  );

  for (i = 0; i < wolf->nparticles; ++i) wolf->field[i * 3 + 0] = wolf->field[i * 3 + 1] = wolf->field[i * 3 + 2] = wolf->potentials[i] = 0.0;

#ifdef PRINT_PARTICLES
  printf("%d:   particles IN:\n", comm_rank);
  wolf_print_particles(wolf->nparticles, wolf->positions, wolf->charges, wolf->field, wolf->potentials);
#endif

  TIMING_SYNC(comm); TIMING_START(t);

  fcs_near_create(&near);

  fcs_near_set_loop(&near, wolf_coulomb_loop_fp);
  fcs_near_set_system(&near, wolf->box_base, wolf->box_a, wolf->box_b, wolf->box_c, wolf->periodicity);
  fcs_near_set_particles(&near, wolf->nparticles, wolf->max_nparticles, wolf->positions, wolf->charges, NULL, wolf->field, wolf->potentials);
  fcs_near_set_max_particle_move(&near, wolf->max_particle_move);
  fcs_near_set_resort(&near, wolf->resort);

  fcs_float acutoff = wolf->alpha * wolf->cutoff;
  fcs_float erfc_part_ri = erfc(acutoff) / wolf->cutoff;

  wcfp.alpha = wolf->alpha;
  wcfp.p_shift = erfc_part_ri;
  wcfp.f_shift = -(erfc_part_ri + 2.0 * wolf->alpha * 0.56418958354775627928034964498 * exp(-acutoff * acutoff)) / wolf->cutoff;

  fcs_near_field_solver(&near, wolf->cutoff, &wcfp, comm);

  if (wolf->resort)
  {
    fcs_near_resort_destroy(&wolf->near_resort);

    fcs_near_resort_create(&wolf->near_resort, &near);
/*    fcs_near_resort_print(wolf->near_resort, comm);*/
  }

  fcs_near_destroy(&near);

  TIMING_SYNC(comm); TIMING_STOP(t);

/*  wolf_virial(wolf->nparticles, wolf->positions, wolf->charges, wolf->field, wolf->virial, comm_size, comm_rank, comm);*/

#ifdef PRINT_PARTICLES
  printf("%d:   particles OUT:\n", comm_rank);
  wolf_print_particles(wolf->nparticles, wolf->positions, wolf->charges, wolf->field, wolf->potentials);
#endif

  TIMING_CMD(
    if (comm_rank == MASTER_RANK)
      printf(TIMING_PRINT_PREFIX "wolf: %f\n", t);
  );
}


void ifcs_wolf_resort_ints(ifcs_wolf_t *wolf, fcs_int *src, fcs_int *dst, fcs_int n, MPI_Comm comm)
{
  if (wolf->near_resort == FCS_NEAR_RESORT_NULL) return;
  
  fcs_near_resort_ints(wolf->near_resort, src, dst, n, comm);
}


void ifcs_wolf_resort_floats(ifcs_wolf_t *wolf, fcs_float *src, fcs_float *dst, fcs_int n, MPI_Comm comm)
{
  if (wolf->near_resort == FCS_NEAR_RESORT_NULL) return;
  
  fcs_near_resort_floats(wolf->near_resort, src, dst, n, comm);
}


void ifcs_wolf_resort_bytes(ifcs_wolf_t *wolf, void *src, void *dst, fcs_int n, MPI_Comm comm)
{
  if (wolf->near_resort == FCS_NEAR_RESORT_NULL) return;
  
  fcs_near_resort_bytes(wolf->near_resort, src, dst, n, comm);
}
