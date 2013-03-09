
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "fcs.h"

#include "common/fcs-common/FCSCommon.h"

#include "common/near/near.h"

#ifndef M_PI
# define M_PI  3.1415926535897932384626433832795029L  /* pi */
#endif


#define ASSERT_FCS(_r_) \
  do { \
    if(_r_) { \
      fcsResult_printResult(_r_); MPI_Finalize(); exit(-1); \
    } \
  } while (0)


fcs_float coulomb_potential(const void *param, fcs_float dist)
{
  return 1.0 / dist;
}


fcs_float coulomb_field(const void *param, fcs_float dist)
{
  return -1.0 / (dist * dist);
}


void coulomb_field_potential(const void *param, fcs_float dist, fcs_float *f, fcs_float *p)
{
  *p = 1.0 / dist;
  *f = -(*p) * (*p);
}


fcs_float coulomb_potential_3diff(const void *param, fcs_float dist, fcs_float dx, fcs_float dy, fcs_float dz)
{
  return 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
}


fcs_float coulomb_field_3diff(const void *param, fcs_float dist, fcs_float dx, fcs_float dy, fcs_float dz)
{
  return -1.0 / (dx * dx + dy * dy + dz * dz);
}


void coulomb_field_potential_3diff(const void *param, fcs_float dist, fcs_float dx, fcs_float dy, fcs_float dz, fcs_float *f, fcs_float *p)
{
  *p = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
  *f = -(*p) * (*p);
}


FCS_NEAR_LOOP_F_P(coulomb_loop_f_p, coulomb_field, coulomb_potential);

FCS_NEAR_LOOP_FP(coulomb_loop_fp, coulomb_field_potential);

#define FCS_NEAR_LOOP2_POTENTIAL(_param_, _dist_)                  (1.0 / (_dist_))
#define FCS_NEAR_LOOP2_FIELD(_param_, _dist_)                      (-1.0 / ((_dist_) * (_dist_)))
FCS_NEAR_LOOP2_F_P(coulomb_loop2_f_p);

#define FCS_NEAR_LOOP2_FIELD_POTENTIAL(_param_, _dist_, _f_, _p_)  do { *(_p_) = (1.0 / (_dist_)); *(_f_) = -(*(_p_)) * (*(_p_)); } while (0)
FCS_NEAR_LOOP2_FP(coulomb_loop2_fp);

FCS_NEAR_LOOP_3DIFF_F_P(coulomb_loop_3diff_f_p, coulomb_field_3diff, coulomb_potential_3diff);

FCS_NEAR_LOOP_3DIFF_FP(coulomb_loop_3diff_fp, coulomb_field_potential_3diff);

#define FCS_NEAR_LOOP2_POTENTIAL_3DIFF(_param_, _dist_, _dx_, _dy_, _dz_)                  (1.0 / sqrt((_dx_) * (_dx_) + (_dy_) * (_dy_) + (_dz_) * (_dz_)))
#define FCS_NEAR_LOOP2_FIELD_3DIFF(_param_, _dist_, _dx_, _dy_, _dz_)                      (-1.0 / ((_dx_) * (_dx_) + (_dy_) * (_dy_) + (_dz_) * (_dz_)))
FCS_NEAR_LOOP2_3DIFF_F_P(coulomb_loop2_3diff_f_p);

#define FCS_NEAR_LOOP2_FIELD_POTENTIAL_3DIFF(_param_, _dist_, _dx_, _dy_, _dz_, _f_, _p_)  do { *(_p_) = (1.0 / sqrt((_dx_) * (_dx_) + (_dy_) * (_dy_) + (_dz_) * (_dz_))); *(_f_) = -(*(_p_)) * (*(_p_)); } while (0)
FCS_NEAR_LOOP2_3DIFF_FP(coulomb_loop2_3diff_fp);


void init_particles_homogen(fcs_int nlocal, fcs_float *xyz, fcs_float *box_base, fcs_float *box_a, fcs_float *box_b, fcs_float *box_c)
{
  fcs_int i;
  fcs_float r[3];


  for (i = 0; i < nlocal; ++i)
  {
    r[0] = (fcs_float) rand() / ((fcs_float) RAND_MAX + 1);
    r[1] = (fcs_float) rand() / ((fcs_float) RAND_MAX + 1);
    r[2] = (fcs_float) rand() / ((fcs_float) RAND_MAX + 1);

    xyz[3 * i + 0] = box_base[0] + box_a[0] * r[0] + box_b[0] * r[1] + box_c[0] * r[2];
    xyz[3 * i + 1] = box_base[1] + box_a[1] * r[0] + box_b[1] * r[1] + box_c[1] * r[2];
    xyz[3 * i + 2] = box_base[2] + box_a[2] * r[0] + box_b[2] * r[1] + box_c[2] * r[2];
  }
}


void print_particles(fcs_int nlocal, fcs_float *xyz)
{
  fcs_int i;


  for (i = 0; i < nlocal; ++i) printf("%" FCS_LMOD_INT "d  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f\n", i, xyz[3 * i + 0], xyz[3 * i + 1], xyz[3 * i + 2]);
}


void print_results(fcs_int nlocal, fcs_float *f, fcs_float *p)
{
  fcs_int i;


  for (i = 0; i < nlocal; ++i) printf("%" FCS_LMOD_INT "d  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f  %" FCS_LMOD_FLOAT "f\n", i, f[3 * i + 0], f[3 * i + 1], f[3 * i + 2], p[i]);
}


#define PRINT_PREFIX  /*"# "*/
/*#define PRINT_PARTICLES*/

#define TRICLINIC
#define PERIODIC
#define RUN_direct
#define RUN_NEAR


int main(int argc, char **argv)
{
  int comm_rank, comm_size;
  MPI_Comm comm = MPI_COMM_WORLD;

  fcs_int i;

  fcs_int nlocal, ntotal, nlocal_max;

  fcs_float *xyz, *q, *p, *f;
  fcs_float v[9];

  fcs_float box_base[] = { 0.0, 0.0, 0.0 };
#ifdef TRICLINIC
  fcs_float box_a[] = { 1.0, 2.0, 3.0 };
  fcs_float box_b[] = { 2.0, 3.0, 1.0 };
  fcs_float box_c[] = { 3.0, 1.0, 2.0 };
#else
  fcs_float box_a[] = { 1.0, 0.0, 0.0 };
  fcs_float box_b[] = { 0.0, 1.0, 0.0 };
  fcs_float box_c[] = { 0.0, 0.0, 1.0 };
#endif
  fcs_int periodicity[] = { 0, 0, 0 };

  fcs_float cutoff = 0.0;
  
  fcs_near_t near;
  
  char* method = "direct";

  double t;

  
  FCS fcs_handle;
  FCSResult fcs_result;
  
  fcs_float e_sum_local = 0.0;
  fcs_float e_sum = 0.0;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  
/*  if(argc == 1)
  {
    printf("No config file was specified!\n");
    MPI_Finalize();
    exit(-1);
  }*/

  srand(2501 * comm_rank);

  nlocal = 100;
  ntotal = comm_size * nlocal;
  nlocal_max = (nlocal * 110) / 100;

  cutoff = 0.1;

#ifdef PERIODIC
  periodicity[0] = periodicity[1] = periodicity[2] = 1;
#endif

  if (comm_rank == 0)
  {
    printf(PRINT_PREFIX "-------------------\n");
    printf(PRINT_PREFIX "Running Direct test\n");
    printf(PRINT_PREFIX "-------------------\n");
    printf(PRINT_PREFIX "  nprocs = %d\n", comm_size);
    printf(PRINT_PREFIX "  ntotal = %" FCS_LMOD_INT "d\n", ntotal);
    printf(PRINT_PREFIX "  box = [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]: [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f] x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f] x [%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f,%" FCS_LMOD_FLOAT "f]\n",
      box_base[0], box_base[1], box_base[2], box_a[0], box_a[1], box_a[2], box_b[0], box_b[1], box_b[2], box_c[0], box_c[1], box_c[2]);
    printf(PRINT_PREFIX "  cutoff = %" FCS_LMOD_FLOAT "f\n", cutoff);
    printf(PRINT_PREFIX "  periodicity = [%" FCS_LMOD_INT "d,%" FCS_LMOD_INT "d,%" FCS_LMOD_INT "d]\n", periodicity[0], periodicity[1], periodicity[2]);
  }

  xyz = malloc(nlocal_max * 3 * sizeof(fcs_float));
  q = malloc(nlocal_max * sizeof(fcs_float));
  f = malloc(nlocal_max * 3 * sizeof(fcs_float));
  p = malloc(nlocal_max * sizeof(fcs_float));

  init_particles_homogen(nlocal, xyz, box_base, box_a, box_b, box_c);

#if 1
  for (i = 0; i < nlocal; ++i) q[i] = 1e-2;
#else
  for (i = 0; i < nlocal; ++i) q[i] = (fcs_float) rand() / (fcs_float) RAND_MAX;
#endif

#ifdef PRINT_PARTICLES
  print_particles(nlocal, xyz);
#endif

  fcs_result = fcs_init(&fcs_handle, method, comm);
  ASSERT_FCS(fcs_result);

#ifdef RUN_direct
  fcs_direct_set_cutoff(fcs_handle, cutoff);
  ASSERT_FCS(fcs_result);

  fcs_result = fcs_set_common(fcs_handle, 1, box_a, box_b, box_c, box_base, periodicity, ntotal);
  ASSERT_FCS(fcs_result);

  fcs_result = fcs_require_virial(fcs_handle, 1);
  ASSERT_FCS(fcs_result);

  fcs_result = fcs_set_periodicity(fcs_handle, periodicity);
  ASSERT_FCS(fcs_result);

  for (i = 0; i < nlocal; ++i) p[i] = f[i * 3 + 0] = f[i * 3 + 1] = f[i * 3 + 2] = 0;

  fcs_result = fcs_tune(fcs_handle, nlocal, nlocal_max, xyz, q);
  ASSERT_FCS(fcs_result);

  MPI_Barrier(comm);
  t = MPI_Wtime();
  fcs_result = fcs_run(fcs_handle, nlocal, nlocal_max, xyz, q, f, p);
  ASSERT_FCS(fcs_result);
  MPI_Barrier(comm);
  t = MPI_Wtime() - t;

  e_sum_local = 0.0;
  for(i = 0; i < nlocal; ++i) e_sum_local += p[i] * q[i];

  MPI_Reduce(&e_sum_local, &e_sum, 1, FCS_MPI_FLOAT, MPI_SUM, 0, comm);

  fcs_get_virial(fcs_handle, v);

/*  print_results(nlocal, f, p);*/

  if (comm_rank == 0)
  {
    printf(PRINT_PREFIX "direct: %f\n", t);
/*    printf("  approx. Madelung's constant: %e\n", e_sum/ntotal * 2.0 * M_PI);*/
    printf(PRINT_PREFIX "  total energy: %.16" FCS_LMOD_FLOAT "e\n", e_sum);
/*    printf("  virial: %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e %" FCS_LMOD_FLOAT "e\n",
      v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);*/
  }
#endif

#ifdef RUN_NEAR
  for (i = 0; i < nlocal; ++i) p[i] = f[i * 3 + 0] = f[i * 3 + 1] = f[i * 3 + 2] = 0;

  fcs_near_create(&near);

#define NEAR_TYPE  0
/*#define NEAR_3DIFF*/

#ifdef NEAR_3DIFF
#if NEAR_TYPE == 0
  fcs_near_set_field_3diff(&near, coulomb_field_3diff);
  fcs_near_set_potential_3diff(&near, coulomb_potential_3diff);
#elif NEAR_TYPE == 1
  fcs_near_set_field_potential_3diff(&near, coulomb_field_potential_3diff);
#elif NEAR_TYPE == 2
  fcs_near_set_loop(&near, coulomb_loop_3diff_f_p);
#elif NEAR_TYPE == 3
  fcs_near_set_loop(&near, coulomb_loop_3diff_fp);
#elif NEAR_TYPE == 4
  fcs_near_set_loop(&near, coulomb_loop2_3diff_f_p);
#elif NEAR_TYPE == 5
  fcs_near_set_loop(&near, coulomb_loop2_3diff_fp);
#endif
#else
#if NEAR_TYPE == 0
  fcs_near_set_field(&near, coulomb_field);
  fcs_near_set_potential(&near, coulomb_potential);
#elif NEAR_TYPE == 1
  fcs_near_set_field_potential(&near, coulomb_field_potential);
#elif NEAR_TYPE == 2
  fcs_near_set_loop(&near, coulomb_loop_f_p);
#elif NEAR_TYPE == 3
  fcs_near_set_loop(&near, coulomb_loop_fp);
#elif NEAR_TYPE == 4
  fcs_near_set_loop(&near, coulomb_loop2_f_p);
#elif NEAR_TYPE == 5
  fcs_near_set_loop(&near, coulomb_loop2_fp);
#endif
#endif

  fcs_near_set_system(&near, box_base, box_a, box_b, box_c, periodicity);

  fcs_near_set_particles(&near, nlocal, nlocal_max, xyz, q, NULL, f, p);

  MPI_Barrier(comm);
  t = MPI_Wtime();
  fcs_near_field_solver(&near, cutoff, NULL, comm);
  MPI_Barrier(comm);
  t = MPI_Wtime() - t;

  fcs_near_destroy(&near);

  e_sum_local = 0.0;
  for(i = 0; i < nlocal; ++i) e_sum_local += p[i] * q[i];
  MPI_Reduce(&e_sum_local, &e_sum, 1, FCS_MPI_FLOAT, MPI_SUM, 0, comm);

/*  print_results(nlocal, f, p);*/

  if (comm_rank == 0)
  {
    printf(PRINT_PREFIX "NEAR: %f\n", t);
    printf(PRINT_PREFIX "  total energy: %.16" FCS_LMOD_FLOAT "e\n", e_sum);
  }
#endif

  fcs_destroy(fcs_handle);

  free(xyz);
  free(q);
  free(f);
  free(p);

  if (comm_rank == 0) printf(PRINT_PREFIX "Direct done.\n");

  MPI_Finalize();

  return 0;
}
