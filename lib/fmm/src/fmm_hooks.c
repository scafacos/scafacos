
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "fmm_hooks.h"


fcs_fmm_hooks_t *fcs_fmm_hooks_create()
{
  printf("fcs_fmm_hooks_create:\n");

  fcs_fmm_hooks_t *hooks = malloc(sizeof(fcs_fmm_hooks_t));

  hooks->demo_value = 0;

  hooks->t_near = hooks->t_far = 0;

  printf("fcs_fmm_hooks_create: return: hooks: %p\n", hooks);

  return hooks;
}


void fcs_fmm_hooks_destroy(fcs_fmm_hooks_t *hooks)
{
  printf("fcs_fmm_hooks_destroy: hooks: %p\n", hooks);

  free(hooks);

  printf("fcs_fmm_hooks_destroy: return\n");
}


void fcs_fmm_hooks_near_start(fcs_fmm_hooks_t *hooks)
{
  printf("fcs_fmm_hooks_near_start: hooks: %p\n", hooks);

  ++hooks->demo_value;

  hooks->t_near_d = MPI_Wtime();

  printf("fcs_fmm_hooks_near_start: return\n");
}


void fcs_fmm_hooks_near_stop(fcs_fmm_hooks_t *hooks)
{
  printf("fcs_fmm_hooks_near_stop: hooks: %p\n", hooks);

  ++hooks->demo_value;

  double t = MPI_Wtime();
  hooks->t_near += t - hooks->t_near_d;
  hooks->t_near_d = t;

  printf("fcs_fmm_hooks_near_stop: return\n");
}


void fcs_fmm_hooks_far_start(fcs_fmm_hooks_t *hooks)
{
  printf("fcs_fmm_hooks_far_start: hooks: %p\n", hooks);

  ++hooks->demo_value;

  hooks->t_far_d = MPI_Wtime();

  printf("fcs_fmm_hooks_far_start: return\n");
}


void fcs_fmm_hooks_far_stop(fcs_fmm_hooks_t *hooks)
{
  printf("fcs_fmm_hooks_far_stop: hooks: %p\n", hooks);

  ++hooks->demo_value;

  double t = MPI_Wtime();
  hooks->t_far += t - hooks->t_far_d;
  hooks->t_far_d = t;

  printf("fcs_fmm_hooks_far_stop: return\n");
}


void fcs_fmm_hooks_print(fcs_fmm_hooks_t *hooks)
{
  printf("fcs_fmm_hooks_print: hooks: %p\n", hooks);

  printf("fcs_fmm_hooks_print: demo_value: %d\n", hooks->demo_value);

  printf("fcs_fmm_hooks_print: t_near: %f\n", hooks->t_near);
  printf("fcs_fmm_hooks_print: t_far: %f\n", hooks->t_far);

  printf("fcs_fmm_hooks_print: return\n");
}
