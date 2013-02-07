/*
 * Copyright (c) 2011-2013 Michael Pippig
 *
 * This file is part of PNFFT.
 *
 * PNFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PNFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PNFFT.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "pnfft.h"
#include "ipnfft.h"

static void timer_reset(
    double *timer);
static int file_exists(
    const char *name); 
static FILE* open_or_create_file_to_append(
    MPI_Comm comm, const char *name);
static void write_info_header(
    MPI_Comm comm, FILE *file);
static int get_matlab_index(
    MPI_Comm comm);
static void write_run_specific_infos(
    MPI_Comm comm, FILE *file, PNX(plan) ths);
static void fprint_average_timer_internal(
    MPI_Comm comm, FILE *file, const char *prefix, double *ths,
    unsigned flags);
static void fprint_average_timer(
    MPI_Comm comm, FILE *file, const PNX(plan) ths, unsigned flags);


double* PNX(get_timer_trafo)(
    PNX(plan) ths
    )
{
  return PNX(timer_copy)(ths->timer_trafo);
}

double* PNX(get_timer_adj)(
    PNX(plan) ths
    )
{
  return PNX(timer_copy)(ths->timer_adj);
}

void PNX(timer_average)(
     double *timer
     )
{
  if(timer[PNFFT_TIMER_ITER] < 1.0)
    return;
  
  for(int t=1; t<PNFFT_TIMER_LENGTH; t++)
    timer[t] /= timer[0];
}

double* PNX(timer_copy)(
    const double *orig
    )
{
  double *copy = PNX(mktimer)();
 
  for(int t=0; t<PNFFT_TIMER_LENGTH; t++)
    copy[t] = orig[t];
  return copy;
}

double* PNX(timer_reduce_max)(
    MPI_Comm comm, double *timer
    )
{
  double *timer_max = PNX(mktimer)();
  MPI_Reduce(timer, timer_max, PNFFT_TIMER_LENGTH, MPI_DOUBLE, MPI_MAX, 0, comm);
  return timer_max;
}

double* PNX(timer_add)(
    const double *sum1, const double *sum2
    )
{
  double *sum = PNX(mktimer)();
  
  for(int t=0; t<PNFFT_TIMER_LENGTH; t++)
    sum[t] = sum1[t] + sum2[t];

  return sum;
}

void PNX(timer_free)(
    double* ths
    )
{
  PNX(rmtimer)(ths); 
}

void PNX(reset_timer)(
    PNX(plan) ths
    )
{
  timer_reset(ths->timer_trafo);
  timer_reset(ths->timer_adj);
}

static void timer_reset(
    double *timer
    )
{
  for(int t=0; t<PNFFT_TIMER_LENGTH; t++)
    timer[t] = 0;
}

double* PNX(mktimer)(
    void
    )
{
  double* timer = (double*) malloc(sizeof(double) * PNFFT_TIMER_LENGTH);
  timer_reset(timer);
  return timer;
}

void PNX(rmtimer)(
    double* timer
    )
{
  if(timer != NULL)
    free(timer);
}


void PNX(write_average_timer)(
    const PNX(plan) ths, const char *name, MPI_Comm comm
    )
{
  int newfile;
  FILE *f;
  
  newfile = !file_exists(name);
  f = open_or_create_file_to_append(comm, name);
  if( newfile )
    write_info_header(comm, f);
  
  write_run_specific_infos(comm, f, ths);
  fprint_average_timer(comm, f, ths, PNFFT_PRINT_TIMER_BASIC);
  
  fclose(f);
}

void PNX(write_average_timer_adv)(
    const PNX(plan) ths, const char *name, MPI_Comm comm
    )
{
  FILE *f;
  PNX(write_average_timer)(ths, name, comm);
  f = open_or_create_file_to_append(comm, name);
  fprint_average_timer(comm, f, ths, PNFFT_PRINT_TIMER_ADV);
  fclose(f);
  PX(write_average_timer_adv)(ths->pfft_forw, name, comm);
  PX(write_average_timer_adv)(ths->pfft_back, name, comm);
  PX(write_average_gctimer_adv)(ths->gcplan, name, comm);
}

void PNX(print_average_timer)(
    const PNX(plan) ths, MPI_Comm comm
    )
{
  write_info_header(comm, stdout);
  write_run_specific_infos(comm, stdout, ths);
  fprint_average_timer(comm, stdout, ths, PNFFT_PRINT_TIMER_BASIC);
}

void PNX(print_average_timer_adv)(
    const PNX(plan) ths, MPI_Comm comm
    )
{
  PNX(print_average_timer)(ths, comm);
  fprint_average_timer(comm, stdout, ths, PNFFT_PRINT_TIMER_ADV);
  PX(print_average_timer_adv)(ths->pfft_forw, comm);
  PX(print_average_timer_adv)(ths->pfft_back, comm);
  PX(print_average_gctimer_adv)(ths->gcplan, comm);
}

static int file_exists(
    const char *name
    )
{
  FILE *f;
  
  f=fopen(name, "r");
  if(f != NULL)
    fclose(f);
  
  return (f != NULL);
}

static FILE* open_or_create_file_to_append(
    MPI_Comm comm, const char *name
    )
{
  FILE *f=fopen(name, "a+");
  if(f==NULL){
    PX(fprintf)(comm, stderr, "Error: Cannot open file %s.\n", name);
    exit(1);
  }
  return f;
}

static void write_info_header(
    MPI_Comm comm, FILE *file
    )
{
  PX(fprintf)(comm, file, "%% N  - NFFT size\n");
  PX(fprintf)(comm, file, "%% n  - FFT size\n");
  PX(fprintf)(comm, file, "%% np - process grid\n");
  PX(fprintf)(comm, file, "%% procs - number of processes\n");
  PX(fprintf)(comm, file, "%% pnfft - PNFFT runtime\n");
  PX(fprintf)(comm, file, "%% pfft  - PFFT runtime\n");
  PX(fprintf)(comm, file, "%% index(i) = log(procs(i)) + 1\n");
}

static int get_matlab_index(
    MPI_Comm comm
    )
{
  int size;
  MPI_Comm_size(comm, &size);
  return lrint(round(log((R) size)/log(2.0)))+1;
}

static void write_run_specific_infos(
    MPI_Comm comm, FILE *file, PNX(plan) ths
    )
{
  
  int size;
  int idx = get_matlab_index(comm);
  MPI_Comm_size(comm, &size);
  
  if(ths->pnfft_flags & PNFFT_WINDOW_GAUSSIAN)
    PX(fprintf)(comm, file, "\n%% pnfft_flags == PNFFT_WINDOW_GAUSSIAN");
  else if(ths->pnfft_flags & PNFFT_WINDOW_BSPLINE)
    PX(fprintf)(comm, file, "\n%% pnfft_flags == PNFFT_WINDOW_BSPLINE");
  else if(ths->pnfft_flags & PNFFT_WINDOW_SINC_POWER)
    PX(fprintf)(comm, file, "\n%% pnfft_flags == PNFFT_WINDOW_SINC_POWER");
  else if(ths->pnfft_flags & PNFFT_WINDOW_BESSEL_I0)
    PX(fprintf)(comm, file, "\n%% pnfft_flags == PNFFT_WINDOW_BESSEL_I0");
  else
    PX(fprintf)(comm, file, "\n%% pnfft_flags == PNFFT_WINDOW_KAISER_BESSEL");

  if(ths->pnfft_flags & PNFFT_PRE_PHI_HAT)
    PX(fprintf)(comm, file, " | PNFFT_PRE_PHI_HAT");
  if(ths->pnfft_flags & PNFFT_FG_PSI)
    PX(fprintf)(comm, file, " | PNFFT_FG_PSI");
  if(ths->pnfft_flags & PNFFT_PRE_LIN_PSI)
    PX(fprintf)(comm, file, " | PNFFT_PRE_LIN_PSI");
  if(ths->pnfft_flags & PNFFT_PRE_QUAD_PSI)
    PX(fprintf)(comm, file, " | PNFFT_PRE_QUAD_PSI");
  if(ths->pnfft_flags & PNFFT_PRE_KUB_PSI)
    PX(fprintf)(comm, file, " | PNFFT_PRE_KUB_PSI");
  if(ths->pnfft_flags & PNFFT_PRE_FG_PSI)
    PX(fprintf)(comm, file, " | PNFFT_PRE_FG_PSI");
  if(ths->pnfft_flags & PNFFT_PRE_PSI)
    PX(fprintf)(comm, file, " | PNFFT_PRE_PSI");
  if(ths->pnfft_flags & PNFFT_PRE_FULL_PSI)
    PX(fprintf)(comm, file, " | PNFFT_PRE_FULL_PSI");

  if(ths->pnfft_flags & PNFFT_FFT_IN_PLACE)
    PX(fprintf)(comm, file, " | PNFFT_FFT_IN_PLACE");
  else
    PX(fprintf)(comm, file, " | PNFFT_FFT_OUT_OF_PLACE");
  if(ths->pnfft_flags & PNFFT_SORT_NODES)
    PX(fprintf)(comm, file, " | PNFFT_SORT_NODES");
  if(ths->pnfft_flags & PNFFT_INTERLACED)
    PX(fprintf)(comm, file, " | PNFFT_INTERLACED");
  if(ths->pnfft_flags & PNFFT_SHIFTED_IN)
    PX(fprintf)(comm, file, " | PNFFT_SHIFTED_IN");
  if(ths->pnfft_flags & PNFFT_SHIFTED_OUT)
    PX(fprintf)(comm, file, " | PNFFT_SHIFTED_OUT");

  if(ths->pnfft_flags & PNFFT_GRAD_IK)
    PX(fprintf)(comm, file, " | PNFFT_GRAD_IK");
  else if(ths->pnfft_flags & PNFFT_GRAD_NONE)
    PX(fprintf)(comm, file, " | PNFFT_GRAD_NONE");
  else
    PX(fprintf)(comm, file, " | PNFFT_GRAD_AD");

  if(ths->pnfft_flags & PNFFT_REAL_F)
    PX(fprintf)(comm, file, " | PNFFT_REAL_F");

  if(ths->pnfft_flags & PNFFT_MALLOC_F_HAT)
    PX(fprintf)(comm, file, " | PNFFT_MALLOC_F_HAT");
  if(ths->pnfft_flags & PNFFT_MALLOC_X)
    PX(fprintf)(comm, file, " | PNFFT_MALLOC_X");
  if(ths->pnfft_flags & PNFFT_MALLOC_F)
    PX(fprintf)(comm, file, " | PNFFT_MALLOC_F");
  if(ths->pnfft_flags & PNFFT_MALLOC_GRAD_F)
    PX(fprintf)(comm, file, " | PNFFT_MALLOC_GRAD_F");
  PX(fprintf)(comm, file, "\n");

  PX(fprintf)(comm, file, "\nindex(%d) = %d;  ", idx, idx);
  PX(fprintf)(comm, file, "procs(%d) = %d;  ", idx, size);
  PX(fprintf)(comm, file, "np_pnfft(%d, 1:3) = [%d %d %d];  ", idx, ths->np[0], ths->np[1], ths->np[2]);
  PX(fprintf)(comm, file, "N_pnfft(%d, 1:%d) = [", idx, ths->d);
  for(int t=0; t<ths->d; t++)
    PX(fprintf)(comm, file, "%td ", ths->N[t]);
  PX(fprintf)(comm, file, "];  ");
  PX(fprintf)(comm, file, "n_pnfft(%d, 1:%d) = [", idx, ths->d);
  for(int t=0; t<ths->d; t++)
    PX(fprintf)(comm, file, "%td ", ths->n[t]);
  PX(fprintf)(comm, file, "];  ");
  PX(fprintf)(comm, file, "m_pnfft(%d) = %d;\n", idx, ths->m);
}

static void fprint_average_timer_internal(
    MPI_Comm comm, FILE *file, const char *prefix, double *timer,
    unsigned flags
    )
{
  int idx = get_matlab_index(comm);
  double *mt;
  mt = PNX(timer_reduce_max)(comm, timer);
  PNX(timer_average)(mt);

  if(flags & PNFFT_PRINT_TIMER_BASIC){
    PX(fprintf)(comm, file, "%s_iter(%d)    = %d;  ", prefix, idx, (int) timer[PNFFT_TIMER_ITER]);
    PX(fprintf)(comm, file, "%s(%d)   = %.3e;\n", prefix, idx, mt[PNFFT_TIMER_WHOLE]);
  } else if(flags & PNFFT_PRINT_TIMER_ADV){
    PX(fprintf)(comm, file, "%s_matrix_D(%d)   = %.3e;  ", prefix, idx, mt[PNFFT_TIMER_MATRIX_D]);
    PX(fprintf)(comm, file, "%s_matrix_F(%d)   = %.3e;\n", prefix, idx, mt[PNFFT_TIMER_MATRIX_F]);
    PX(fprintf)(comm, file, "%s_matrix_B(%d)   = %.3e;  ", prefix, idx, mt[PNFFT_TIMER_MATRIX_B]);
    PX(fprintf)(comm, file, "%s_gcells(%d)     = %.3e;\n", prefix, idx, mt[PNFFT_TIMER_GCELLS]);
    PX(fprintf)(comm, file, "%s_sort_nodes(%d) = %.3e;  ", prefix, idx, mt[PNFFT_TIMER_SORT_NODES]);
    PX(fprintf)(comm, file, "%s_loop_B(%d)     = %.3e;\n", prefix, idx, mt[PNFFT_TIMER_LOOP_B]);
    PX(fprintf)(comm, file, "%s_shift_in(%d)   = %.3e;  ", prefix, idx, mt[PNFFT_TIMER_SHIFT_INPUT]);
    PX(fprintf)(comm, file, "%s_shift_out(%d)  = %.3e;\n", prefix, idx, mt[PNFFT_TIMER_SHIFT_OUTPUT]);
  }
}

static void fprint_average_timer(
    MPI_Comm comm, FILE *file, const PNX(plan) ths, unsigned flags
    )
{
  fprint_average_timer_internal(comm, file, "pnfft_trf", ths->timer_trafo, flags);
  fprint_average_timer_internal(comm, file, "pnfft_adj", ths->timer_adj, flags);
}
  


