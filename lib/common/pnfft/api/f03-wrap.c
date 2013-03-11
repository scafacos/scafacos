/* Generated automatically.  DO NOT EDIT! */

#include "pnfft.h"
#include "ipnfft.h"

PNFFT_EXTERN int PNX(create_procmesh_2d_f03)(MPI_Fint f_comm, int np0, int np1, MPI_Fint * f_comm_cart_2d);
PNFFT_EXTERN int PNX(create_procmesh_f03)(int rnk, MPI_Fint f_comm, const int * np, MPI_Fint * f_comm_cart);
PNFFT_EXTERN void PNX(local_size_3d_f03)(const INT * N, MPI_Fint f_comm_cart, INT * local_N, INT * local_N_start, R * lower_border, R * upper_border);
PNFFT_EXTERN void PNX(local_size_adv_f03)(int d, const INT * N, MPI_Fint f_comm_cart, INT * local_N, INT * local_N_start, R * lower_border, R * upper_border);
PNFFT_EXTERN void PNX(local_size_guru_f03)(int d, const INT * N, const INT * n, const R * x_max, int m, MPI_Fint f_comm_cart, INT * local_N, INT * local_N_start, R * lower_border, R * upper_border);
PNFFT_EXTERN PNX(plan) PNX(init_3d_f03)(const INT * N, INT local_M, MPI_Fint f_comm_cart);
PNFFT_EXTERN PNX(plan) PNX(init_adv_f03)(int d, const INT * N, INT local_M, unsigned pnfft_flags, unsigned fftw_flags, MPI_Fint f_comm_cart);
PNFFT_EXTERN PNX(plan) PNX(init_guru_f03)(int d, const INT * N, const INT * n, const R * x_max, INT local_M, int m, unsigned pnfft_flags, unsigned fftw_flags, MPI_Fint f_comm_cart);
PNFFT_EXTERN void PNX(vpr_complex_f03)(C * data, INT N, const char * name, MPI_Fint f_comm);
PNFFT_EXTERN void PNX(vpr_real_f03)(R * data, INT N, const char * name, MPI_Fint f_comm);
PNFFT_EXTERN void PNX(apr_complex_3d_f03)(C * data, INT * local_N, INT * local_N_start, const char * name, MPI_Fint f_comm);
PNFFT_EXTERN double * PNX(timer_reduce_max_f03)(MPI_Fint f_comm, double * timer);
PNFFT_EXTERN void PNX(print_average_timer_f03)(const PNX(plan) ths, MPI_Fint f_comm);
PNFFT_EXTERN void PNX(print_average_timer_adv_f03)(const PNX(plan) ths, MPI_Fint f_comm);
PNFFT_EXTERN void PNX(write_average_timer_f03)(const PNX(plan) ths, const char * name, MPI_Fint f_comm);
PNFFT_EXTERN void PNX(write_average_timer_adv_f03)(const PNX(plan) ths, const char * name, MPI_Fint f_comm);

int PNX(create_procmesh_2d_f03)(MPI_Fint f_comm, int np0, int np1, MPI_Fint * f_comm_cart_2d)
{
  MPI_Comm comm, comm_cart_2d;

  comm = MPI_Comm_f2c(f_comm);
  int ret = PNX(create_procmesh_2d)(comm, np0, np1, &comm_cart_2d);
  *f_comm_cart_2d = MPI_Comm_c2f(comm_cart_2d);
  return ret;
}

int PNX(create_procmesh_f03)(int rnk, MPI_Fint f_comm, const int * np, MPI_Fint * f_comm_cart)
{
  MPI_Comm comm, comm_cart;

  comm = MPI_Comm_f2c(f_comm);
  int ret = PNX(create_procmesh)(rnk, comm, np, &comm_cart);
  *f_comm_cart = MPI_Comm_c2f(comm_cart);
  return ret;
}

void PNX(local_size_3d_f03)(const INT * N, MPI_Fint f_comm_cart, INT * local_N, INT * local_N_start, R * lower_border, R * upper_border)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PNX(local_size_3d)(N, comm_cart, local_N, local_N_start, lower_border, upper_border);
}

void PNX(local_size_adv_f03)(int d, const INT * N, MPI_Fint f_comm_cart, INT * local_N, INT * local_N_start, R * lower_border, R * upper_border)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PNX(local_size_adv)(d, N, comm_cart, local_N, local_N_start, lower_border, upper_border);
}

void PNX(local_size_guru_f03)(int d, const INT * N, const INT * n, const R * x_max, int m, MPI_Fint f_comm_cart, INT * local_N, INT * local_N_start, R * lower_border, R * upper_border)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PNX(local_size_guru)(d, N, n, x_max, m, comm_cart, local_N, local_N_start, lower_border, upper_border);
}

PNX(plan) PNX(init_3d_f03)(const INT * N, INT local_M, MPI_Fint f_comm_cart)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PNX(plan) ret = PNX(init_3d)(N, local_M, comm_cart);
  return ret;
}

PNX(plan) PNX(init_adv_f03)(int d, const INT * N, INT local_M, unsigned pnfft_flags, unsigned fftw_flags, MPI_Fint f_comm_cart)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PNX(plan) ret = PNX(init_adv)(d, N, local_M, pnfft_flags, fftw_flags, comm_cart);
  return ret;
}

PNX(plan) PNX(init_guru_f03)(int d, const INT * N, const INT * n, const R * x_max, INT local_M, int m, unsigned pnfft_flags, unsigned fftw_flags, MPI_Fint f_comm_cart)
{
  MPI_Comm comm_cart;

  comm_cart = MPI_Comm_f2c(f_comm_cart);
  PNX(plan) ret = PNX(init_guru)(d, N, n, x_max, local_M, m, pnfft_flags, fftw_flags, comm_cart);
  return ret;
}

void PNX(vpr_complex_f03)(C * data, INT N, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PNX(vpr_complex)(data, N, name, comm);
}

void PNX(vpr_real_f03)(R * data, INT N, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PNX(vpr_real)(data, N, name, comm);
}

void PNX(apr_complex_3d_f03)(C * data, INT * local_N, INT * local_N_start, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PNX(apr_complex_3d)(data, local_N, local_N_start, name, comm);
}

double * PNX(timer_reduce_max_f03)(MPI_Fint f_comm, double * timer)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  double * ret = PNX(timer_reduce_max)(comm, timer);
  return ret;
}

void PNX(print_average_timer_f03)(const PNX(plan) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PNX(print_average_timer)(ths, comm);
}

void PNX(print_average_timer_adv_f03)(const PNX(plan) ths, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PNX(print_average_timer_adv)(ths, comm);
}

void PNX(write_average_timer_f03)(const PNX(plan) ths, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PNX(write_average_timer)(ths, name, comm);
}

void PNX(write_average_timer_adv_f03)(const PNX(plan) ths, const char * name, MPI_Fint f_comm)
{
  MPI_Comm comm;

  comm = MPI_Comm_f2c(f_comm);
  PNX(write_average_timer_adv)(ths, name, comm);
}
