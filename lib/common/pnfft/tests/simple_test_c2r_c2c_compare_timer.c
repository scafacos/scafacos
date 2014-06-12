#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>

#include <fftw3.h>
#include <pfft.h>
#include <pnfft.h>

#define LEN 1000

void message(char *string, MPI_Comm comm);

static void init_equispaced_x(const ptrdiff_t *N, const double *lo, const double *up, double *x);
static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, double complex *data);

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *N, int *np,
    ptrdiff_t *local_M,
    int *repetitions);


int main(int argc, char **argv)
{
  int np[3];
  ptrdiff_t N[3];
  ptrdiff_t local_M;
  double err;
  MPI_Comm comm_cart_3d;
  int repetitions;
  char str[LEN];

  ptrdiff_t local_N_c2c[3], local_N_start_c2c[3];
  double lower_border_c2c[3], upper_border_c2c[3];
  pnfft_plan plan_c2c;
  pnfft_complex *f_hat_c2c, *f_c2c;
  double *x_c2c;

  ptrdiff_t local_N_c2r[3], local_N_start_c2r[3];
  double lower_border_c2r[3], upper_border_c2r[3];
  pnfft_plan plan_c2r;
  pnfft_complex *f_hat_c2r;
  double *x_c2r, *f_c2r;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pnfft_init();

  np[0] = 2; np[1] = 2; np[2] = 8;
  N[0] = 64; N[1] = 64; N[2] = 128;
  local_M = N[0]*N[1]*N[2]/(np[0]*np[1]*np[2]);

  repetitions=5;

  /* set parameters by command line */
  init_parameters(argc, argv, N, np, &local_M, &repetitions);

   /* Print infos */
  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of parallel NFFT\n");
  pfft_printf(MPI_COMM_WORLD, "* for  N[0] x N[1] x N[2] = %td x %td x %td Fourier coefficients)\n", N[0], N[1], N[2]);
  pfft_printf(MPI_COMM_WORLD, "* at   local_M = %td nodes per process\n", local_M);
  pfft_printf(MPI_COMM_WORLD, "* on   np[0] x np[1] x np[2] = %td x %td x %td processes\n", np[0], np[1], np[2]);
  pfft_printf(MPI_COMM_WORLD, "* with %d repetitions\n", repetitions);
  pfft_printf(MPI_COMM_WORLD, "*******************************************************************************************************\n\n");

  /* create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pnfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh of size %d x %d x %d does not fit to number of allocated processes.\n", np[0], np[1], np[2]);
    pfft_fprintf(MPI_COMM_WORLD, stderr, "       Please allocate %d processes (mpiexec -np %d ...) or change the procmesh (with -pnfft_np * * *).\n", np[0]*np[1]*np[2], np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  pnfft_local_size_3d(N, comm_cart_3d, PNFFT_TRANSPOSED_NONE,
      local_N_c2c, local_N_start_c2c, lower_border_c2c, upper_border_c2c);
  pnfft_local_size_3d_c2r(N, comm_cart_3d, PNFFT_TRANSPOSED_NONE,
      local_N_c2r, local_N_start_c2r, lower_border_c2r, upper_border_c2r);

  /* Plan parallel NFFT */
  plan_c2c = pnfft_init_adv(3, N, local_M,
      PNFFT_TRANSPOSED_NONE| PNFFT_WINDOW_SINC_POWER| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
      comm_cart_3d);
  plan_c2r = pnfft_init_adv_c2r(3, N, local_M,
      PNFFT_TRANSPOSED_NONE| PNFFT_WINDOW_SINC_POWER| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
      comm_cart_3d);


  f_hat_c2c = pnfft_get_f_hat(plan_c2c);
  f_c2c     = pnfft_get_f(plan_c2c);
  x_c2c     = pnfft_get_x(plan_c2c);
  f_hat_c2r = pnfft_get_f_hat(plan_c2r);
  f_c2r     = pnfft_get_f(plan_c2r);
  x_c2r     = pnfft_get_x(plan_c2r);

  /* Initialize Fourier coefficients with random numbers */
  init_input(N, local_N_c2c, local_N_start_c2c, f_hat_c2c);
  init_input(N, local_N_c2r, local_N_start_c2r, f_hat_c2r);

  /* Initialize nodes with random numbers */
  pnfft_init_x_3d(lower_border_c2c, upper_border_c2c, local_M, x_c2c);
//  init_equispaced_x(N, lower_border_c2c, upper_border_c2c, x_c2c);
  for (int k=0; k<local_M*3; k++)
    x_c2r[k] = x_c2c[k];

  /* execute parallel NFFT */
  for (int k=0; k<repetitions; k++) {
    MPI_Barrier(MPI_COMM_WORLD);
    pnfft_trafo(plan_c2c);
    MPI_Barrier(MPI_COMM_WORLD);
    pnfft_adj(plan_c2c);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  snprintf(str, LEN, "./c2c_timings-%td_%td_%td-%d_%d_%d-%d.txt", N[0], N[1], N[2], np[0], np[1], np[2], repetitions);
  pnfft_write_average_timer_adv(plan_c2c, str, comm_cart_3d);


  /* execute parallel NFFT */
  for (int k=0; k<repetitions; k++) {
    MPI_Barrier(MPI_COMM_WORLD);
    pnfft_trafo(plan_c2r);
    MPI_Barrier(MPI_COMM_WORLD);
    pnfft_adj(plan_c2r);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  snprintf(str, LEN, "./c2r_timings-%td_%td_%td-%d_%d_%d-%d.txt", N[0], N[1], N[2], np[0], np[1], np[2], repetitions);
  pnfft_write_average_timer_adv(plan_c2r, str, comm_cart_3d);



  /* free mem and finalize */
  pnfft_finalize(plan_c2c, PNFFT_FREE_X | PNFFT_FREE_F_HAT | PNFFT_FREE_F);
  pnfft_finalize(plan_c2r, PNFFT_FREE_X | PNFFT_FREE_F_HAT | PNFFT_FREE_F);
  MPI_Comm_free(&comm_cart_3d);

  /* Finalize MPI */
  MPI_Finalize();
  return 0;
}

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *N, int *np, 
    ptrdiff_t *local_M,
    int *repetitions
    )
{
  pfft_get_args(argc, argv, "-pnfft_N", 3, PFFT_PTRDIFF_T, N);
  pfft_get_args(argc, argv, "-pnfft_np", 3, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pnfft_local_M", 1, PFFT_PTRDIFF_T, local_M);
  pfft_get_args(argc, argv, "-pnfft_repetitions", 1, PFFT_INT, repetitions);
}

#define DATA_INIT(i) (( (double)1000 ) / ( (double)( (i) == 0 ? 1 : i) ))

static pfft_complex semirandom(const ptrdiff_t *N, ptrdiff_t k0, ptrdiff_t k1, ptrdiff_t k2) {
  if ((k0 == -N[0]/2) || (k1 == -N[1]/2) || (k2 == -N[2]/2))
    return 0 + 0*I;
  if (k0 == (N[0]+1)/2)
    return semirandom(N, -k0, k1, k2);
  if (k1 == (N[1]+1)/2)
    return semirandom(N, k0, -k1, k2);
  if (k2 == (N[2]+1)/2)
    return semirandom(N, k0, k1, -k2);

  ptrdiff_t l0 = (k0) ? k0+N[0] : 1;
  ptrdiff_t l1 = (k1) ? k1-2*N[1] : 1;
  ptrdiff_t l2 = (k2) ? k2+N[2] : 1;

  double re = DATA_INIT(l0+l1*l2) + 3*DATA_INIT(l0*l1+l2) + DATA_INIT(l2*l0+l1);
  double im = 3*DATA_INIT(l0+l1*l1) + DATA_INIT(l1+l2*l2) - DATA_INIT(l2+l0*l0);

  return re + im*I;
}

static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, pfft_complex *data)
{
  int m = 0;
  for(ptrdiff_t k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
    for(ptrdiff_t k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
      for(ptrdiff_t k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++, m++)
        data[m] = semirandom(N, k0, k1, k2) + conj(semirandom(N, -k0, -k1, -k2));
}

static void init_equispaced_x(
    const ptrdiff_t *N, const double *lo, const double *up,
    double *x
    )
{
  /* enter your code here and call make when you are finished */

  ptrdiff_t local_N[3], local_N_start[3], m=0;
  for(int t=0; t<3; t++){
    local_N[t] = ((up[t]-lo[t]) * N[t]);
    local_N_start[t] = lo[t] * N[t];
  }

  m=0;
  for(ptrdiff_t k0=local_N_start[0]; k0<local_N_start[0] + local_N[0]; k0++)
    for(ptrdiff_t k1=local_N_start[1]; k1<local_N_start[1] + local_N[1]; k1++)
      for(ptrdiff_t k2=local_N_start[2]; k2<local_N_start[2] + local_N[2]; k2++, m++){
        x[3*m+0] = (double) k0 / N[0];
        x[3*m+1] = (double) k1 / N[1];
        x[3*m+2] = (double) k2 / N[2];
      }
}

