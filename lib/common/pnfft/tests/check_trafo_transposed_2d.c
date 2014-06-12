#include <stdlib.h>
#include <complex.h>
#include <pnfft.h>

static void pnfft_perform_guru(
    const ptrdiff_t *N, const ptrdiff_t *n, ptrdiff_t local_M,
    int m, const double *x_max, unsigned window_flag,
    const int *np, MPI_Comm comm,
    pnfft_complex **f, double *f_hat_sum);

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *N, ptrdiff_t *n, ptrdiff_t *local_M,
    int *m, int *window,
    double *x_max, int *np);
static void compare_f(
    const pnfft_complex *f_pnfft, const pnfft_complex *f_nfft, ptrdiff_t local_M,
    double f_hat_sum, const char *name, MPI_Comm comm_cart_2d);


int main(int argc, char **argv){
  int np[2], m, window;
  unsigned window_flag;
  ptrdiff_t N[3], n[3], local_M;
  double f_hat_sum, x_max[3];
  pnfft_complex *f1, *f2;
  
  MPI_Init(&argc, &argv);
  pnfft_init();
  
  /* set default values */
  N[0] = N[1] = N[2] = 16;
  n[0] = n[1] = n[2] = 0;
  local_M = 0;
  m = 6;
  window = 4;
  x_max[0] = x_max[1] = x_max[2] = 0.5;
  np[0]=2; np[1]=1;
  
  /* set parameters by command line */
  init_parameters(argc, argv, N, n, &local_M, &m, &window, x_max, np);

  /* if M or n are set to zero, we choose nice values */
  local_M = (local_M==0) ? N[0]*N[1]*N[2]/(np[0]*np[1]) : local_M;
  for(int t=0; t<3; t++)
    n[t] = (n[t]==0) ? 2*N[t] : n[t];

  switch(window){
    case 0: window_flag = PNFFT_WINDOW_GAUSSIAN; break;
    case 1: window_flag = PNFFT_WINDOW_BSPLINE; break;
    case 2: window_flag = PNFFT_WINDOW_SINC_POWER; break;
    case 3: window_flag = PNFFT_WINDOW_BESSEL_I0; break;
    default: window_flag = PNFFT_WINDOW_KAISER_BESSEL;
  }

  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of parallel NFFT\n");
  pfft_printf(MPI_COMM_WORLD, "* for  N[0] x N[1] x N[2] = %td x %td x %td Fourier coefficients (change with -pnfft_N * * *)\n", N[0], N[1], N[2]);
  pfft_printf(MPI_COMM_WORLD, "* at   local_M = %td nodes per process (change with -pnfft_local_M *)\n", local_M);
  pfft_printf(MPI_COMM_WORLD, "* with n[0] x n[1] x n[2]= %td x %td x %td FFT grid size (change with -pnfft_n * * *),\n", n[0], n[1], n[2]);
  pfft_printf(MPI_COMM_WORLD, "*      m = %td real space cutoff (change with -pnfft_m *),\n", m);
  pfft_printf(MPI_COMM_WORLD, "*      window = %d window function ", window);
  switch(window){
    case 0: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_GAUSSIAN) "); break;
    case 1: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_BSPLINE) "); break;
    case 2: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_SINC_POWER) "); break;
    case 3: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_BESSEL_I0) "); break;
    default: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_KAISER_BESSEL) "); break;
  }
  pfft_printf(MPI_COMM_WORLD, "(change with -pnfft_window *),\n");
  pfft_printf(MPI_COMM_WORLD, "* on   np[0] x np[1] = %td x %td processes (change with -pnfft_np * *)\n", np[0], np[1]);
  pfft_printf(MPI_COMM_WORLD, "*******************************************************************************************************\n\n");


  /* calculate parallel NFFT */
  pnfft_perform_guru(N, n, local_M, m,   x_max, window_flag, np, MPI_COMM_WORLD,
      &f1, &f_hat_sum);

  /* calculate parallel NFFT with higher accuracy */
  pnfft_perform_guru(N, n, local_M, m+2, x_max, window_flag, np, MPI_COMM_WORLD,
      &f2, &f_hat_sum);

  /* calculate error of PNFFT */
  compare_f(f1, f2, local_M, f_hat_sum, "* Results in", MPI_COMM_WORLD);

  /* free mem and finalize */
  pnfft_free(f1); pnfft_free(f2);
  pnfft_cleanup();
  MPI_Finalize();
  return 0;
}


static void pnfft_perform_guru(
    const ptrdiff_t *N, const ptrdiff_t *n, ptrdiff_t local_M,
    int m, const double *x_max, unsigned window_flag,
    const int *np, MPI_Comm comm,
    pnfft_complex **f, double *f_hat_sum
    )
{
  int myrank;
  ptrdiff_t local_N[3], local_N_start[3];
  double lower_border[3], upper_border[3];
  double local_sum = 0;
  MPI_Comm comm_cart_2d;
  pnfft_complex *f_hat;
  double *x;
  pnfft_plan pnfft;

  /* create three-dimensional process grid of size np[0] x np[1], if possible */
  if( pnfft_create_procmesh(2, comm, np, &comm_cart_2d) ){
    pfft_fprintf(comm, stderr, "Error: Procmesh of size %d x %d does not fit to number of allocated processes.\n", np[0], np[1]);
    pfft_fprintf(comm, stderr, "       Please allocate %d processes (mpiexec -np %d ...) or change the procmesh (with -pnfft_np * *).\n", np[0]*np[1], np[0]*np[1]);
    MPI_Finalize();
    return;
  }

  MPI_Comm_rank(comm, &myrank);

  /* get parameters of data distribution */
  pnfft_local_size_guru(3, N, n, x_max, m, comm_cart_2d, PNFFT_TRANSPOSED_F_HAT,
      local_N, local_N_start, lower_border, upper_border);

  /* plan parallel NFFT */
  pnfft = pnfft_init_guru(3, N, n, x_max, local_M, m,
      PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F| PNFFT_TRANSPOSED_F_HAT| window_flag, PFFT_ESTIMATE,
      comm_cart_2d);

  /* get data pointers */
  f_hat = pnfft_get_f_hat(pnfft);
  *f    = pnfft_get_f(pnfft);
  x     = pnfft_get_x(pnfft);

  /* initialize Fourier coefficients */
  pnfft_init_f_hat_3d(N, local_N, local_N_start, PNFFT_TRANSPOSED_F_HAT,
      f_hat);

  /* initialize nonequispaced nodes */
  pnfft_init_x_3d_adv(lower_border, upper_border, x_max, local_M,
      x);

  if(myrank==0)
    x[0] = -0.1;
  else
    x[0] = 0.1;

  /* execute parallel NFFT */
  pnfft_trafo(pnfft);

  /* calculate norm of Fourier coefficients for calculation of relative error */ 
  for(ptrdiff_t k=0; k<local_N[0]*local_N[1]*local_N[2]; k++)
    local_sum += cabs(f_hat[k]);
  MPI_Allreduce(&local_sum, f_hat_sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart_2d);

  /* free mem and finalize, do not free nfft.f */
  pnfft_finalize(pnfft, PNFFT_FREE_X | PNFFT_FREE_F_HAT);
  MPI_Comm_free(&comm_cart_2d);
}


static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *N, ptrdiff_t *n, ptrdiff_t *local_M,
    int *m, int *window,
    double *x_max, int *np
    )
{
  pfft_get_args(argc, argv, "-pnfft_local_M", 1, PFFT_PTRDIFF_T, local_M);
  pfft_get_args(argc, argv, "-pnfft_N", 3, PFFT_PTRDIFF_T, N);
  pfft_get_args(argc, argv, "-pnfft_n", 3, PFFT_PTRDIFF_T, n);
  pfft_get_args(argc, argv, "-pnfft_np", 2, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pnfft_m", 1, PFFT_INT, m);
  pfft_get_args(argc, argv, "-pnfft_window", 1, PFFT_INT, window);
  pfft_get_args(argc, argv, "-pnfft_x_max", 3, PFFT_DOUBLE, x_max);
}


static void compare_f(
    const pnfft_complex *f_pnfft, const pnfft_complex *f_nfft, ptrdiff_t local_M,
    double f_hat_sum, const char *name, MPI_Comm comm
    )
{
  double error = 0, error_max;

  for(ptrdiff_t j=0; j<local_M; j++)
    if( cabs(f_pnfft[j]-f_nfft[j]) > error)
      error = cabs(f_pnfft[j]-f_nfft[j]);

  MPI_Reduce(&error, &error_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  pfft_printf(comm, "%s absolute error = %6.2e\n", name, error_max);
  pfft_printf(comm, "%s relative error = %6.2e\n", name, error_max/f_hat_sum);
}

