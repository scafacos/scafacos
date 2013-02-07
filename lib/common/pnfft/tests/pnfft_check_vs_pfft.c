#include <stdlib.h>
#include <complex.h>
#include <pnfft.h>

double pnfft_bessel_i0(double x);

static void init_equispaced_x(
    const ptrdiff_t *N, const double *lo, const double *up,
    double *x);
static void init_input_c2c_3d(
    const ptrdiff_t *n, const ptrdiff_t *local_n, const ptrdiff_t *local_start,
    pnfft_complex *data);
static double check_output_c2c_3d(
    const ptrdiff_t *n, const ptrdiff_t *local_n, const ptrdiff_t *local_start,
    const pnfft_complex *data, MPI_Comm comm);
static void vpr_complex(
    MPI_Comm comm, ptrdiff_t num, pnfft_complex *vec, const char *info);

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *N, ptrdiff_t *n,
    int *m, int *window, int *np);

int main(int argc, char **argv){
  int np[3], m, window;
  unsigned window_flag;
  ptrdiff_t N[3], n[3], local_M;
  ptrdiff_t local_N[3], local_N_start[3];
  double lower_border[3], upper_border[3];
  double f_hat_sum, local_f_hat_sum, err;
  MPI_Comm comm_cart_3d;
  pnfft_complex *f_hat, *f;
  double *x, x_max[3];
  pnfft_plan pnfft;
  
  MPI_Init(&argc, &argv);
  pnfft_init();
  
  /* set default values */
  N[0] = N[1] = N[2] = 16;
  n[0] = n[1] = n[2] = 0;
  m = 6;
  x_max[0] = x_max[1] = x_max[2] = 0.5;
  window = 0;
  np[0]=2; np[1]=2; np[2]=2;
  
  /* set parameters by command line */
  init_parameters(argc, argv, N, n, &m, &window, np);

  /* if n is set to zero, we choose nice values */
  for(int t=0; t<3; t++)
    n[t] = (n[t]==0) ? 2*N[t] : n[t];

  switch(window){
    case 1: window_flag = PNFFT_WINDOW_GAUSSIAN; break;
    case 2: window_flag = PNFFT_WINDOW_BSPLINE; break;
    case 3: window_flag = PNFFT_WINDOW_SINC_POWER; break;
    case 4: window_flag = PNFFT_WINDOW_BESSEL_I0; break;
    default: window_flag = PNFFT_WINDOW_KAISER_BESSEL;
  }

  /* Print infos */
  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of parallel NFFT\n");
  pfft_printf(MPI_COMM_WORLD, "* for  N[0] x N[1] x N[2] = %td x %td x %td Fourier coefficients (change with -pnfft_N * * *)\n", N[0], N[1], N[2]);
  pfft_printf(MPI_COMM_WORLD, "* at   M = %td nodes (equal to N[0] x N[1] x N[2]\n", N[0]*N[1]*N[2]);
  pfft_printf(MPI_COMM_WORLD, "* with n[0] x n[1] x n[2] = %td x %td x %td FFT grid size (change with -pnfft_n * * *),\n", n[0], n[1], n[2]);
  pfft_printf(MPI_COMM_WORLD, "*      m = %d real space cutoff (change with -pnfft_m *),\n", m);
  pfft_printf(MPI_COMM_WORLD, "*      window = %d window function ", window);
  switch(window){
    case 1: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_GAUSSIAN) "); break;
    case 2: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_BSPLINE) "); break;
    case 3: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_SINC_POWER) "); break;
    case 4: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_BESSEL_I0) "); break;
    default: pfft_printf(MPI_COMM_WORLD, "(PNFFT_WINDOW_KAISER_BESSEL) "); break;
  }
  pfft_printf(MPI_COMM_WORLD, "(change with -pnfft_window *),\n");
  pfft_printf(MPI_COMM_WORLD, "* on   np[0] x np[1] x np[2] = %td x %td x %td processes (change with -pnfft_np * * *)\n", np[0], np[1], np[2]);
  pfft_printf(MPI_COMM_WORLD, "*******************************************************************************************************\n\n");

  /* create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pnfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh of size %d x %d x %d does not fit to number of allocated processes.\n", np[0], np[1], np[2]);
    pfft_fprintf(MPI_COMM_WORLD, stderr, "       Please allocate %d processes (mpiexec -np %d ...) or change the procmesh (with -pnfft_np * * *).\n", np[0]*np[1]*np[2], np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }

  /* get parameters of data distribution */
  pnfft_local_size_3d(N, comm_cart_3d,
      local_N, local_N_start, lower_border, upper_border);

  local_M = local_N[0]*local_N[1]*local_N[2];

  /* plan parallel NFFT */
  pnfft = pnfft_init_guru(3, N, n, x_max, local_M, m,
                         PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F| window_flag, PFFT_ESTIMATE,
//       PNFFT_PRE_KUB_PSI| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F| window_flag, PFFT_ESTIMATE,
      comm_cart_3d);

  double t=0.1;
  pfft_printf(MPI_COMM_WORLD, "Psi(%f) = %.16e\n", t, pnfft_psi(pnfft, 0, t));
  ptrdiff_t kt=5;
  pfft_printf(MPI_COMM_WORLD, "Phi_hat(%td) = %.16e\n", kt, pnfft_phi_hat(pnfft, 0, kt));

  t=0.1;  pfft_printf(MPI_COMM_WORLD, "pnfft_bessel_i0(%f) = %.16e\n", t, pnfft_bessel_i0(t));
  t=0.0;  pfft_printf(MPI_COMM_WORLD, "pnfft_bessel_i0(%f) = %.16e\n", t, pnfft_bessel_i0(t));
  t=10.0; pfft_printf(MPI_COMM_WORLD, "pnfft_bessel_i0(%f) = %.16e\n", t, pnfft_bessel_i0(t));
  t=15.0; pfft_printf(MPI_COMM_WORLD, "pnfft_bessel_i0(%f) = %.16e\n", t, pnfft_bessel_i0(t));
  t=16.0; pfft_printf(MPI_COMM_WORLD, "pnfft_bessel_i0(%f) = %.16e\n", t, pnfft_bessel_i0(t));
  t=20.0; pfft_printf(MPI_COMM_WORLD, "pnfft_bessel_i0(%f) = %.16e\n", t, pnfft_bessel_i0(t));
  t=2000; pfft_printf(MPI_COMM_WORLD, "pnfft_bessel_i0(%f) = %.16e\n", t, pnfft_bessel_i0(t));

  /* get data pointers */
  f_hat = pnfft_get_f_hat(pnfft);
  f     = pnfft_get_f(pnfft);
  x     = pnfft_get_x(pnfft);

  /* initialize Fourier coefficients */
  init_input_c2c_3d(N, local_N, local_N_start,
      f_hat);

  /* initialize nonequispaced nodes */
  init_equispaced_x(N, lower_border, upper_border,
      x);

  /* print input Fourier coefficents */
  vpr_complex(comm_cart_3d, 8, f_hat,
      "Input Fourier coefficients on process 1:");

  /* execute parallel NFFT */
  pnfft_trafo(pnfft);

  /* print NFFT results */
  vpr_complex(comm_cart_3d, 8, f,
      "PNFFT Results on process 1:");

  /* execute parallel adjoint NFFT */
  pnfft_adj(pnfft);

  /* scale data */
  for(ptrdiff_t k=0; k < local_N[0] * local_N[1] * local_N[2]; k++)
    f_hat[k] /= (N[0]*N[1]*N[2]);

  /* print output Fourier coefficents */
  vpr_complex(comm_cart_3d, 8, f_hat,
      "Fourier coefficients after one forward and backward PNFFT on process 1:");

  /* calculate norm of Fourier coefficients for calculation of relative error */ 
  local_f_hat_sum = 0.0;
  for(ptrdiff_t k=0; k<local_N[0]*local_N[1]*local_N[2]; k++)
    local_f_hat_sum += cabs(f_hat[k]);
  MPI_Allreduce(&local_f_hat_sum, &f_hat_sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart_3d);

  err = check_output_c2c_3d(N, local_N, local_N_start, f_hat, comm_cart_3d);
  pfft_printf(comm_cart_3d, "\nError after one forward and backward PNFFT of size n=(%td, %td, %td) with equispaced nodes:\n", N[0], N[1], N[2]); 
  pfft_printf(comm_cart_3d, "absolute maximum error = %6.2e\n", err);
  pfft_printf(comm_cart_3d, "relative maximum error = %6.2e\n", err/f_hat_sum);


  /* free mem and finalize */
  pnfft_finalize(pnfft, PNFFT_FREE_X | PNFFT_FREE_F_HAT| PNFFT_FREE_F);
  MPI_Comm_free(&comm_cart_3d);
  pnfft_cleanup();
  MPI_Finalize();
  return 0;
}

static void init_input_c2c_3d(
    const ptrdiff_t *n, const ptrdiff_t *local_n, const ptrdiff_t *local_start,
    pnfft_complex *data
    )
{
  ptrdiff_t m=0, glob_ind;

  for(ptrdiff_t k0=local_start[0]; k0<local_start[0]+local_n[0]; k0++)
    for(ptrdiff_t k1=local_start[1]; k1<local_start[1]+local_n[1]; k1++)
      for(ptrdiff_t k2=local_start[2]; k2<local_start[2]+local_n[2]; k2++, m++){
        glob_ind = (k0+n[0]/2)*n[1]*n[2] + (k1+n[1]/2)*n[2] + (k2+n[2]/2);
        data[m] = 1000.0/(2*glob_ind+1) + 1000.0/(2*glob_ind+2)*I;
      }
}

static double check_output_c2c_3d(
    const ptrdiff_t *n, const ptrdiff_t *local_n, const ptrdiff_t *local_start,
    const pnfft_complex *data, MPI_Comm comm
    )
{
  ptrdiff_t m=0, glob_ind;
  double err=0, max_err=0, glob_max_err, im, re;

  for(ptrdiff_t k0=local_start[0]; k0<local_start[0]+local_n[0]; k0++)
    for(ptrdiff_t k1=local_start[1]; k1<local_start[1]+local_n[1]; k1++)
      for(ptrdiff_t k2=local_start[2]; k2<local_start[2]+local_n[2]; k2++, m++){
        glob_ind = (k0+n[0]/2)*n[1]*n[2] + (k1+n[1]/2)*n[2] + (k2+n[2]/2);
        re = creal(data[m]) - 1000.0/(2*glob_ind+1);
        im = cimag(data[m]) - 1000.0/(2*glob_ind+2);
        err = sqrt(re*re + im*im);
        if( err > max_err )
          max_err = err;
      }

  MPI_Allreduce(&max_err, &glob_max_err, 1, MPI_DOUBLE, MPI_MAX, comm);

  return glob_max_err;
}

static void vpr_complex(
    MPI_Comm comm, ptrdiff_t num, pnfft_complex *vec, const char *info
    )
{
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  if(myrank==0){
    printf("%s\n", info);
    for(ptrdiff_t k=0; k<num; k++){
      if(k%4 == 0)
        printf("%4td. ", k);
      printf("%.2e + %.2e I\t", creal(vec[k]), cimag(vec[k]));
      if(k%4 == 3)
        printf("\n");
    }
  }
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

static void init_parameters(
    int argc, char **argv,
    ptrdiff_t *N, ptrdiff_t *n,
    int *m, int *window, int *np
    )
{
  pfft_get_args(argc, argv, "-pnfft_N", 3, PFFT_PTRDIFF_T, N);
  pfft_get_args(argc, argv, "-pnfft_n", 3, PFFT_PTRDIFF_T, n);
  pfft_get_args(argc, argv, "-pnfft_np", 3, PFFT_INT, np);
  pfft_get_args(argc, argv, "-pnfft_m", 1, PFFT_INT, m);
  pfft_get_args(argc, argv, "-pnfft_window", 1, PFFT_INT, window);
}

