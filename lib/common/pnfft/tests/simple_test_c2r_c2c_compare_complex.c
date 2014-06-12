#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>

#include <fftw3.h>
#include <pfft.h>
#include <pnfft.h>

void message(char *string, MPI_Comm comm);

static void init_equispaced_x(const ptrdiff_t *N, const double *lo, const double *up, double *x);
static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, const int is_complex, double *data);

int main(int argc, char **argv)
{
  int np[3];
  ptrdiff_t N[3];
  ptrdiff_t local_M;

  ptrdiff_t local_N_c2c[3], local_N_start_c2c[3];
  double lower_border_c2c[3], upper_border_c2c[3];
  pnfft_plan plan_c2c;
  pnfft_complex *f_hat_c2c, *f_c2c;
  double *x_c2c;
  MPI_Comm comm_3d_c2c;

  ptrdiff_t local_N_c2r[3], local_N_start_c2r[3];
  double lower_border_c2r[3], upper_border_c2r[3];
  pnfft_plan plan_c2r;
  pnfft_complex *f_hat_c2r;
  double *x_c2r, *f_c2r;
  MPI_Comm comm_3d_c2r;

  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pnfft_init();

  // considerations:
  // because of the test layout, the resulting blocks in frequency domain
  // need to have the same size for c2c and c2r transform. To achieve that,
  // N[2] needs to be even and divisible by np[2] and np2_c2r;
  // ceil( (N[2]/2 + 1)/np2_c2r ) == N[2]/np[2];
  // thus we choose N[2] = 24 = 6*4 = np[2]*np2_c2r,
  // because N[2]/2+1 = 13 and ceil(13/4) = 4 = 24/6;
  //
  // this should have worked, but didn't.
  // So I decided to just print the results into files and compare those


  np[0] = 1; np[1] = 1; np[2] = 8;
  N[0] = 4; N[1] = 4; N[2] = 8;
  local_M = N[0]*N[1]*N[2]/(np[0]*np[1]*np[2]);

   /* Print infos */
  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of parallel NFFT\n");
  pfft_printf(MPI_COMM_WORLD, "* for  N[0] x N[1] x N[2] = %td x %td x %td Fourier coefficients)\n", N[0], N[1], N[2]);
  pfft_printf(MPI_COMM_WORLD, "* at   local_M = %td nodes per process\n", local_M);
  pfft_printf(MPI_COMM_WORLD, "* on   np[0] x np[1] x np[2] = %td x %td x %td processes\n", np[0], np[1], np[2]);
  pfft_printf(MPI_COMM_WORLD, "*******************************************************************************************************\n\n");

  /* create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pnfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_3d_c2c) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh of size %d x %d x %d does not fit to number of allocated processes.\n", np[0], np[1], np[2]);
    pfft_fprintf(MPI_COMM_WORLD, stderr, "       Please allocate %d processes (mpiexec -np %d ...) or change the procmesh (with -pnfft_np * * *).\n", np[0]*np[1]*np[2], np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }
  if( pnfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_3d_c2r) ) {
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: building c2r procmesh\n");
    MPI_Finalize();
    return 1;
  }

  /* Get parameters of data distribution */
  pnfft_local_size_3d(N, comm_3d_c2c, PNFFT_TRANSPOSED_NONE,
      local_N_c2c, local_N_start_c2c, lower_border_c2c, upper_border_c2c);
  pnfft_local_size_3d_c2r(N, comm_3d_c2r, PNFFT_TRANSPOSED_NONE,
      local_N_c2r, local_N_start_c2r, lower_border_c2r, upper_border_c2r);

  /* Plan parallel NFFT */
  plan_c2c = pnfft_init_adv(3, N, local_M,
      PNFFT_TRANSPOSED_NONE| PNFFT_WINDOW_SINC_POWER| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
      comm_3d_c2c);
  plan_c2r = pnfft_init_adv_c2r(3, N, local_M,
      PNFFT_TRANSPOSED_NONE| PNFFT_WINDOW_SINC_POWER| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
      comm_3d_c2r);

  f_c2c     = pnfft_get_f(plan_c2c);
  f_hat_c2c = pnfft_get_f_hat(plan_c2c);
  x_c2c     = pnfft_get_x(plan_c2c);
  f_c2r     = pnfft_get_f(plan_c2r);
  f_hat_c2r = pnfft_get_f_hat(plan_c2r);
  x_c2r     = pnfft_get_x(plan_c2r);

  /* Initialize Fourier coefficients with random numbers */
  init_input(N, local_N_c2c, local_N_start_c2c, 1, f_c2c);
  init_input(N, local_N_c2c, local_N_start_c2c, 0, f_c2r);

  /* Initialize nodes with random numbers */
  init_equispaced_x(N, lower_border_c2c, upper_border_c2c, x_c2c);
  //pnfft_init_x_3d(lower_border_c2c, upper_border_c2c, local_M, x_c2c);
  for (int k=0; k<local_M*3; k++)
    x_c2r[k] = x_c2c[k];

  pnfft_adj(plan_c2c);
  pnfft_adj(plan_c2r);

  for(ptrdiff_t l=0; l < local_N_c2c[0] * local_N_c2c[1] * local_N_c2c[2]; l++)
    f_hat_c2c[l] /= (N[0]*N[1]*N[2]);
  for(ptrdiff_t l=0; l < local_N_c2r[0] * local_N_c2r[1] * local_N_c2r[2]; l++)
    f_hat_c2r[l] /= (N[0]*N[1]*N[2]);

  // c2c_results-N_start-....txt
#define LEN 50
  char str[LEN];
  FILE *file_c2c, *file_c2r;
  snprintf(str, LEN, "./c2c_results-N_start-%03td_%03td_%03td.txt",
      N[0]/2 + local_N_start_c2c[0], N[1]/2 + local_N_start_c2c[1], N[2]/2 + local_N_start_c2c[2]);
  if ((file_c2c = fopen(str, "w+")) != NULL) {
    for(ptrdiff_t k2=0; k2 < local_N_c2c[2]; k2++) {
      for(ptrdiff_t k1=0; k1 < local_N_c2c[1]; k1++) {
        for(ptrdiff_t k0=0; k0 < local_N_c2c[0]; k0++) {
          pnfft_complex tmp = f_hat_c2c[k0*local_N_c2c[1]*local_N_c2c[2]+k1*local_N_c2c[2]+k2];
          fprintf(file_c2c, "% 03td\t% 03td\t% 03td\t% .2e\t% .2e\n",
              k0 + local_N_start_c2c[0], k1 + local_N_start_c2c[1], k2 + local_N_start_c2c[2],
              creal(tmp), cimag(tmp));
        }
      }
    }
    fsync(file_c2c);
    fclose(file_c2c);
  }

  snprintf(str, LEN, "./c2r_results-N_start-%03td_%03td_%03td.txt",
      N[0]/2 + local_N_start_c2c[0], N[1]/2 + local_N_start_c2c[1], N[2]/2 + local_N_start_c2c[2]);
  if ((file_c2r = fopen(str, "w+")) != NULL) {
    for(ptrdiff_t k2=0; k2 < local_N_c2r[2]; k2++) {
      for(ptrdiff_t k1=0; k1 < local_N_c2r[1]; k1++) {
        for(ptrdiff_t k0=0; k0 < local_N_c2r[0]; k0++) {
          pnfft_complex tmp = f_hat_c2r[k0*local_N_c2r[1]*local_N_c2r[2]+k1*local_N_c2r[2]+k2];
          fprintf(file_c2r, "% 03td\t% 03td\t% 03td\t% .2e\t% .2e\n",
              k0 + local_N_start_c2r[0], k1 + local_N_start_c2r[1], k2 + local_N_start_c2r[2],
              creal(tmp), cimag(tmp));
        }
      }
    }
    fsync(file_c2r);
    fclose(file_c2r);
  }

  pnfft_print_average_timer_adv(plan_c2c, comm_3d_c2c);
  pnfft_print_average_timer_adv(plan_c2r, comm_3d_c2r);


  /* free mem and finalize */
  pnfft_finalize(plan_c2c, PNFFT_FREE_X | PNFFT_FREE_F_HAT | PNFFT_FREE_F);
  pnfft_finalize(plan_c2r, PNFFT_FREE_X | PNFFT_FREE_F_HAT | PNFFT_FREE_F);
  MPI_Comm_free(&comm_3d_c2c);
  MPI_Comm_free(&comm_3d_c2r);

  /* Finalize MPI */
  MPI_Finalize();
  return 0;
}

#define DATA_INIT(i) (( (double)1000 ) / ( (double)( (i) == 0 ? 1 : i) ))

static double semirandom(const ptrdiff_t *N, ptrdiff_t k0, ptrdiff_t k1, ptrdiff_t k2) {
  if (((k0 == -N[0]/2) || (k0 == 0)) &&
      ((k1 == -N[1]/2) || (k1 == 0)) &&
      ((k2 == -N[2]/2) || (k2 == 0)))
    return 0;
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

  return re;
}

static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, const int is_complex, double *data)
{
  int m = 0;
  for(ptrdiff_t k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
    for(ptrdiff_t k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
      for(ptrdiff_t k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++, m++) {
        if (is_complex) {
          data[2*m+0] = semirandom(N, k0, k1, k2);
          data[2*m+1] = 0;
        } else {
          data[m] = semirandom(N, k0, k1, k2);
        }
      }
}

static void init_equispaced_x(
    const ptrdiff_t *N, const double *lo, const double *up,
    double *x
    )
{
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

