#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>

#include <fftw3.h>
#include <pfft.h>
#include <pnfft.h>

void message(char *string, MPI_Comm comm);

static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, double complex *data);
static double compare_complex(const ptrdiff_t *local_N_c2r,
    const double complex *data_c2c, const double complex *data_c2r,
    MPI_Comm comm);
static double check_output_c2c(const ptrdiff_t *N, 
    const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, 
    const double complex *data, MPI_Comm comm);

int main(int argc, char **argv)
{
  int np[3];
  ptrdiff_t N[3];
  ptrdiff_t local_M;

  ptrdiff_t local_N_c2c[3], local_N_start_c2c[3];
  double lower_border_c2c[3], upper_border_c2c[3];
  MPI_Comm comm_cart_3d_c2c;
  pnfft_plan plan_c2c;
  pnfft_complex *f_hat_c2c, *f_c2c;
  double *x_c2c;

  ptrdiff_t local_N_c2r[3], local_N_start_c2r[3];
  double lower_border_c2r[3], upper_border_c2r[3];
  MPI_Comm comm_3d_c2r, comm_cart_3d_c2r;
  pnfft_plan plan_c2r;
  pnfft_complex *f_hat_c2r, *f_c2r;
  double *x_c2r;

  int is_c2r;
  
  /* Initialize MPI and PFFT */
  MPI_Init(&argc, &argv);
  pnfft_init();
  
  np[0] = 2; np[1] = 2; np[2] = 4;
  // N[2] = 35 is chosen so that the data distribution is the same for c2c and c2r on the first 2 processes of np[2]
  N[0] = N[1] = 8; N[2] = 35;
  local_M = N[0]*N[1]*N[2]/(np[0]*np[1]*np[2]);
  
   /* Print infos */
  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of parallel NFFT\n");
  pfft_printf(MPI_COMM_WORLD, "* for  N[0] x N[1] x N[2] = %td x %td x %td Fourier coefficients)\n", N[0], N[1], N[2]);
  pfft_printf(MPI_COMM_WORLD, "* at   local_M = %td nodes per process\n", local_M);
  pfft_printf(MPI_COMM_WORLD, "* on   np[0] x np[1] x np[2] = %td x %td x %td processes\n", np[0], np[1], np[2]);
  pfft_printf(MPI_COMM_WORLD, "*******************************************************************************************************\n\n");

  /* create three-dimensional process grid of size np[0] x np[1] x np[2], if possible */
  if( pnfft_create_procmesh(3, MPI_COMM_WORLD, np, &comm_cart_3d_c2c) ){
    pfft_fprintf(MPI_COMM_WORLD, stderr, "Error: Procmesh of size %d x %d x %d does not fit to number of allocated processes.\n", np[0], np[1], np[2]);
    pfft_fprintf(MPI_COMM_WORLD, stderr, "       Please allocate %d processes (mpiexec -np %d ...) or change the procmesh (with -pnfft_np * * *).\n", np[0]*np[1]*np[2], np[0]*np[1]*np[2]);
    MPI_Finalize();
    return 1;
  }
 
  /* Get parameters of data distribution */
  pnfft_local_size_3d(N, comm_cart_3d_c2c, PNFFT_TRANSPOSED_NONE,
      local_N_c2c, local_N_start_c2c, lower_border_c2c, upper_border_c2c);

  printf("c2c: locN: %d, %d, %d, locNstart: %d, %d, %d\n", local_N_c2c[0], local_N_c2c[1], local_N_c2c[2], local_N_start_c2c[0], local_N_start_c2c[1], local_N_start_c2c[2]);

  np[2] /= 2;
//  if (local_N_start_c2c[2] < 0) {
  if (0) {
    is_c2r = 1;
    MPI_Comm_split(MPI_COMM_WORLD, 1 , 1, &comm_3d_c2r);
  } else {
    is_c2r = 0;
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 1, &comm_3d_c2r);
  }

  /* Plan parallel NFFT */
//  plan_c2c = pnfft_init_3d(N, local_M, comm_cart_3d_c2c);
  plan_c2c = pnfft_init_adv(3, N, local_M,
      PNFFT_TRANSPOSED_NONE| PNFFT_WINDOW_SINC_POWER| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
      comm_cart_3d_c2c);

  f_hat_c2c = pnfft_get_f_hat(plan_c2c);
  f_c2c     = pnfft_get_f(plan_c2c);
  x_c2c     = pnfft_get_x(plan_c2c);

  /* Initialize Fourier coefficients with random numbers */
  init_input(N, local_N_c2c, local_N_start_c2c, f_hat_c2c);
  
  /* Initialize nodes with random numbers */
  pnfft_init_x_3d(lower_border_c2c, upper_border_c2c, local_M, x_c2c);


  if ((local_N_start_c2c[2] == -17) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);
  if ((local_N_start_c2c[2] == -8) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);
  if ((local_N_start_c2c[2] == 1) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);
  if ((local_N_start_c2c[2] == 10) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);




  if (is_c2r) {
    pnfft_create_procmesh(3, comm_3d_c2r, np, &comm_cart_3d_c2r);

    pnfft_local_size_3d_c2r(N, comm_cart_3d_c2r, PNFFT_TRANSPOSED_NONE,
        local_N_c2r, local_N_start_c2r, lower_border_c2r, upper_border_c2r);

    printf("c2r: locN: %d, %d, %d, locNstart: %d, %d, %d\n", local_N_c2r[0], local_N_c2r[1], local_N_c2r[2], local_N_start_c2r[0], local_N_start_c2r[1], local_N_start_c2r[2]);

    plan_c2r = pnfft_init_adv_c2r(3, N, local_M,
        PNFFT_TRANSPOSED_NONE| PNFFT_WINDOW_SINC_POWER| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
        comm_cart_3d_c2r);

    f_hat_c2r = pnfft_get_f_hat(plan_c2r);
    f_c2r     = pnfft_get_f(plan_c2r);
    x_c2r     = pnfft_get_x(plan_c2r);

    init_input(N, local_N_c2r, local_N_start_c2r, f_hat_c2r);

    for (int k=0; k<local_M*3; k++)
      x_c2r = x_c2c;
 
    double err = compare_complex(local_N_c2r, f_hat_c2c, f_hat_c2r, comm_cart_3d_c2r);
    pfft_printf(MPI_COMM_WORLD, "max error between c2c and c2r before trafo: %6.2e\n", err);
  }


  /* execute parallel NFFT */
  pnfft_trafo(plan_c2c);
  
  /* Print transformed data */
//  pnfft_vpr_real(x, 3*local_M, "PNFFT, x", MPI_COMM_WORLD);
//  pnfft_vpr_complex(f, local_M, "PNFFT, f", MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  /* execute parallel adjoint NFFT */
  pnfft_adj(plan_c2c);
  
  /* Print result of adjoint NFFT */
//  pnfft_apr_complex_3d(
//      f_hat, local_N, local_N_start, 0, "PNFFT^H, f_hat", MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  for(ptrdiff_t l=0; l < local_N_c2c[0] * local_N_c2c[1] * local_N_c2c[2]; l++)
    f_hat_c2c[l] /= (N[0]*N[1]*N[2]);

  double err;
  if(is_c2r) {
    err = compare_complex(local_N_c2r, f_hat_c2c, f_hat_c2r, comm_cart_3d_c2r);
    pfft_printf(MPI_COMM_WORLD, "max error between c2c and c2r after trafo: %6.2e\n", err);
  }
//  err = check_output_c2c(N, local_N_c2c, local_N_start_c2c, f_hat_c2c, comm_cart_3d_c2c);
//  pfft_printf(MPI_COMM_WORLD, "max error between c2c and c2r after trafo: %6.2e\n", err);


  if ((local_N_start_c2c[2] == -17) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);
  if ((local_N_start_c2c[2] == -8) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);
  if ((local_N_start_c2c[2] == 1) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);
  if ((local_N_start_c2c[2] == 10) && (local_N_start_c2c[1] == -4) && (local_N_start_c2c[0] == -4))
    for(ptrdiff_t l=0; l < local_N_c2c[2]; l++)
      printf("f_hat_c2c[%6d]: % 6.2e  + I * % 6.2e\n", l, creal(f_hat_c2c[l]), cimag(f_hat_c2c[l]));
  MPI_Barrier(MPI_COMM_WORLD);





  /* free mem and finalize */
  pnfft_finalize(plan_c2c, PNFFT_FREE_X | PNFFT_FREE_F_HAT | PNFFT_FREE_F);
  MPI_Comm_free(&comm_cart_3d_c2c);
  if (is_c2r) {
    pnfft_finalize(plan_c2r, PNFFT_FREE_X | PNFFT_FREE_F_HAT | PNFFT_FREE_F);
    MPI_Comm_free(&comm_3d_c2r);
    MPI_Comm_free(&comm_cart_3d_c2r);
  }

  /* Finalize MPI */
  MPI_Finalize();
  return 0;
}

#define DATA_INIT(i) (( (double)1000 ) / ( (double)( (i) == 0 ? 1 : i) ))

static void init_input(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, double complex *data)
{
  int m = 0;
  for(ptrdiff_t k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
    for(ptrdiff_t k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
      for(ptrdiff_t k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++, m++) {
        if ((N[2]%2 == 0) && (k2 == -N[2]/2))
            data[m] = 0;
        data[m] = DATA_INIT(k0+k1)+abs(k2) + I*DATA_INIT(k0+k1)*k2;
      }
}

static double compare_complex(const ptrdiff_t *local_N_c2r,
    const double complex *data_c2c, const double complex *data_c2r,
    MPI_Comm comm)
{
  double err = 0, max_err = 0;
  double glob_max_err, re, im;

  ptrdiff_t loc_n = pfft_prod_INT(3, local_N_c2r);

  for (int k=0; k<loc_n; k++) {
    re = creal(data_c2c[k]) - creal(data_c2r[k]);
    im = cimag(data_c2c[k]) - cimag(data_c2r[k]);
    err = sqrt(re*re + im*im);
    if (err > max_err)
      max_err = err;

    MPI_Barrier(comm);

    if (err > 1)
      printf("k: %6d, data_c2c[k]: % 6.2e  + I * % 6.2e\t, data_c2r[k]: % 6.2e  + I * % 6.2e\n", k, creal(data_c2c[k]), cimag(data_c2c[k]), creal(data_c2r[k]), cimag(data_c2r[k]));
  }

  MPI_Allreduce(&max_err, &glob_max_err, 1, MPI_DOUBLE, MPI_MAX, comm);
  return glob_max_err;
}

static double check_output_c2c(const ptrdiff_t *N, const ptrdiff_t *local_N, const ptrdiff_t *local_N_start, const double complex *data, MPI_Comm comm) {
  double complex want;
  double glob_max_err, re, im;
  double err = 0, max_err = 0;

  int m = 0;
  for(ptrdiff_t k0=local_N_start[0]; k0<local_N_start[0]+local_N[0]; k0++)
    for(ptrdiff_t k1=local_N_start[1]; k1<local_N_start[1]+local_N[1]; k1++)
      for(ptrdiff_t k2=local_N_start[2]; k2<local_N_start[2]+local_N[2]; k2++, m++) {
        if ((N[2]%2 == 0) && (k2 == -N[2]/2))
            want = 0;
        want = DATA_INIT(k0+k1)+abs(k2) + I*DATA_INIT(k0+k1)*k2;

        re = creal(data[m]) - creal(want);
        im = cimag(data[m]) - cimag(want);
        err = sqrt(re*re + im*im);
        if (err > max_err)
          max_err = err;

//        MPI_Barrier(comm);

        if (err > 1)
          printf("m: %6d, data[m]: % 6.2e  + I * % 6.2e\t, data_want[m]: % 6.2e  + I * % 6.2e\n", m, creal(data[m]), cimag(data[m]), creal(want), cimag(want));
      }
  MPI_Allreduce(&max_err, &glob_max_err, 1, MPI_DOUBLE, MPI_MAX, comm);
  return glob_max_err;
}

