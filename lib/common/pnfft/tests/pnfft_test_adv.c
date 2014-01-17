#include <stdlib.h>
#include <complex.h>
#include <pnfft.h>

static void init_equispaced_x(
    const ptrdiff_t *N, const double *lo, const double *up,
    double *x);
static void init_input_complex_3d(
    const ptrdiff_t *n, const ptrdiff_t *local_n, const ptrdiff_t *local_start,
    pnfft_complex *data);
static double check_output_c2c_3d(
    const ptrdiff_t *n, const ptrdiff_t *local_n, const ptrdiff_t *local_start,
    const pnfft_complex *data, MPI_Comm comm);
static void vpr_complex(
    MPI_Comm comm, ptrdiff_t num, pnfft_complex *vec, const char *info);


int main(int argc, char **argv){
  int np[3];
  ptrdiff_t N[3], local_M;
  ptrdiff_t local_N[3], local_N_start[3];
  double lower_border[3], upper_border[3];
  double err;
  MPI_Comm comm_cart_3d;
  pnfft_complex *f_hat, *f;
  double *x;
  pnfft_plan pnfft;
  
  MPI_Init(&argc, &argv);
  pnfft_init();
  
  /* Set default values */
  N[0] = N[1] = N[2] = 16;
  np[0]=2; np[1]=2; np[2]=2;
  local_M = N[0]*N[1]*N[2]/(np[0]*np[1]*np[2]);
  
  /* Print infos */
  pfft_printf(MPI_COMM_WORLD, "******************************************************************************************************\n");
  pfft_printf(MPI_COMM_WORLD, "* Computation of parallel NFFT\n");
  pfft_printf(MPI_COMM_WORLD, "* for  N[0] x N[1] x N[2] = %td x %td x %td Fourier coefficients)\n", N[0], N[1], N[2]);
  pfft_printf(MPI_COMM_WORLD, "* at   local_M = %td nodes per process\n", local_M);
  pfft_printf(MPI_COMM_WORLD, "* on   np[0] x np[1] x np[2] = %td x %td x %td processes\n", np[0], np[1], np[2]);
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
      local_N, local_N_start, lower_border, upper_border);

  /* Plan parallel NFFT */
//   pnfft = pnfft_init_3d(N local_M, comm_cart_3d);
  pnfft = pnfft_init_adv(3, N, local_M,
//       PNFFT_WINDOW_BESSEL_I0| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
//       PNFFT_WINDOW_GAUSSIAN| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
      PNFFT_WINDOW_SINC_POWER| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
//       PNFFT_WINDOW_KAISER_BESSEL| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
//       PNFFT_WINDOW_BSPLINE| PNFFT_MALLOC_X| PNFFT_MALLOC_F_HAT| PNFFT_MALLOC_F, PFFT_ESTIMATE,
      comm_cart_3d);

  /* get data pointers */
  f_hat = pnfft_get_f_hat(pnfft);
  f     = pnfft_get_f(pnfft);
  x     = pnfft_get_x(pnfft);

  /* Initialize Fourier coefficients */
  init_input_complex_3d(N, local_N, local_N_start,
      f_hat);

  /* Initialize nonequispaced nodes */
  init_equispaced_x(N, lower_border, upper_border,
      x);

  /* Print input Fourier coefficents */
  vpr_complex(comm_cart_3d, 8, f_hat,
      "Input Fourier coefficients on process 1:");

  /* Execute parallel NFFT */
  pnfft_trafo(pnfft);

  /* Print NFFT results */
  vpr_complex(comm_cart_3d, 8, f,
      "PNFFT Results on process 1:");

  /* Execute parallel adjoint NFFT */
  pnfft_adj(pnfft);

  /* Scale data */
  for(ptrdiff_t l=0; l < local_N[0] * local_N[1] * local_N[2]; l++)
    f_hat[l] /= (N[0]*N[1]*N[2]);

  /* Print output Fourier coefficents */
  vpr_complex(comm_cart_3d, 8, f_hat,
      "Fourier coefficients after one forward and backward PNFFT on process 1:");

  err = check_output_c2c_3d(N, local_N, local_N_start, f_hat, comm_cart_3d);
  pfft_printf(comm_cart_3d, "\nError after one forward and backward PNFFT of size n=(%td, %td, %td) with equispaced nodes:\n", N[0], N[1], N[2]); 
  pfft_printf(comm_cart_3d, "maximum error = %6.2e;\n", err);


  /* free mem and finalize */
  pnfft_finalize(pnfft, PNFFT_FREE_X | PNFFT_FREE_F_HAT| PNFFT_FREE_F);
  MPI_Comm_free(&comm_cart_3d);
  pnfft_cleanup();
  MPI_Finalize();
  return 0;
}

static void init_input_complex_3d(
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


