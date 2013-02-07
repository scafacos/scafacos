#include <stdlib.h>
#include <complex.h>
#include <pnfft.h>

static void init_random_x(
    const double *lo, const double *up, ptrdiff_t local_M,
    double *x);
static double random_number_less_than_one(
    void);
static void vpr_complex(
    MPI_Comm comm, ptrdiff_t num, pnfft_complex *vec, const char *info);


int main(int argc, char **argv){
  int np[3];
  ptrdiff_t N[3], local_M;
  ptrdiff_t local_N[3], local_N_start[3];
  double lower_border[3], upper_border[3];
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
  pfft_printf(MPI_COMM_WORLD, "* for  N[0] x N[1] x N[2] = %td x %td x %td Fourier coefficients\n", N[0], N[1], N[2]);
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
  pnfft_local_size_3d(N, comm_cart_3d,
      local_N, local_N_start, lower_border, upper_border);

  /* Plan parallel NFFT */
  pnfft = pnfft_init_3d(N, local_M, comm_cart_3d);

  /* get data pointers */
  f_hat = pnfft_get_f_hat(pnfft);
  f     = pnfft_get_f(pnfft);
  x     = pnfft_get_x(pnfft);

  /* Initialize Fourier coefficients */
  pnfft_init_f_hat_3d(N, local_N, local_N_start,
      f_hat);

  /* Initialize nonequispaced nodes */
  init_random_x(lower_border, upper_border, local_M,
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

  /* free mem and finalize */
  pnfft_finalize(pnfft, PNFFT_FREE_X | PNFFT_FREE_F_HAT| PNFFT_FREE_F);
  MPI_Comm_free(&comm_cart_3d);
  pnfft_cleanup();
  MPI_Finalize();
  return 0;
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


static void init_random_x(
    const double *lo, const double *up, ptrdiff_t local_M,
    double *x
    )
{
  double tmp;
  
  for (ptrdiff_t j=0; j<local_M; j++){
    for(int t=0; t<3; t++){
      do{
        tmp = random_number_less_than_one();
        tmp = (up[t]-lo[t]) * tmp + lo[t];
      }
      while( (tmp < -0.5) || (0.5 <= tmp) );
      x[3*j+t] = tmp;
    }
  }
}


static double random_number_less_than_one(
    void
    )
{
  double tmp;
  
  do
    tmp = ( 1.0 * rand()) / RAND_MAX;
  while(tmp>=1.0);
  
  return tmp;
}

