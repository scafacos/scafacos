
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "sl_pepckeys.h"
#include "sl_pepcparts.h"

#include "fortran2c_types.h"


typedef FINT_TYPE_C finteger_t;
#define finteger_mpi  FINT_TYPE_MPI
#define finteger_fmt  FINT_TYPE_FMT

typedef pepckeys_slint_t slint_t;
#define slint_fmt pepckeys_sl_int_type_fmt


void slcheck_fortran2c_types_(double *);

#pragma weak slcheck_fortran2c_types_ = slcheck_fortran2c_types
void slcheck_fortran2c_types(double *f2c_sizes)
{
  int error = 0;

  if (f2c_sizes[0] != sizeof(FINT_TYPE_C)) { fprintf(stderr, "WARNING: fortran integer = %d vs. FINT_TYPE_C = %d\n", (int) f2c_sizes[0], (int) sizeof(FINT_TYPE_C)); ++error; }
  if (f2c_sizes[1] != sizeof(FINT8_TYPE_C)) { fprintf(stderr, "WARNING: fortran integer*8 = %d vs. FINT8_TYPE_C = %d\n", (int) f2c_sizes[1], (int) sizeof(FINT8_TYPE_C)); ++error; }
  if (f2c_sizes[2] != sizeof(FREAL8_TYPE_C)) { fprintf(stderr, "WARNING: fortran real*8 = %d vs. FREAL8_TYPE_C = %d\n", (int) f2c_sizes[2], (int) sizeof(FREAL8_TYPE_C)); ++error; }

  if (error) fprintf(stderr, "WARNING: There seems to be a problem between Fortran and C data types. Please adjust file 'fortran2c_types.h'!\n");
}


void receive_stats(finteger_t nmax, int *scounts, int *sdispls, int *rcounts, int *rdispls, pepckeys_sldata0_t *work, int size, int rank, MPI_Comm comm)
{
  slint_t i, j;
  double w, sweights[size], rweights[size];


  printf("%d: slsort_parts: total receive count: %d (nmax: %" finteger_fmt ")%s\n", rank, rdispls[size - 1] + rcounts[size - 1], nmax, (rdispls[size - 1] + rcounts[size - 1] > nmax)?" ERROR!!!":"");

  for (j = 0; j < size; ++j)
  {
    sweights[j] = 0;
    for (i = sdispls[j]; i < sdispls[j]+scounts[j]; ++i) sweights[j] += work[i];
  }
  MPI_Alltoall(sweights, 1, MPI_DOUBLE, rweights, 1, MPI_DOUBLE, comm);
  w = 0;
  for (i = 0; i < size; ++i) w += rweights[i];

/*  printf("%d: slsort_parts: total receive weight: %f\n", rank, w);*/

  double cw[2], cw_min[2], cw_max[2], cw_sum[2];
  cw[0] = rdispls[size - 1] + rcounts[size - 1];
  cw[1] = w;
  MPI_Reduce(cw, cw_min, 2, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(cw, cw_max, 2, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(cw, cw_sum, 2, MPI_DOUBLE, MPI_SUM, 0, comm);

  if (rank == 0) printf("%d: slsort_parts: receive stats: %d  %d  %d  /  %f  %f  %f\n", rank, (int) cw_min[0], (int) cw_max[0], (int) (cw_sum[0] / size), cw_min[1], cw_max[1], cw_sum[1] / size);
}


void border_stats(slint_t nkeys, pepckeys_slkey_t *keys, int size, int rank, MPI_Comm comm)
{
  pepcparts_sl_key_type_c lmm[2], gmm[2 * size];
  
  slint_t nb, i;
  double b, bsum, bmin, bmax;

  if (nkeys > 0)
  {
    lmm[0] = keys[0];
    lmm[1] = keys[nkeys - 1];

  } else lmm[0] = lmm[1] = 0;

  MPI_Gather(lmm, 2 * pepcparts_sl_key_size_mpi, pepcparts_sl_key_type_mpi, gmm, 2 * pepcparts_sl_key_size_mpi, pepcparts_sl_key_type_mpi, 0, comm);
  
  if (rank == 0)
  {
    nb = 0;
    bsum = 0.0;
    bmin = 100.0;
    bmax = 0.0;

/*    for (i = 0; i < size; ++i)
    {
      if (i > 0) printf("%" pepckeys_sl_int_type_fmt "  %llX -> %llX\n", i, gmm[i * 2 + 0], gmm[i * 2 - 1] ^ gmm[i * 2 + 0]);
      else printf("%" pepckeys_sl_int_type_fmt "  %llX\n", i, gmm[i * 2 + 0]);
      printf("%" pepckeys_sl_int_type_fmt "  %llX\n", i, gmm[i * 2 + 1]);
    }*/

/*    printf("%d: borders:\n", rank);*/
    for (i = 0; i < size - 1; ++i)
    {
      if (gmm[i * 2 + 1] > 0 && gmm[(i + 1) * 2] > 0 && (gmm[i * 2 + 1] ^ gmm[(i + 1) * 2]) > 0)
      {
        ++nb;
        b = log((double) (gmm[i * 2 + 1] ^ gmm[(i + 1) * 2])) / log(2.0);
/*        printf("%" pepckeys_sl_int_type_fmt "  %" pepcparts_sl_key_type_fmt "  %f\n", i, gmm[i * 2 + 1] ^ gmm[(i + 1) * 2], b);*/

      } else b = 0;
      
      bsum += b;
      if (b < bmin) bmin = b;
      if (b > bmax) bmax = b;
    }

    printf("%d: borders: avg: %f, min: %f, max: %f\n", rank, (nb)?(bsum / (double) nb):-1.0, bmin, bmax);
  }
}
