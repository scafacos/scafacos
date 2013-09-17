
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define Z_MOP(_mop_)  do { _mop_ } while (0)
#define Z_NOP()       Z_MOP()

#include "config_pepc_sort.h"

#include "sl_pepckeys.h"


typedef pepckeys_slint_t slint_t;
#define slint_fmt pepckeys_sl_int_type_fmt

typedef FINT_TYPE_C finteger_t;
#define finteger_mpi  FINT_TYPE_MPI
#define finteger_fmt  FINT_TYPE_FMT

/*#define MAX_IMBALANCE  0.01*/

#define PART_MINMAX

/*#define MPI_PARTITION_RADIX_2GROUPS*/
/*#define MPI_PARTITION_RADIX_NGROUPS  2*/
/*#define MPI_PARTITION_SAMPLE*/

/*#define SORT_INSTEAD_OF_MERGE*/

/*#define VALIDATE*/

#define DO_TIMING_SYNC

#define SORT_RHIGH   -1
#define SORT_RLOW    -1
#define SORT_RWIDTH  -1

#define PART_RHIGH   62
#define PART_RLOW    -1
#define PART_RWIDTH  3


#ifdef DO_TIMING
# define TIMING_DECL(_decl_)       _decl_
# define TIMING_CMD(_cmd_)         Z_MOP(_cmd_)
#else
# define TIMING_DECL(_decl_)
# define TIMING_CMD(_cmd_)         Z_NOP()
#endif
#ifdef DO_TIMING_SYNC
# define TIMING_SYNC(_c_)          TIMING_CMD(MPI_Barrier(_c_);)
#else
# define TIMING_SYNC(_c_)          Z_NOP()
#endif
#define TIMING_START(_t_)          TIMING_CMD(((_t_) = MPI_Wtime());)
#define TIMING_STOP(_t_)           TIMING_CMD(((_t_) = MPI_Wtime() - (_t_));)
#define TIMING_STOP_ADD(_t_, _r_)  TIMING_CMD(((_r_) += MPI_Wtime() - (_t_));)


void receive_stats(finteger_t nmax, int *scounts, int *sdispls, int *rcounts, int *rdispls, pepckeys_sldata0_t *work, int size, int rank, MPI_Comm comm);
void border_stats(slint_t nkeys, pepckeys_slkey_t *keys, int size, int rank, MPI_Comm comm);


void slsort_keys(finteger_t *nin,                                       /* IN */
                 finteger_t *nmax,                                      /* IN */
                 pepckeys_slkey_t *keys,                                /* INOUT */
                 pepckeys_sldata0_t *work,                              /* INOUT */
                 finteger_t *balance_weight,                            /* IN */
                 double *max_imbalance,                                 /* IN */
                 finteger_t *nout,                                      /* OUT */
                 finteger_t *indxl, finteger_t *irnkl,                  /* OUT */
                 finteger_t *fscounts, finteger_t *frcounts,            /* OUT */
                 finteger_t *fsdispls, finteger_t *frdispls,            /* OUT */
                 pepckeys_slkey_t *keys2,                               /* SCRATCH */
                 finteger_t *irnkl2,                                    /* SCRATCH */
                 finteger_t *fsize, finteger_t *frank, MPI_Fint *fcomm) /* IN */
{
  int size = *fsize;
  int rank = *frank;
  MPI_Comm comm = MPI_Comm_f2c(*fcomm);

  slint_t i;

  pepckeys_elements_t s0, s1;
  pepckeys_partcond_t pc;

#ifdef MAX_IMBALANCE
  double imba = MAX_IMBALANCE;
#else
  double imba = *max_imbalance;
#endif

  int scounts[size], rcounts[size], sdispls[size], rdispls[size];

#if defined(MPI_PARTITION_RADIX_2GROUPS) || defined(MPI_PARTITION_RADIX_NGROUPS)
  const slint_t max_nsubs = 4;

  slint_t nsubs;
  MPI_Comm sub_comms[max_nsubs];
  int sub_sizes[max_nsubs], sub_ranks[max_nsubs];
#endif

#ifdef VALIDATE
  slint_t o, l;
#endif

#ifdef DO_TIMING
  double ttotal, tinitpre, tpresort, tpartition, talltoall, talltoallv, tinitpost, tpostmerge;
#endif

#if defined(DO_INFO)
  slint_t ntotal = *nin * size;
  const slint_t ndims = 4;
  slint_t dims[ndims], pos[ndims];
#endif


  TIMING_SYNC(comm); TIMING_START(ttotal);

  INFO_CMD(
    pepckeys_mpi_get_grid_properties(ndims, dims, pos, size, rank, comm);
    if (rank == 0)
      printf(INFO_PRINT_PREFIX "# np: %d, grid: %dx%dx%dx%d, n: %" finteger_fmt ", nmax: %" finteger_fmt ", ntotal: %" slint_fmt ", partitioning: %s, minmax: %s, weighted: %" finteger_fmt ", imba: %f\n",
        size, (int) dims[3], (int) dims[2], (int) dims[1], (int) dims[0], *nin, *nmax, ntotal,
#ifdef MPI_PARTITION_SAMPLE
        "rs",
#else
        "ep",
#endif
#ifdef PART_MINMAX
        "yes",
#else
        "no",
#endif
        *balance_weight, imba);
  );

  DEBUG_CMD(
    printf(DEBUG_PRINT_PREFIX "%d: slsort_keys with %d processes\n", rank, size);
    printf(DEBUG_PRINT_PREFIX "%d:  nin: %" finteger_fmt "\n", rank, *nin);
    printf(DEBUG_PRINT_PREFIX "%d:  nmax: %" finteger_fmt "\n", rank, *nmax);
    printf(DEBUG_PRINT_PREFIX "%d:  balance_weight: %" finteger_fmt "\n", rank, *balance_weight);
    printf(DEBUG_PRINT_PREFIX "%d:  max_imbalance: %f\n", rank, *max_imbalance);

    printf(DEBUG_PRINT_PREFIX "%d:  sizeof(integer) = %d\n", rank, (int) sizeof(finteger_t));
    printf(DEBUG_PRINT_PREFIX "%d:  sizeof(integer*8) = %d\n", rank, (int) sizeof(FINT8_TYPE_C));
    printf(DEBUG_PRINT_PREFIX "%d:  sizeof(pepckeys_key) = %d\n", rank, (int) sizeof(pepckeys_slkey_t));
    printf(DEBUG_PRINT_PREFIX "%d:  sizeof(pepckeys_index) = %d\n", rank, (int) sizeof(pepckeys_slindex_t));
  );

  pepckeys_SL_DEFCON(mpi.rank) = rank;

  pepckeys_mpi_datatypes_init();

  TIMING_SYNC(comm); TIMING_START(tinitpre);

  /* init pre indexes */
  for (i = 0; i < *nin; ++i) indxl[i] = i + 1;

  pepckeys_elem_set_size(&s0, *nin);
  pepckeys_elem_set_max_size(&s0, *nmax);
  pepckeys_elem_set_keys(&s0, keys);
  pepckeys_elem_set_indices(&s0, indxl);
  pepckeys_elem_set_data(&s0, work);

  TIMING_SYNC(comm); TIMING_STOP(tinitpre);

  /* pre sort local */
  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 1. sort keys\n", rank););

  TIMING_SYNC(comm); TIMING_START(tpresort);
  pepckeys_sort_radix(&s0, NULL, SORT_RHIGH, SORT_RLOW, SORT_RWIDTH);
  TIMING_SYNC(comm); TIMING_STOP(tpresort);

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 1. sort keys done\n", rank););


#ifdef VALIDATE
  l = pepckeys_elements_validate_order(&s0, 1);
#endif


  /* partitioning */
  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 2. partitioning\n", rank););

#define REDUCTION  0.25

#if defined(PART_MINMAX) && !defined(MPI_PARTITION_SAMPLE)
  if (*balance_weight == 0)
  {
    pc.pcm = SLPC_COUNTS_MM;
    pc.count_min = -(1.0 - imba);
    pc.count_max = -(1.0 + imba);

    pepckeys_SL_DEFCON(mseg.border_update_count_reduction) = REDUCTION;

  } else
  {
    pc.pcm = SLPC_COUNTS_MM|SLPC_WEIGHTS_MM;
    pc.count_min = 0;
    pc.count_max = *nmax;
    pc.weight_min = -(1.0 - imba);
    pc.weight_max = -(1.0 + imba);

    pepckeys_SL_DEFCON(mseg.border_update_weight_reduction) = REDUCTION;
  }
#else
  if (*balance_weight == 0)
  {
    pc.pcm = SLPC_COUNTS_LH;
    pc.count_low = ((double) (rank + 0) / (double) size) - (0.5 * imba / size);
    pc.count_low = -((pc.count_low > 0)?pc.count_low:0);
    pc.count_high = ((double) (rank + 1) / (double) size) + (0.5 * imba / size);
    pc.count_high = -((pc.count_high < 1)?pc.count_high:1);

  } else
  {
    pc.pcm = SLPC_WEIGHTS_LH;
    pc.weight_low = ((double) (rank + 0) / (double) size) - (0.5 * imba / size);
    pc.weight_low = -((pc.weight_low > 0)?pc.weight_low:0);
    pc.weight_high = ((double) (rank + 1) / (double) size) + (0.5 * imba / size);
    pc.weight_high = -((pc.weight_high < 1)?pc.weight_high:1);

    pc.pcm |= SLPC_COUNTS_MM;
    pc.count_min = 0;
    pc.count_max = *nmax;
  }
#endif

/*  pepckeys_mseg_finalize_mode = SL_MSEG_FM_ALLORNOTHING;*/

  TIMING_SYNC(comm); TIMING_START(tpartition);

#if defined(MPI_PARTITION_SAMPLE)

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_parts: 2. partitioning: regular sampling\n", rank););

/*  pepckeys_mss_root = -1;*/

  pepckeys_mpi_partition_sample_regular(&s0, &pc, scounts, NULL, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_2GROUPS)

  nsubs = 2;

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_parts: 2. partitioning: 2groups\n", rank););

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_2groups(&s0, &pc, sub_comms[1], NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_NGROUPS)

  nsubs = (MPI_PARTITION_RADIX_NGROUPS <= max_nsubs)?MPI_PARTITION_RADIX_NGROUPS:max_nsubs;

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_parts: 2. partitioning: ngroups (%" slint_fmt ")\n", rank, nsubs););

  pepckeys_elem_set_indices(&k0, indxl2); /* ..._ngroups requires indices that can be modified (indxl2 is not used otherwise) */

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_ngroups(&s0, &pc, nsubs, sub_comms, NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#else

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_parts: 2. partitioning: direct\n", rank););

  pepckeys_mpi_partition_exact_radix(&s0, &pc, PART_RHIGH, PART_RLOW, PART_RWIDTH, SL_SORTED_IN, scounts, NULL, size, rank, comm);

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: average finish round: %f\n", rank, pepckeys_SL_DEFCON(mseg.info_finish_rounds_avg)););

#endif

  TIMING_SYNC(comm); TIMING_STOP(tpartition);

  pepckeys_counts2displs(size, scounts, sdispls);

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 2. partitioning done\n", rank););


  /* alltoallv keys */
  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 3. alltoallv\n", rank););

  TIMING_SYNC(comm); TIMING_START(talltoall);
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  TIMING_SYNC(comm); TIMING_STOP(talltoall);
  for (rdispls[0] = 0, i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  TIMING_SYNC(comm); TIMING_START(talltoallv);
  MPI_Alltoallv(keys, scounts, sdispls, pepckeys_sl_key_type_mpi, keys2, rcounts, rdispls, pepckeys_sl_key_type_mpi, comm);
  TIMING_SYNC(comm); TIMING_STOP(talltoallv);

  *nout = rdispls[size - 1] + rcounts[size - 1];

  for (i = 0; i < size; ++i)
  {
    fscounts[i] = scounts[i];
    frcounts[i] = rcounts[i];
    fsdispls[i] = sdispls[i];
    frdispls[i] = rdispls[i];
  }

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 3. alltoallv done\n", rank););


  /* post sort local (or mergep) */
  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 4. merge keys\n", rank););

#ifdef SORT_INSTEAD_OF_MERGE

  TIMING_SYNC(comm); TIMING_START(tinitpost);

  /* init post indexes */
  for (i = 0; i < *nout; ++i) irnkl2[i] = i;

  pepckeys_elem_set_size(&s0, *nout);
  pepckeys_elem_set_max_size(&s0, *nmax);
  pepckeys_elem_set_keys(&s0, keys2);
  pepckeys_elem_set_indices(&s0, irnkl2);
  pepckeys_elem_set_data(&s0, work);

  pepckeys_elem_set_size(&s1, *nout);
  pepckeys_elem_set_max_size(&s1, *nmax);
  pepckeys_elem_set_keys(&s1, keys2);
  pepckeys_elem_set_indices(&s1, irnkl2);
  pepckeys_elem_set_data(&s1, work);

  TIMING_SYNC(comm); TIMING_STOP(tinitpost);

  TIMING_SYNC(comm); TIMING_START(tpostmerge);
  pepckeys_sort_radix(&s0, NULL, -1, -1, -1);
  TIMING_SYNC(comm); TIMING_STOP(tpostmerge);

#else

  TIMING_SYNC(comm); TIMING_START(tinitpost);

  /* init post indexes */
  for (i = 0; i < *nout; ++i) irnkl[i] = i;

  memcpy(keys, keys2, *nout * sizeof(pepckeys_slkey_t));

  pepckeys_elem_set_size(&s0, *nout);
  pepckeys_elem_set_max_size(&s0, *nmax);
  pepckeys_elem_set_keys(&s0, keys);
  pepckeys_elem_set_indices(&s0, irnkl);
  pepckeys_elem_set_data(&s0, work);

  pepckeys_elem_set_size(&s1, *nout);
  pepckeys_elem_set_max_size(&s1, *nmax);
  pepckeys_elem_set_keys(&s1, keys2);
  pepckeys_elem_set_indices(&s1, irnkl2);
  pepckeys_elem_set_data(&s1, work);

  TIMING_SYNC(comm); TIMING_STOP(tinitpost);

  TIMING_SYNC(comm); TIMING_START(tpostmerge);
  pepckeys_mergep_heap_idx(&s0, &s1, size, frdispls, frcounts);
  TIMING_SYNC(comm); TIMING_STOP(tpostmerge);

#endif

  for (i = 0; i < s0.size; ++i) irnkl[irnkl2[i]] = i + 1;

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: slsort_keys: 4. merge keys done\n", rank););


#ifdef VALIDATE
  o = pepckeys_mpi_elements_validate_order(&s1, 1, size, rank, comm);
  INFO_CMD(
    if (rank == 0) printf(INFO_PRINT_PREFIX "slsort_keys: global order: %s - local order: %s\n", (!o)?"success":"FAILED", (!l)?"success":"failed");
  );
#endif


  pepckeys_mpi_datatypes_release();

  TIMING_SYNC(comm); TIMING_STOP(ttotal);

  DEBUG_CMD(printf(DEBUG_PRINT_PREFIX "%d: nout: %" FINT_TYPE_FMT "\n", rank, *nout););

  TIMING_CMD(
    if (rank == 0)
    {
      printf(TIMING_PRINT_PREFIX "slsort_keys: %f\n", ttotal);
      printf(TIMING_PRINT_PREFIX "slsort_keys: initpre: %f\n", tinitpre);
      printf(TIMING_PRINT_PREFIX "slsort_keys: presort: %f\n", tpresort);
      printf(TIMING_PRINT_PREFIX "slsort_keys: partition: %f\n", tpartition);
      printf(TIMING_PRINT_PREFIX "slsort_keys: alltoall: %f\n", talltoall);
      printf(TIMING_PRINT_PREFIX "slsort_keys: alltoallv: %f\n", talltoallv);
      printf(TIMING_PRINT_PREFIX "slsort_keys: initpost: %f\n", tinitpost);
      printf(TIMING_PRINT_PREFIX "slsort_keys: postmerge: %f\n", tpostmerge);
    }
  );
}


void slsort_keys_(finteger_t *nin,                                       /* IN */
                  finteger_t *nmax,                                      /* IN */
                  pepckeys_slkey_t *keys,                                /* INOUT */
                  pepckeys_sldata0_t *work,                              /* INOUT */
                  finteger_t *balance_weight,                            /* IN */
                  double *max_imbalance,                                 /* IN */
                  finteger_t *nout,                                      /* OUT */
                  finteger_t *indxl, finteger_t *irnkl,                  /* OUT */
                  finteger_t *fscounts, finteger_t *frcounts,            /* OUT */
                  finteger_t *fsdispls, finteger_t *frdispls,            /* OUT */
                  pepckeys_slkey_t *keys2,                               /* SCRATCH */
                  finteger_t *irnkl2,                                    /* SCRATCH */
                  finteger_t *fsize, finteger_t *frank, MPI_Fint* fcomm) /* IN */
{
  slsort_keys(nin, nmax, keys, work, balance_weight, max_imbalance, nout, indxl, irnkl, fscounts, frcounts, fsdispls, frdispls, keys2, irnkl2, fsize, frank, fcomm);
}
