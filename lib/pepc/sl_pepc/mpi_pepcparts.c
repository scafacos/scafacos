
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "sl_pepckeys.h"
#include "sl_pepcparts.h"

#include "fortran2c_types.h"


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

/*#define THREE_STEPS_VERSION*/
/*#define MERGE_AND_UNPACK*/

/*#define VERBOSE*/
/*#define VALIDATE*/
/*#define TIMING*/
/*#define TIMING_ROW*/

/*#define BORDER_STATS*/
/*#define RECEIVE_STATS*/


#define SORT_RHIGH   -1
#define SORT_RLOW    -1
#define SORT_RWIDTH  -1

#define PART_RHIGH   62
#define PART_RLOW    -1
#define PART_RWIDTH  3


#if defined(MPI_PARTITION_RADIX_2GROUPS) || defined(MPI_PARTITION_RADIX_NGROUPS) || defined(MPI_PARTITION_SAMPLE)
# undef THREE_STEPS_VERSION
#endif

#ifdef VERBOSE
# define VERBOSE_MOP(_mop_)  do { _mop_ ; } while (0)
#else
# define VERBOSE_MOP(_mop_)  do { } while (0)
#endif

#ifdef TIMING
# define TSTART(tid)  tid = MPI_Wtime()
# define TSTOP(tid)   tid = MPI_Wtime() - tid
#else
# define TSTART(tid)  do { } while (0)
# define TSTOP(tid)   do { } while (0)
#endif


void receive_stats(finteger_t nmax, int *scounts, int *sdispls, int *rcounts, int *rdispls, pepckeys_sldata0_t *work, int size, int rank, MPI_Comm comm);
void border_stats(slint_t nkeys, pepckeys_slkey_t *keys, int size, int rank, MPI_Comm comm);

void slsort_parts_(finteger_t *, finteger_t *, pepcparts_slkey_t *, pepcparts_sldata0_t *, pepcparts_sldata1_t *, pepcparts_sldata2_t *, pepcparts_sldata3_t *, pepcparts_sldata4_t *, pepcparts_sldata5_t *,
                   pepcparts_sldata6_t *, pepcparts_sldata7_t *, pepcparts_sldata8_t *, pepcparts_sldata9_t *, pepcparts_sldata10_t *, pepcparts_sldata11_t *, pepcparts_sldata12_t *,
                   finteger_t *, double *, finteger_t *, finteger_t *, finteger_t *, finteger_t *, finteger_t *, finteger_t *,
                   void *, void *, pepcparts_slkey_t *, pepcparts_sldata8_t *, finteger_t *, finteger_t *, finteger_t *);

#pragma weak slsort_parts_ = slsort_parts
void slsort_parts(finteger_t *n,                                                             /* INOUT */
                  finteger_t *nmax,                                                          /* IN */
                  pepcparts_slkey_t *keys,                                                   /* INOUT */
                  pepcparts_sldata0_t *x,                                                    /* INOUT */
                  pepcparts_sldata1_t *y,                                                    /* INOUT */
                  pepcparts_sldata2_t *z,                                                    /* INOUT */
                  pepcparts_sldata3_t *ux,                                                   /* INOUT */
                  pepcparts_sldata4_t *uy,                                                   /* INOUT */
                  pepcparts_sldata5_t *uz,                                                   /* INOUT */
                  pepcparts_sldata6_t *q,                                                    /* INOUT */
                  pepcparts_sldata7_t *m,                                                    /* INOUT */
                  pepcparts_sldata8_t *work,                                                 /* INOUT */
                  pepcparts_sldata9_t  *ex,                                                  /* INOUT */
                  pepcparts_sldata10_t *ey,                                                  /* INOUT */
                  pepcparts_sldata11_t *ez,                                                  /* INOUT */
                  pepcparts_sldata12_t *pelabel,                                             /* INOUT */
                  finteger_t *balance_weight,                                                /* IN */
                  double *max_imbalance,                                                     /* IN */
                  finteger_t *indxl, finteger_t *irnkl,                                      /* OUT */
                  finteger_t *fscounts, finteger_t *frcounts,                                /* OUT */
                  finteger_t *fsdispls, finteger_t *frdispls,                                /* OUT */
                  void *parts0, void *parts1,                                                /* SCRATCH */
                  pepcparts_slkey_t *keys2, pepcparts_sldata8_t *work2, finteger_t *irnkl2,  /* SCRATCH */
                  finteger_t *fsize, finteger_t *frank)                                      /* IN */
{
#ifndef NOT_sl_pepcparts

  int size = *fsize;
  int rank = *frank;
  MPI_Comm comm = MPI_COMM_WORLD;

  slint_t i;

  pepckeys_elements_t k0, k1;
  pepckeys_partcond_t pc;

  pepcparts_elements_t d0;
  pepcparts_packed_elements_t pd0;
  MPI_Datatype pdt;

  finteger_t nin = *n;

#ifdef MAX_IMBALANCE
  double imba = MAX_IMBALANCE;
#else
  double imba = *max_imbalance;
#endif

  int scounts[size], rcounts[size], sdispls[size], rdispls[size];

#define indxl2  irnkl2

#if defined(MPI_PARTITION_RADIX_2GROUPS) || defined(MPI_PARTITION_RADIX_NGROUPS)
  const slint_t max_nsubs = 4;

  slint_t nsubs;
  MPI_Comm sub_comms[max_nsubs];
  int sub_sizes[max_nsubs], sub_ranks[max_nsubs];
#endif

#ifdef VALIDATE
  slint_t o, l;
#endif

#ifdef TIMING
  double ttotal, tinitindxl, tcopy, tsort, tpartition, tpack, talltoall, talltoallv, tmakeindices;
#ifdef THREE_STEPS_VERSION
  double tunpack;
#else
  double tmergeunpack;
# ifndef MERGE_AND_UNPACK
  double tunpack;
# endif
#endif
#endif

#if defined(VERBOSE) || defined(TIMING)
  slint_t ntotal = *n * size;
  const pepcparts_slint_t ndims = 4;
  pepcparts_slint_t dims[ndims], pos[ndims];
#endif


#if defined(VERBOSE) || defined(TIMING)
  pepcparts_mpi_get_grid_properties(ndims, dims, pos, size, rank, comm);
  if (rank == 0)
    printf("# np: %d, grid: %dx%dx%dx%d, n: %" finteger_fmt ", nmax: %" finteger_fmt ", ntotal: %" slint_fmt ", partitioning: %s, minmax: %s, weighted: %" finteger_fmt ", imba: %f\n",
      size, (int) dims[3], (int) dims[2], (int) dims[1], (int) dims[0], *n, *nmax, ntotal,
# ifdef MPI_PARTITION_SAMPLE
      "rs",
# else
      "ep",
# endif
# ifdef PART_MINMAX
      "yes",
# else
      "no",
# endif
      *balance_weight, imba);
#endif

  TSTART(ttotal);

#ifdef VERBOSE
  printf("%d: slsort_parts with %d processes\n", rank, size);
  printf("%d:  n: %" finteger_fmt "\n", rank, *n);
  printf("%d:  nmax: %" finteger_fmt "\n", rank, *nmax);
  printf("%d:  balance_weight: %" finteger_fmt "\n", rank, *balance_weight);
  printf("%d:  max_imbalance: %f\n", rank, *max_imbalance);

  printf("%d:  sizeof(integer) = %d\n", rank, (int) sizeof(finteger_t));
  printf("%d:  sizeof(integer*8) = %d\n", rank, (int) sizeof(FINT8_TYPE_C));
  printf("%d:  sizeof(pepckeys_key) = %d\n", rank, (int) sizeof(pepckeys_slkey_t));
  printf("%d:  sizeof(pepckeys_index) = %d\n", rank, (int) sizeof(pepckeys_slindex_t));
  printf("%d:  sizeof(pepcparts_key) = %d\n", rank, (int) sizeof(pepcparts_slkey_t));
  printf("%d:  sizeof(pepcparts_index) = %d\n", rank, (int) sizeof(pepcparts_slindex_t));
#endif

  pepckeys_SL_DEFCON(mpi.rank) = rank;
  pepcparts_SL_DEFCON(mpi.rank) = rank;

  pepckeys_mpi_datatypes_init();
  pepcparts_mpi_datatypes_init();

  /* init pre indexes */
  TSTART(tinitindxl);
  for (i = 0; i < *n; ++i) indxl[i] = i;
  TSTOP(tinitindxl);

  pepckeys_elem_set_size(&k0, *n);
  pepckeys_elem_set_max_size(&k0, *nmax);
  pepckeys_elem_set_keys(&k0, keys);
  pepckeys_elem_set_indices(&k0, indxl);
  pepckeys_elem_set_data(&k0, work);

  pepckeys_elem_set_size(&k1, *n);
  pepckeys_elem_set_max_size(&k1, *nmax);
  pepckeys_elem_set_keys(&k1, keys2);
  pepckeys_elem_set_indices(&k1, indxl2); /* just a dummy array, indxl2 is never used */
  pepckeys_elem_set_data(&k1, work2);

  TSTART(tcopy);
  pepckeys_elements_ncopy(&k0, &k1, *n);
  TSTOP(tcopy);

#ifdef THREE_STEPS_VERSION

#ifdef VALIDATE
  l = 0;
#endif

#else /* THREE_STEPS_VERSION */

  /* pre sort keys (+ indices and work) */
  VERBOSE_MOP(printf("%d: slsort_parts: 1. sort keys\n", rank));

  TSTART(tsort);
  pepckeys_sort_radix(&k0, NULL, SORT_RHIGH, SORT_RLOW, SORT_RWIDTH);
  TSTOP(tsort);

/*  for (i = 0; i < k0.size; ++i) printf("%d  %d  %f\n", rank, (int) i, work[i]);*/

  VERBOSE_MOP(printf("%d: slsort_parts: 1. sort keys done\n", rank));


#ifdef VALIDATE
  VERBOSE_MOP(printf("%d: slsort_parts: 1. not sort keys!\n", rank));

  l = pepckeys_elements_validate_order(&k0, 1);
#endif

#endif /* THREE_STEPS_VERSION */


  /* partitioning */
  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning\n", rank));

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

/*  pepckeys_SL_DEFCON(mseg.finalize_mode) = SL_MSEG_FM_ALLORNOTHING;*/

  TSTART(tpartition);

#if defined(MPI_PARTITION_SAMPLE)

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: regular sampling\n", rank));

/*  pepckeys_SL_DEFCON(mss.root) = -1;*/

  pepckeys_mpi_partition_sample_regular(&k0, &pc, scounts, NULL, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_2GROUPS)

  nsubs = 2;

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: 2groups\n", rank));

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_2groups(&k0, &pc, sub_comms[1], NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_NGROUPS)

  nsubs = (MPI_PARTITION_RADIX_NGROUPS <= max_nsubs)?MPI_PARTITION_RADIX_NGROUPS:max_nsubs;

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: ngroups (%" slint_fmt ")\n", rank, nsubs));

  pepckeys_elem_set_indices(&k0, indxl2); /* ..._ngroups requires indices that can be modified (indxl2 is not used otherwise) */

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_ngroups(&k0, &pc, nsubs, sub_comms, NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#else

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: direct\n", rank));

#ifdef THREE_STEPS_VERSION
  pepckeys_mpi_partition_exact_radix(&k0, &pc, PART_RHIGH, PART_RLOW, PART_RWIDTH, SL_SORTED_OUT, scounts, NULL, size, rank, comm);
#else
  pepckeys_mpi_partition_exact_radix(&k0, &pc, PART_RHIGH, PART_RLOW, PART_RWIDTH, SL_SORTED_IN, scounts, NULL, size, rank, comm);
#endif

  VERBOSE_MOP(printf("%d: average finish round: %f\n", rank, pepckeys_mseg_info_finish_rounds_avg));

#endif

  TSTOP(tpartition);

  pepckeys_counts2displs(size, scounts, sdispls);

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning done\n", rank));


  /* pack */
  VERBOSE_MOP(printf("%d: slsort_parts: 3. pack\n", rank));

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys2);
/*  pepcparts_elem_set_indices(&d0, indxl); */  /* indices are not required for packaging */
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work2, ex, ey, ez, pelabel);

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts0);

  TSTART(tpack);
  pepcparts_elements_pack_indexed(&d0, &pd0, (pepcparts_slindex_t *) indxl, NULL);
  TSTOP(tpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 3. pack done\n", rank));


  /* alltoallv */
  VERBOSE_MOP(printf("%d: slsort_parts: 4. alltoallv\n", rank));

  pepcparts_mpi_elements_packed_datatype_create(&pdt, 0);

  TSTART(talltoall);
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  TSTOP(talltoall);

  for (rdispls[0] = 0, i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];

#ifdef RECEIVE_STATS
#endif

  TSTART(talltoallv);
  MPI_Alltoallv(parts0, scounts, sdispls, pdt, parts1, rcounts, rdispls, pdt, comm);
  TSTOP(talltoallv);

  pepcparts_mpi_elements_packed_datatype_destroy(&pdt);

  *n = rdispls[size - 1] + rcounts[size - 1];

  for (i = 0; i < size; ++i)
  {
    fscounts[i] = scounts[i];
    frcounts[i] = rcounts[i];
    fsdispls[i] = sdispls[i];
    frdispls[i] = rdispls[i];
  }

  VERBOSE_MOP(printf("%d: slsort_parts: 4. alltoallv done\n", rank));


#ifdef THREE_STEPS_VERSION

  VERBOSE_MOP(printf("%d: slsort_parts: 5a. sort keys\n", rank));

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts1);
  
  pepcparts_elements_unpack_keys(&pd0, keys2);

  /* init post indexes */
  for (i = 0; i < *n; ++i) irnkl2[i] = i;
  
  pepckeys_elem_set_size(&k0, *n);
  pepckeys_elem_set_max_size(&k0, *nmax);
  pepckeys_elem_set_keys(&k0, keys2);
  pepckeys_elem_set_indices(&k0, irnkl2);
  pepckeys_elem_set_data(&k0, work2);

  TSTART(tsort);
  pepckeys_sort_radix(&k0, NULL, SORT_RHIGH, SORT_RLOW, SORT_RWIDTH);
  TSTOP(tsort);

  VERBOSE_MOP(printf("%d: slsort_parts: 5b. sort keys done\n", rank));

  VERBOSE_MOP(printf("%d: slsort_parts: 5b. unpack\n", rank));

  for (i = 0; i < *n; ++i) irnkl[irnkl2[i]] = i;

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys);
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);

  TSTART(tunpack);
  pepcparts_elements_unpack_indexed(&pd0, &d0, NULL, irnkl);
  TSTOP(tunpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 5b. unpack data done\n", rank));

  TSTART(tmakeindices);
  for (i = 0; i < nin; ++i) ++indxl[i];
  for (i = 0; i < *n; ++i) ++irnkl[i];
  TSTOP(tmakeindices);

#else /* THREE_STEPS_VERSION */

#ifdef MERGE_AND_UNPACK

  /* fused merge and unpack (indices created during merge) */
  VERBOSE_MOP(printf("%d: slsort_parts: 5. merge and unpack\n", rank));

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts1);

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys);
  pepcparts_elem_set_indices(&d0, irnkl2);
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);

  TSTART(tmergeunpack);
  pepcparts_mergep_heap_unpack_idx(&pd0, &d0, size, frdispls, frcounts);
  TSTOP(tmergeunpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 5. merge and unpack done\n", rank));

  TSTART(tmakeindices);
  for (i = 0; i < nin; ++i) ++indxl[i];
  for (i = 0; i < *n; ++i) irnkl[irnkl2[i]] = i + 1;
  TSTOP(tmakeindices);

#else

  /* fused merge and unpack of keys (indices created during merge) */
  VERBOSE_MOP(printf("%d: slsort_parts: 5a. merge and unpack indices\n", rank));

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts1);

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
/*  pepcparts_elem_set_keys(&d0, keys);*/
  pepcparts_elem_set_indices(&d0, irnkl2);
/*  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);*/ 

  TSTART(tmergeunpack);
  pepcparts_mergep_heap_unpack_idxonly(&pd0, &d0, size, frdispls, frcounts);
  TSTOP(tmergeunpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 5a. merge and unpack indices done\n", rank));

  /* unpack */
  VERBOSE_MOP(printf("%d: slsort_parts: 5b. unpack\n", rank));

  for (i = 0; i < *n; ++i) irnkl[irnkl2[i]] = i;

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys);
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);

  TSTART(tunpack);
  pepcparts_elements_unpack_indexed(&pd0, &d0, NULL, irnkl);
  TSTOP(tunpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 5b. unpack data done\n", rank));

  TSTART(tmakeindices);
  for (i = 0; i < nin; ++i) ++indxl[i];
  for (i = 0; i < *n; ++i) ++irnkl[i];
  TSTOP(tmakeindices);

#endif

#endif /* THREE_STEPS_VERSION */


#ifdef VALIDATE
  o = pepcparts_mpi_elements_validate_order(&d0, 1, size, rank, comm);
 #ifdef VERBOSE
  printf("%d: slsort_parts: global order: %s - local order: %s\n", rank, (!o)?"success":"FAILED", (!l)?"success":"failed");
  printf("%d: %" pepcparts_sl_key_type_fmt " - %" pepcparts_sl_key_type_fmt "\n", rank, d0.keys[0], d0.keys[d0.size - 1]);
 #endif
#endif


  pepcparts_mpi_datatypes_release();

  TSTOP(ttotal);

  VERBOSE_MOP(printf("%d: out: n = %" FINT_TYPE_FMT "\n", rank, *n));

#ifdef BORDER_STATS
#endif

#ifdef TIMING
  if (rank == 0)
  {
# ifndef TIMING_ROW
    printf("%d: slsort_parts: %f\n", rank, ttotal);
    printf("%d: slsort_parts: initindxl: %f\n", rank, tinitindxl);
    printf("%d: slsort_parts: copy: %f\n", rank, tcopy);
#ifndef THREE_STEPS_VERSION
    printf("%d: slsort_parts: presort: %f\n", rank, tsort);
#endif
    printf("%d: slsort_parts: partition: %f\n", rank, tpartition);
    printf("%d: slsort_parts: pack: %f\n", rank, tpack);
    printf("%d: slsort_parts: alltoall: %f\n", rank, talltoall);
    printf("%d: slsort_parts: alltoallv: %f\n", rank, talltoallv);
#ifdef THREE_STEPS_VERSION
    printf("%d: slsort_parts: postsort: %f\n", rank, tsort);
    printf("%d: slsort_parts: unpack: %f\n", rank, tunpack);
#else
    printf("%d: slsort_parts: mergeunpack: %f\n", rank, tmergeunpack);
#  ifndef MERGE_AND_UNPACK
    printf("%d: slsort_parts: unpack: %f\n", rank, tunpack);
#  endif
#endif
    printf("%d: slsort_parts: makeindices: %f\n", rank, tmakeindices);
# else
    printf("%" slint_fmt "  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", ntotal, ttotal, tinitindxl, tcopy, tsort, tpartition, tpack, talltoall, talltoallv, tmergeunpack, tunpack, tmakeindices);
# endif
  }
#endif

#else

  fprintf(stderr, "ERROR: slsort_parts is not available. Use 'choose_sort = 3' (sl_sort_keys)!\n");

#endif /* NOT_sl_pepcparts */
}
