/*
  Copyright (C) 2011, 2012, 2013 Michael Hofmann
  
  This file is part of ScaFaCoS.
  
  ScaFaCoS is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ScaFaCoS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser Public License for more details.
  
  You should have received a copy of the GNU Lesser Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* GRIDSORT_FRONT_TPROC_NAME */
/* GRIDSORT_FRONT_TPROC_GHOST */
/* GRIDSORT_FRONT_TPROC_PERIODIC */
/* GRIDSORT_FRONT_TPROC_TRICLINIC */
/* GRIDSORT_FRONT_TPROC_BOUNDS */
/* GRIDSORT_FRONT_TPROC_ZSLICES */
/* GRIDSORT_FRONT_TPROC_ZONLY */
/* GRIDSORT_FRONT_TPROC_RANK_CACHE */


#define CONCAT_(_a_, _b_)  _a_##_b_
#define CONCAT(_a_, _b_)   CONCAT_(_a_, _b_)

#ifdef GRIDSORT_FRONT_TPROC_TRICLINIC
# define GET_COORD(_d_, _p_, _gb_, _gf_)  (int) fcs_floor(((_p_)[0] - (_gb_)[0]) * (_gf_)[0] + ((_p_)[1] - (_gb_)[1]) * (_gf_)[1] + ((_p_)[2] - (_gb_)[2]) * (_gf_)[2])
#else
# define GET_COORD(_d_, _p_, _gb_, _gf_)  (int) fcs_floor(((_p_)[_d_] - (_gb_)[_d_]) * (_gf_)[_d_])
#endif

#ifdef GRIDSORT_FRONT_TPROC_BOUNDS
# ifdef GRIDSORT_FRONT_TPROC_GHOST
#  ifdef GRIDSORT_FRONT_TPROC_PERIODIC
#   ifdef GRIDSORT_FRONT_TPROC_TRICLINIC
#    define BOUNDS_XYZ2COORDS(_c_, _xyz_, _base_, _b_, _cd_, _p_, _gd_)  bounds_xyz2coords_ghost_periodic_tricl((_c_), (_xyz_), (_base_), (_b_), (_cd_), (_p_), (_gd_))
#   else
#    define BOUNDS_XYZ2COORDS(_c_, _xyz_, _base_, _b_, _cd_, _p_, _gd_)  bounds_xyz2coords_ghost_periodic((_c_), (_xyz_), (_base_), (_b_), (_cd_), (_p_), (_gd_))
#   endif
#  else
#   ifdef GRIDSORT_FRONT_TPROC_TRICLINIC
#    define BOUNDS_XYZ2COORDS(_c_, _xyz_, _base_, _b_, _cd_, _p_, _gd_)  bounds_xyz2coords_ghost_tricl((_c_), (_xyz_), (_base_), (_b_), (_cd_), (_gd_))
#   else
#    define BOUNDS_XYZ2COORDS(_c_, _xyz_, _base_, _b_, _cd_, _p_, _gd_)  bounds_xyz2coords_ghost((_c_), (_xyz_), (_base_), (_b_), (_cd_), (_gd_))
#   endif
#  endif
# else
#  ifdef GRIDSORT_FRONT_TPROC_TRICLINIC
#   define BOUNDS_XYZ2COORDS(_c_, _xyz_, _base_, _b_, _cd_, _p_, _gd_)   bounds_xyz2coords_tricl((_c_), (_xyz_), (_base_), (_b_), (_cd_), (_gd_))
#  else
#   define BOUNDS_XYZ2COORDS(_c_, _xyz_, _base_, _b_, _cd_, _p_, _gd_)   bounds_xyz2coords((_c_), (_xyz_), (_base_), (_b_), (_cd_), (_gd_))
#  endif
# endif
#endif

#if defined(GRIDSORT_FRONT_TPROC_GHOST) || defined(GRIDSORT_FRONT_TPROC_ZSLICES)
# define DUPLICATE
#endif

#ifdef GRIDSORT_FRONT_TPROC_RANK_CACHE
# define CART2RANK(_comm_, _coords_, _rank_, _cache_, _pos_)  do { \
  if ((_cache_)[(_pos_)] == MPI_UNDEFINED) MPI_Cart_rank((_comm_), (_coords_), &(_cache_)[(_pos_)]); \
  *(_rank_) = (_cache_)[(_pos_)]; \
} while (0)
#else
# define CART2RANK(_comm_, _coords_, _rank_, _cache_, _pos_)  do { \
  MPI_Cart_rank((_comm_), (_coords_), (_rank_)); \
} while (0)
#endif


#ifndef DUPLICATE
static int CONCAT(GRIDSORT_PREFIX, GRIDSORT_FRONT_TPROC_NAME)(GRIDSORT_ELEM_BUF_T *s, GRIDSORT_ELEM_INDEX_T x, void *data)
#else
static void CONCAT(GRIDSORT_PREFIX, GRIDSORT_FRONT_TPROC_NAME)(GRIDSORT_ELEM_BUF_T *s, GRIDSORT_ELEM_INDEX_T x, void *data, GRIDSORT_INT_T *nprocs, int *procs, GRIDSORT_ELEM_BUF_T *sd)
#endif
{
  void **data_ptrs = data;

  MPI_Comm cart_comm;
  fcs_float *grid_data;

  int base_coords[3];

  int *cart_dims;

#ifdef DUPLICATE
  fcs_int nparts;
  int low_coords[3], high_coords[3], coords[3];
#else
  int rank;
#endif

#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  fcs_int *periodicity;
# ifdef DUPLICATE
#  ifdef GRIDSORT_KEY
  fcs_gridsort_index_t key[3];
#  endif
# endif
#endif

#ifdef GRIDSORT_FRONT_TPROC_BOUNDS
  fcs_float *bounds;
#endif

#if defined(GRIDSORT_FRONT_TPROC_ZSLICES)
  int *cart_coords;
#endif

#ifdef GRIDSORT_FRONT_TPROC_RANK_CACHE
  int *rank_cache;
  fcs_int *rank_cache_sizes;
#endif

#ifdef GRIDSORT_XYZ_PTR
  fcs_float *xyz_ptr;
#endif

#if 0
#define XSTR(_s_)  STR(_s_)
#define STR(_s_)   #_s_
  printf("GRIDSORT_FRONT_TPROC_NAME: " XSTR(GRIDSORT_FRONT_TPROC_NAME) "\n");
#undef XSTR
#undef STR
#endif

  cart_comm = *((MPI_Comm *) data_ptrs[GRID_TPROC_DATA_COMM]);
  grid_data = data_ptrs[GRID_TPROC_DATA_GRID_DATA];
  cart_dims = data_ptrs[GRID_TPROC_DATA_CART_DIMS];
#if defined(GRIDSORT_FRONT_TPROC_ZSLICES)
  cart_coords = data_ptrs[GRID_TPROC_DATA_CART_COORDS];
#endif

#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  periodicity = data_ptrs[GRID_TPROC_DATA_PERIODICITY];
#endif

#ifdef GRIDSORT_FRONT_TPROC_BOUNDS
  bounds = data_ptrs[GRID_TPROC_DATA_BOUNDS];
#endif

#ifdef GRIDSORT_FRONT_TPROC_RANK_CACHE
  rank_cache = data_ptrs[GRID_TPROC_DATA_RANK_CACHE];
  rank_cache_sizes = data_ptrs[GRID_TPROC_DATA_RANK_CACHE_SIZES];
#endif

#ifdef GRIDSORT_XYZ_PTR
  xyz_ptr = data_ptrs[GRID_TPROC_DATA_XYZ];
#endif

#ifdef GRIDSORT_FRONT_TPROC_BOUNDS
  BOUNDS_XYZ2COORDS(base_coords, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_BASE], bounds, cart_dims, periodicity, grid_data);
# ifdef GRIDSORT_FRONT_TPROC_GHOST
  BOUNDS_XYZ2COORDS(low_coords, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_LOW], bounds, cart_dims, periodicity, grid_data);
  BOUNDS_XYZ2COORDS(high_coords, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_HIGH], bounds, cart_dims, periodicity, grid_data);
# endif
#else /* GRIDSORT_FRONT_TPROC_BOUNDS */
# ifdef GRIDSORT_FRONT_TPROC_ZONLY
  base_coords[0] = cart_coords[0];
  base_coords[1] = cart_coords[1];
# else /* GRIDSORT_FRONT_TPROC_ZONLY */
  base_coords[0] = GET_COORD(0, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_BASE], &grid_data[GRID_DATA_INVERT_0]);
  base_coords[1] = GET_COORD(1, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_BASE], &grid_data[GRID_DATA_INVERT_1]);
# endif /* GRIDSORT_FRONT_TPROC_ZONLY */
  base_coords[2] = GET_COORD(2, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_BASE], &grid_data[GRID_DATA_INVERT_2]);
# ifdef DUPLICATE
#  ifdef GRIDSORT_FRONT_TPROC_ZONLY
  low_coords[0] = cart_coords[0];
  low_coords[1] = cart_coords[1];
#  else /* GRIDSORT_FRONT_TPROC_ZONLY */
  low_coords[0] = GET_COORD(0, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_LOW], &grid_data[GRID_DATA_INVERT_0]);
  low_coords[1] = GET_COORD(1, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_LOW], &grid_data[GRID_DATA_INVERT_1]);
#  endif /* GRIDSORT_FRONT_TPROC_ZONLY */
#  ifdef GRIDSORT_FRONT_TPROC_ZSLICES
  low_coords[2] = GET_COORD(2, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_ZSLICES_LOW], &grid_data[GRID_DATA_INVERT_2]);
#  else /* GRIDSORT_FRONT_TPROC_ZSLICES */
  low_coords[2] = GET_COORD(2, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_LOW], &grid_data[GRID_DATA_INVERT_2]);
#  endif /* GRIDSORT_FRONT_TPROC_ZSLICES */
#  ifdef GRIDSORT_FRONT_TPROC_ZONLY
  high_coords[0] = cart_coords[0];
  high_coords[1] = cart_coords[1];
#  else /* GRIDSORT_FRONT_TPROC_ZONLY */
  high_coords[0] = GET_COORD(0, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_HIGH], &grid_data[GRID_DATA_INVERT_0]);
  high_coords[1] = GET_COORD(1, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_HIGH], &grid_data[GRID_DATA_INVERT_1]);
#  endif /* GRIDSORT_FRONT_TPROC_ZONLY */
#  ifdef GRIDSORT_FRONT_TPROC_ZSLICES
  high_coords[2] = GET_COORD(2, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_ZSLICES_HIGH], &grid_data[GRID_DATA_INVERT_2]);
#  else /* GRIDSORT_FRONT_TPROC_ZSLICES */
  high_coords[2] = GET_COORD(2, GRIDSORT_XYZ(s, xyz_ptr, x), &grid_data[GRID_DATA_HIGH], &grid_data[GRID_DATA_INVERT_2]);
#  endif /* GRIDSORT_FRONT_TPROC_ZSLICES */
# endif /* DUPLICATE */
#endif /* GRIDSORT_FRONT_TPROC_BOUNDS */

#ifndef GRIDSORT_FRONT_TPROC_ZONLY
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  if (!periodicity[0])
  {
#endif
    base_coords[0] = z_minmax(0, base_coords[0], cart_dims[0] - 1);
#ifdef GRIDSORT_FRONT_TPROC_GHOST
    low_coords[0] = z_minmax(0, low_coords[0], cart_dims[0] - 1);
    high_coords[0] = z_minmax(0, high_coords[0], cart_dims[0] - 1);
#endif
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  }
#endif
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  if (!periodicity[1])
  {
#endif
    base_coords[1] = z_minmax(0, base_coords[1], cart_dims[1] - 1);
#ifdef GRIDSORT_FRONT_TPROC_GHOST
    low_coords[1] = z_minmax(0, low_coords[1], cart_dims[1] - 1);
    high_coords[1] = z_minmax(0, high_coords[1], cart_dims[1] - 1);
#endif
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  }
#endif
#endif /* GRIDSORT_FRONT_TPROC_ZONLY */
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  if (!periodicity[2])
  {
#endif
    base_coords[2] = z_minmax(0, base_coords[2], cart_dims[2] - 1);
#if defined(GRIDSORT_FRONT_TPROC_GHOST) || defined(GRIDSORT_FRONT_TPROC_ZSLICES)
    low_coords[2] = z_minmax(0, low_coords[2], cart_dims[2] - 1);
    high_coords[2] = z_minmax(0, high_coords[2], cart_dims[2] - 1);
#endif
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
  }
#endif

/*  printf("%d: %f,%f,%f - %d,%d,%d"
#ifdef DUPLICATE
  " - %d,%d,%d - %d,%d,%d"
#endif
  "\n"
  , (int) x, GRIDSORT_XYZ(s, xyz, x)[0], GRIDSORT_XYZ(s, xyz, x)[1], GRIDSORT_XYZ(s, xyz, x)[2],
    base_coords[0], base_coords[1], base_coords[2]
#ifdef DUPLICATE
  , low_coords[0], low_coords[1], low_coords[2], high_coords[0], high_coords[1], high_coords[2]
#endif
    );*/

#ifndef DUPLICATE

  CART2RANK(cart_comm, base_coords, &rank, rank_cache,
    base_coords[0] + rank_cache_sizes[0] * (base_coords[1] + rank_cache_sizes[1] * base_coords[2]));

  return rank;

#else /* DUPLICATE */

  nparts = 0;

#ifndef GRIDSORT_FRONT_TPROC_ZONLY
  for (coords[0] = low_coords[0]; coords[0] <= high_coords[0]; ++coords[0])
  {
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
# ifdef GRIDSORT_KEY
    key[0] = GRIDSORT_GHOST_BASE;

    if (coords[0] < 0)             key[0] |= GRIDSORT_PERIODIC_SET((cart_dims[0] - 1 - coords[0]) / cart_dims[0], 0);
    if (coords[0] >= cart_dims[0]) key[0] |= GRIDSORT_PERIODIC_SET(coords[0] / cart_dims[0], 1);
# endif
#endif

    for (coords[1] = low_coords[1]; coords[1] <= high_coords[1]; ++coords[1])
    {
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
# ifdef GRIDSORT_KEY
      key[1] = key[0];

      if (coords[1] < 0)             key[1] |= GRIDSORT_PERIODIC_SET((cart_dims[1] - 1 - coords[1]) / cart_dims[1], 2);
      if (coords[1] >= cart_dims[1]) key[1] |= GRIDSORT_PERIODIC_SET(coords[1] / cart_dims[1], 3);
# endif
#endif
#else /* GRIDSORT_FRONT_TPROC_ZONLY */
  coords[0] = base_coords[0];
  coords[1] = base_coords[1];
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
# ifdef GRIDSORT_KEY
  key[1] = GRIDSORT_GHOST_BASE;
# endif
#endif
#endif /* GRIDSORT_FRONT_TPROC_ZONLY */
      for (coords[2] = low_coords[2]; coords[2] <= high_coords[2]; ++coords[2])
      {
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
# ifdef GRIDSORT_KEY
        key[2] = key[1];

        if (coords[2] < 0)             key[2] |= GRIDSORT_PERIODIC_SET((cart_dims[2] - 1 - coords[2]) / cart_dims[2], 4);
        if (coords[2] >= cart_dims[2]) key[2] |= GRIDSORT_PERIODIC_SET(coords[2] / cart_dims[2], 5);
# endif
#endif

        CART2RANK(cart_comm, coords, &procs[nparts], rank_cache,
          coords[0] + rank_cache_sizes[0] * (coords[1] + rank_cache_sizes[1] * coords[2]));

        if (sd)
        {
#ifdef GRIDSORT_KEY
          if (coords[0] == base_coords[0] && coords[1] == base_coords[1] && coords[2] == base_coords[2]) GRIDSORT_KEY(sd, nparts) = GRIDSORT_KEY(s, x);
#ifdef GRIDSORT_FRONT_TPROC_PERIODIC
          else GRIDSORT_KEY(sd, nparts) = key[2];
#else
          else GRIDSORT_KEY(sd, nparts) = GRIDSORT_GHOST_BASE;
#endif
#endif

#ifdef GRIDSORT_DATA0
          GRIDSORT_DATA0(sd, nparts)[0] = GRIDSORT_DATA0(s, x)[0];
          GRIDSORT_DATA0(sd, nparts)[1] = GRIDSORT_DATA0(s, x)[1];
          GRIDSORT_DATA0(sd, nparts)[2] = GRIDSORT_DATA0(s, x)[2];
#endif
#ifdef GRIDSORT_DATA1
          GRIDSORT_DATA1(sd, nparts) = GRIDSORT_DATA1(s, x);
#endif
        }

        ++nparts;
      }
#ifndef GRIDSORT_FRONT_TPROC_ZONLY
    }
  }
#endif /* GRIDSORT_FRONT_TPROC_ZONLY */

  *nprocs = nparts;

#endif /* DUPLICATE */
}


#ifdef DEFINE_GRIDSORT_FRONT_TPROC_EXDEF
# ifndef DUPLICATE
DEFINE_GRIDSORT_FRONT_TPROC_EXDEF(CONCAT(CONCAT(GRIDSORT_PREFIX, GRIDSORT_FRONT_TPROC_NAME), _exdef), CONCAT(GRIDSORT_PREFIX, GRIDSORT_FRONT_TPROC_NAME), static);
# else
DEFINE_GRIDSORT_FRONT_TPROCS_MOD_EXDEF(CONCAT(CONCAT(GRIDSORT_PREFIX, GRIDSORT_FRONT_TPROC_NAME), _exdef), CONCAT(GRIDSORT_PREFIX, GRIDSORT_FRONT_TPROC_NAME), static);
# endif
#endif

#undef GET_COORD
#undef BOUNDS_XYZ2COORDS
#undef DUPLICATE
#undef CART2RANK

#undef GS_CONCAT_
#undef GS_CONCAT

#undef GRIDSORT_FRONT_TPROC_NAME
