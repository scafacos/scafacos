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

/* BOUNDS_XYZ2COORDS_NAME */
/* BOUNDS_XYZ2COORDS_GHOST */
/* BOUNDS_XYZ2COORDS_PERIODIC */
/* BOUNDS_XYZ2COORDS_TRICLINIC */


#ifdef BOUNDS_XYZ2COORDS_TRICLINIC
# define GET_COORD(_d_, _p_, _gb_, _gf_)  (((_p_)[0] - (_gb_)[0]) * (_gf_)[0] + ((_p_)[1] - (_gb_)[1]) * (_gf_)[1] + ((_p_)[2] - (_gb_)[2]) * (_gf_)[2])
#else
# define GET_COORD(_d_, _p_, _gb_, _gf_)  (((_p_)[_d_] - (_gb_)[_d_]) * (_gf_)[_d_])
#endif


static void BOUNDS_XYZ2COORDS_NAME(int *coord, fcs_float *xyz, fcs_float *base, fcs_float *bounds, int *cart_dims
#if defined(BOUNDS_XYZ2COORDS_GHOST) && defined(BOUNDS_XYZ2COORDS_PERIODIC)
  , fcs_int *periodicity
#endif
  , fcs_float *grid_data
  )
{
  fcs_float v;
  fcs_int l, h, m;

#ifdef BOUNDS_XYZ2COORDS_PERIODIC
  fcs_int p;
#endif


  v = GET_COORD(0, xyz, base, &grid_data[GRID_DATA_INVERT_0]);
  l = 0;
  h = cart_dims[0];

#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  p = 0;
  if (periodicity[0])
  {
    while (v < 0.0) { v += 1.0; --p; }
    while (v >= 1.0) { v -= 1.0; ++p; }

  } else
# endif
  {
    if (v < 0.0) h = l++ - 1;
    if (v >= 1.0) l = h--;
  }
#endif

  while (l <= h) { m = (l + h) / 2; if (v < bounds[m]) h = m - 1; else l = m + 1; } 
  coord[0] = l - 1;
#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  coord[0] += p * cart_dims[0];
# endif
#endif


  v = GET_COORD(1, xyz, base, &grid_data[GRID_DATA_INVERT_1]);
  l = cart_dims[0] + 1;
  h = cart_dims[0] + cart_dims[1] + 1;

#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  p = 0;
  if (periodicity[1])
  {
    while (v < 0.0) { v += 1; --p; }
    while (v >= 1.0) { v -= 1; ++p; }

  } else
# endif
  {
    if (v < 0.0) h = l++ - 1;
    if (v >= 1.0) l = h--;
  }
#endif

  while (l <= h) { m = (l + h) / 2; if (v < bounds[m]) h = m - 1; else l = m + 1; }
  coord[1] = l - 1 - (cart_dims[0] + 1);
#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  coord[1] += p * cart_dims[1];
# endif
#endif


  v = GET_COORD(2, xyz, base, &grid_data[GRID_DATA_INVERT_2]);
  l = cart_dims[0] + cart_dims[1] + 2;
  h = cart_dims[0] + cart_dims[1] + cart_dims[2] + 2;

#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  p = 0;
  if (periodicity[2])
  {
    while (v < 0.0) { v += 1; --p; }
    while (v >= 1.0) { v -= 1; ++p; }

  } else
# endif
  {
    if (v < 0.0) h = l++ - 1;
    if (v >= 1.0) l = h--;
  }
#endif

  while (l <= h) { m = (l + h) / 2; if (v < bounds[m]) h = m - 1; else l = m + 1; }
  coord[2] = l - 1 - (cart_dims[0] + cart_dims[1] + 2);
#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  coord[2] += p * cart_dims[2];
# endif
#endif

/*  printf("%f, %f, %f - %f, %f, %f -> %d, %d, %d\n",
    xyz[0], xyz[1], xyz[2], base[0], base[1], base[2], coord[0], coord[1], coord[2]);*/
}


#undef GET_COORD

#undef BOUNDS_XYZ2COORDS_NAME
