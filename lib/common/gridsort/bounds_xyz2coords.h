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


static void BOUNDS_XYZ2COORDS_NAME(int *coord, fcs_float *xyz, fcs_float *base, fcs_float *bounds, int *cart_dims
#if defined(BOUNDS_XYZ2COORDS_GHOST) && defined(BOUNDS_XYZ2COORDS_PERIODIC)
  , fcs_int *periodicity, fcs_float *box_size
#endif
  )
{
  fcs_float v;
  fcs_int l, h, m;

#ifdef BOUNDS_XYZ2COORDS_PERIODIC
  fcs_int p;
#endif


  v = xyz[0] - base[0];
  l = 0;
  h = cart_dims[0];

#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  p = 0;
  if (periodicity[0])
  {
    while (v < bounds[l]) { v += box_size[0]; --p; }
    while (v >= bounds[h]) { v -= box_size[0]; ++p; }

  } else
# endif
  {
    if (v < bounds[l]) h = l++ - 1;
    if (v >= bounds[h]) l = h--;
  }
#endif

  while (l <= h) { m = (l + h) / 2; if (v < bounds[m]) h = m - 1; else l = m + 1; } 
  coord[0] = l - 1;
#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  coord[0] += p * cart_dims[0];
# endif
#endif


  v = xyz[1] - base[1];
  l = cart_dims[0] + 1;
  h = cart_dims[0] + cart_dims[1] + 1;

#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  p = 0;
  if (periodicity[1])
  {
    while (v < bounds[l]) { v += box_size[1]; --p; }
    while (v >= bounds[h]) { v -= box_size[1]; ++p; }

  } else
# endif
  {
    if (v < bounds[l]) h = l++ - 1;
    if (v >= bounds[h]) l = h--;
  }
#endif

  while (l <= h) { m = (l + h) / 2; if (v < bounds[m]) h = m - 1; else l = m + 1; }
  coord[1] = l - 1 - (cart_dims[0] + 1);
#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  coord[1] += p * cart_dims[1];
# endif
#endif


  v = xyz[2] - base[2];
  l = cart_dims[0] + cart_dims[1] + 2;
  h = cart_dims[0] + cart_dims[1] + cart_dims[2] + 2;

#ifdef BOUNDS_XYZ2COORDS_GHOST
# ifdef BOUNDS_XYZ2COORDS_PERIODIC
  p = 0;
  if (periodicity[2])
  {
    while (v < bounds[l]) { v += box_size[2]; --p; }
    while (v >= bounds[h]) { v -= box_size[2]; ++p; }

  } else
# endif
  {
    if (v < bounds[l]) h = l++ - 1;
    if (v >= bounds[h]) l = h--;
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


#undef BOUNDS_XYZ2COORDS_NAME
