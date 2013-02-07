/*
 *  stencil.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "stencil.h"
#include "cuboid.h"
#include "galerkin.h"

void find_min_max( double a, double b, double c, double d, double e, double f, double g, double* min, double* max )
{

  /* 1. Fall: -2a - 2c - 4e */
  *min = -2*a - 2*c - 4*e;
  *max = -2*a - 2*c - 4*e;

  /* 2. Fall: -2b - 2c - 4f */
  *min = (*min > (-2*b - 2*c - 4*f)) ? (-2*b - 2*c - 4*f) : *min;
  *max = (*max < (-2*b - 2*c - 4*f)) ? (-2*b - 2*c - 4*f) : *max;

  /* 3. Fall: -2a - 2b - 4d */
  *min = (*min > (-2*a - 2*b - 4*d)) ? (-2*a - 2*b - 4*d) : *min;
  *max = (*max < (-2*a - 2*b - 4*d)) ? (-2*a - 2*b - 4*d) : *max;

  /* 4. Fall: 2a + 2b + 2c - 4d - 4e - 4f + 8g */
  *min = (*min > (2*a + 2*b + 2*c - 4*d - 4*e - 4*f + 8*g)) ? (2*a + 2*b + 2*c - 4*d - 4*e - 4*f + 8*g) : *min;
  *max = (*max < (2*a + 2*b + 2*c - 4*d - 4*e - 4*f + 8*g)) ? (2*a + 2*b + 2*c - 4*d - 4*e - 4*f + 8*g) : *max;

  /* 5. Fall: -2a + 2b + 2c + 4d + 4e - 4f - 8g */
  *min = (*min > (-2*a + 2*b + 2*c + 4*d + 4*e - 4*f - 8*g)) ? (-2*a + 2*b + 2*c + 4*d + 4*e - 4*f - 8*g) : *min;
  *max = (*max < (-2*a + 2*b + 2*c + 4*d + 4*e - 4*f - 8*g)) ? (-2*a + 2*b + 2*c + 4*d + 4*e - 4*f - 8*g) : *max;

  /* 6. Fall: 2a - 2b + 2c + 4d - 4e + 4f - 8g */
  *min = (*min > (2*a - 2*b + 2*c + 4*d - 4*e + 4*f - 8*g)) ? (2*a - 2*b + 2*c + 4*d - 4*e + 4*f - 8*g) : *min;
  *max = (*max < (2*a - 2*b + 2*c + 4*d - 4*e + 4*f - 8*g)) ? (2*a - 2*b + 2*c + 4*d - 4*e + 4*f - 8*g) : *max;

  /* 7. Fall: 2a + 2b - 2c - 4d + 4e + 4f - 8g */
  *min = (*min > (2*a + 2*b - 2*c - 4*d + 4*e + 4*f - 8*g)) ? (2*a + 2*b - 2*c - 4*d + 4*e + 4*f - 8*g) : *min;
  *max = (*max < (2*a + 2*b - 2*c - 4*d + 4*e + 4*f - 8*g)) ? (2*a + 2*b - 2*c - 4*d + 4*e + 4*f - 8*g) : *max;

  /* 8. Fall: -2a + 2b - 2c + 4d - 4e + 4f + 8g */
  *min = (*min > (-2*a + 2*b - 2*c + 4*d - 4*e + 4*f + 8*g)) ? (-2*a + 2*b - 2*c + 4*d - 4*e + 4*f + 8*g) : *min;
  *max = (*max < (-2*a + 2*b - 2*c + 4*d - 4*e + 4*f + 8*g)) ? (-2*a + 2*b - 2*c + 4*d - 4*e + 4*f + 8*g) : *max;

  /* 9. Fall: 2a - 2b - 2c + 4d + 4e - 4f + 8g */
  *min = (*min > (2*a - 2*b - 2*c + 4*d + 4*e - 4*f + 8*g)) ? (2*a - 2*b - 2*c + 4*d + 4*e - 4*f + 8*g) : *min;
  *max = (*max < (2*a - 2*b - 2*c + 4*d + 4*e - 4*f + 8*g)) ? (2*a - 2*b - 2*c + 4*d + 4*e - 4*f + 8*g) : *max;

  /* 10. Fall: -2a - 2b + 2c - 4d + 4e + 4f + 8g */
  *min = (*min > (-2*a - 2*b + 2*c - 4*d + 4*e + 4*f + 8*g)) ? (-2*a - 2*b + 2*c - 4*d + 4*e + 4*f + 8*g) : *min;
  *max = (*max < (-2*a - 2*b + 2*c - 4*d + 4*e + 4*f + 8*g)) ? (-2*a - 2*b + 2*c - 4*d + 4*e + 4*f + 8*g) : *max;
}


double calculate_omega( double a, double b, double c, double d, double e, double f, double g, double min, double max )
{
  double omega;
  double numerator, denominator;

  /* calculate numerator */
  numerator = 2 * ( 2*(a+b+c) + 4*(d+e+f) + 8*g );

  /* calculate denominator */
  denominator = 2 * ( 2*(a+b+c) + 4*(d+e+f) + 8*g ) + min + max;

  omega = numerator/denominator;

  return omega;
}


void set_stencil( mg_data* data, int level, int dim )
{
  /* calculates the stencil for the next level */

  /* loop counter */
  int i, j, k;
  int count;

  double a, b, c, d, e, f, g;
  double min, max;
  double*** Agalerkin;
  Agalerkin = galerkin( data, level, dim );

  /* set coefficients */
  a = Agalerkin[0][1][1];
  b = Agalerkin[1][0][1];
  c = Agalerkin[1][1][0];
  d = Agalerkin[0][0][1];
  e = Agalerkin[0][1][0];
  f = Agalerkin[1][0][0];
  g = Agalerkin[0][0][0];

  find_min_max( a, b, c, d, e, f, g, &min, &max );
  data[level+1].omega = calculate_omega( a, b, c, d, e, f, g, min, max );

  /* find out new stencil size */
  data[level+1].size = 0;
  for( i = 0; i < dim; i++ )
    for( j = 0; j < dim; j++ )
      for( k = 0; k < dim; k++ )
        if( Agalerkin[i][j][k] != 0 )
          data[level+1].size++;

  /* allocate space for stencil data */
   data[level+1].values    = (double*) realloc(data[level+1].values, data[level+1].size*sizeof(double) );
   data[level+1].x_offsets = (int*) realloc(data[level+1].x_offsets, data[level+1].size*sizeof(int) );
   data[level+1].y_offsets = (int*) realloc(data[level+1].y_offsets, data[level+1].size*sizeof(int) );
   data[level+1].z_offsets = (int*) realloc(data[level+1].z_offsets, data[level+1].size*sizeof(int) );

  /* set stencil data */
  data[level+1].values[0]    = Agalerkin[1][1][1];
  data[level+1].x_offsets[0] = 0;
  data[level+1].y_offsets[0] = 0;
  data[level+1].z_offsets[0] = 0;

  count = 1;
  for( i = 0; i < dim; i++ )
    for( j = 0; j < dim; j++ )
      for( k = 0; k < dim; k++ )
        if( Agalerkin[i][j][k] != 0 && !( (i == 1) && (j == 1) && (k == 1) ) ){
          data[level+1].values[count]    = Agalerkin[i][j][k];
          data[level+1].x_offsets[count] = i-1;
          data[level+1].y_offsets[count] = j-1;
          data[level+1].z_offsets[count] = k-1;
          count++;
        }

  cuboid_free( Agalerkin, dim, dim, dim );
}


