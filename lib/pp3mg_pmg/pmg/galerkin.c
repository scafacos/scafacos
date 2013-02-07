/*
 *  galerkin.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "galerkin.h"
#include "cuboid.h"
#include "rectangle.h"

double** calculate_R( double* r, int dim )
{
  /* loop counter */
  int i, j;

  double** R1;
  double** R2;
  double** R3;
  double** R;

  R1 = rectangle_alloc( dim, dim );
  R2 = rectangle_alloc( dim, dim );
  R3 = rectangle_alloc( dim, dim );
  R  = rectangle_alloc( dim, dim );

  /* set R1 */
  R1[0][0] = r[1];
  R1[0][1] = 0;
  R1[0][2] = 0;
  R1[1][0] = r[0];
  R1[1][1] = r[2];
  R1[1][2] = 0;
  R1[2][0] = 0;
  R1[2][1] = r[1];
  R1[2][2] = 0;

  /* set R2 */
  R2[0][0] = r[0];
  R2[0][1] = r[2];
  R2[0][2] = 0;
  R2[1][0] = 0;
  R2[1][1] = r[1];
  R2[1][2] = 0;
  R2[2][0] = 0;
  R2[2][1] = r[0];
  R2[2][2] = r[2];

  /* set R3 */
  R3[0][0] = 0;
  R3[0][1] = r[1];
  R3[0][2] = 0;
  R3[1][0] = 0;
  R3[1][1] = r[0];
  R3[1][2] = r[2];
  R3[2][0] = 0;
  R3[2][1] = 0;
  R3[2][2] = r[1];

  /* calculate R */
  for( i = 0; i < dim; i++ )
    for( j = 0; j < dim; j++ )
      R[i][j] = r[0]*R1[i][j] + r[1]*R2[i][j] + r[2]*R3[i][j];

  rectangle_free( R1, dim, dim );
  rectangle_free( R2, dim, dim );
  rectangle_free( R3, dim, dim );

  return R;

}

void part( double*** M, double** R, int k, int l, int dim, int flag, double* partM )
{
  int i, j;

  double sum;
  double* tmp;

  tmp   = (double*)malloc( dim*sizeof(double) );

  if( flag == 0 ){
    for( i = 0; i < dim; i++ )
      tmp[i] = M[k][i][l];
  }
  else if( flag == 1 ){
    for( i = 0; i < dim; i++ )
      tmp[i] = M[i][k][l];
  }
  else if( flag == 2 ){
    for( i = 0; i < dim; i++ )
      tmp[i] = M[k][l][i];
  }

  for( i = 0; i < dim; i++ ){
    sum = 0.0;
    for( j = 0; j < dim; j++ )
      sum += tmp[j]*R[j][i];
    partM[i] = sum;
  }

  free( tmp );
}

double*** galerkin( mg_data* data, int level, int dim )
{
  /* loop counter */
  int i, j, k;

  double r[] = { 0.25, 0.5, 0.25 };
  double** R;

  double*** A;
  double*** B;
  double*** C;
  double*** D;

  double* tmpvec;
  tmpvec = (double*)malloc( dim*sizeof(double) );

  R = calculate_R( r, dim );

  A = cuboid_alloc( dim, dim, dim );
  B = cuboid_alloc( dim, dim, dim );
  C = cuboid_alloc( dim, dim, dim );
  D = cuboid_alloc( dim, dim, dim );

  /* init A, B, C, D */
  for( i = 0; i < dim; i++ )
    for( j = 0; j < dim; j++ )
      for( k = 0; k < dim; k++ ){
        A[i][j][k] = 0.0;
        B[i][j][k] = 0.0;
        C[i][j][k] = 0.0;
        D[i][j][k] = 0.0;
    }

  /* init A with stencil data */
  for( i = 0; i < data[level].size; i++ )
    A[1+data[level].x_offsets[i]][1+data[level].y_offsets[i]][1+data[level].z_offsets[i]] = data[level].values[i];

  /* calculate B */
  for( i = 0; i < dim; i++ )
    for( j = 0; j < dim; j++ ){
      part( A, R, i, j, dim, 0, tmpvec );
      for( k = 0; k < dim; k++ )
        B[i][k][j] = tmpvec[k];
    }

  /* calculate C */
  for( i = 0; i < dim; i++ )
    for( j = 0; j < dim; j++ ){
      part( B, R, i, j, dim, 1, tmpvec );
      for( k = 0; k < dim; k++ )
        C[k][i][j] = tmpvec[k];
    }

  if( dim == 2 ) {
    double*** tmp = C;
    C = D;
    D = tmp;
    goto out;
  }

  /* calculate D */
  for( i = 0; i < dim; i++ )
    for( j = 0; j < dim; j++ ){
      part( C, R, i, j, dim, 2, tmpvec );
      for( k = 0; k < dim; k++ )
        D[j][i][k] = tmpvec[k];
    }

// #ifdef DEBUG
//   printf( "Agalerkin (D):\n" );
// 
//   for( i = 0; i < dim; i++ ){
//     for( j = 0; j < dim; j++ ){
//       for( k = 0; k < dim; k++ )
//         printf( "%3.4f ", D[k][j][i] );
//     printf( "\n" );
//     }
//     printf( "\n\n" );
//   }
// #endif

out:
  /* free pointers */
  free( tmpvec );
  rectangle_free( R, dim, dim );
  cuboid_free( A, dim, dim, dim );
  cuboid_free( B, dim, dim, dim );
  cuboid_free( C, dim, dim, dim );

  return D;

}

