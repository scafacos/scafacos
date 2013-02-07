/*
 *  lueqf.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <mpi.h>
#include "lueqf.h"

double lueqf_res( double*** v, double*** f, double*** r, mg_data* data, int level )
{
  int i, j, k, count;
  double res = 0.0;
  double stencilsum;

  for (i=1;i<data[level].m_l-1;i++) {
    for (j=1;j<data[level].n_l-1;j++) {
      for (k=1;k<data[level].o_l-1;k++) {
	stencilsum = 0.0;
        for( count = 0; count < data[level].size; count++ )
	  stencilsum += ( data[level].values[count] *
			  v[i+data[level].x_offsets[count]]
			   [j+data[level].y_offsets[count]]
			   [k+data[level].z_offsets[count]] );
	r[i][j][k] = f[i][j][k] - stencilsum;
	res = res + r[i][j][k] * r[i][j][k];
      }
    }
  }


  return(res);
}

void lueqf_invd( double ***r, mg_data* data, int level )
{
  int i, j, k;
  double alpha;

  alpha = 1.0 / data[level].values[0];

  for (i=1;i<data[level].m_l-1;i++) {
    for (j=1;j<data[level].n_l-1;j++) {
      for (k=1;k<data[level].o_l-1;k++) {
	r[i][j][k] = alpha * r[i][j][k];
      }
    }
  }

  return;
}

