/*
 *  restrict.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "restrict.h"

void restrict_inj(double ***fine, double ***coarse, mg_data* data, int level )
{
  int i, j, k;

  for (i = data[level].x_off; i < data[level+1].m_l-1; i++) {
    for (j = data[level].y_off; j < data[level+1].n_l-1; j++) {
      for (k = data[level].z_off; k < data[level+1].o_l-1; k++) {
        coarse[i][j][k] = fine[2*i-data[level].x_off][2*j-data[level].y_off][2*k-data[level].z_off];
      }
    }
  }

  return;
}

void restrict_fw(double ***fine, double ***coarse, mg_data* data, int level )
{
  int i, j, k;
  int xoff = data[level].x_off, yoff = data[level].y_off, zoff = data[level].z_off;

  for ( i = 1; i < data[level+1].m_l-1; i++ ) {
    for ( j = 1; j < data[level+1].n_l-1; j++ ) {
      for ( k = 1; k < data[level+1].o_l-1; k++ ) {
        coarse[i][j][k] = 
	  0.015625 * fine[2*i-xoff-1][2*j-yoff-1][2*k-zoff-1] +
	  0.031250 * fine[2*i-xoff-1][2*j-yoff-1][2*k-zoff  ] + 
	  0.015625 * fine[2*i-xoff-1][2*j-yoff-1][2*k-zoff+1] +

	  0.031250 * fine[2*i-xoff-1][2*j-yoff  ][2*k-zoff-1] +
	  0.062500 * fine[2*i-xoff-1][2*j-yoff  ][2*k-zoff  ] +
	  0.031250 * fine[2*i-xoff-1][2*j-yoff  ][2*k-zoff+1] +

	  0.015625 * fine[2*i-xoff-1][2*j-yoff+1][2*k-zoff-1] +
	  0.031250 * fine[2*i-xoff-1][2*j-yoff+1][2*k-zoff  ] +
	  0.015625 * fine[2*i-xoff-1][2*j-yoff+1][2*k-zoff+1] +


	  0.031250 * fine[2*i-xoff  ][2*j-yoff-1][2*k-zoff-1] +
	  0.062500 * fine[2*i-xoff  ][2*j-yoff-1][2*k-zoff  ] + 
	  0.031250 * fine[2*i-xoff  ][2*j-yoff-1][2*k-zoff+1] +

	  0.062500 * fine[2*i-xoff  ][2*j-yoff  ][2*k-zoff-1] +
	  0.125000 * fine[2*i-xoff  ][2*j-yoff  ][2*k-zoff  ] +
	  0.062500 * fine[2*i-xoff  ][2*j-yoff  ][2*k-zoff+1] +

	  0.031250 * fine[2*i-xoff  ][2*j-yoff+1][2*k-zoff-1] +
	  0.062500 * fine[2*i-xoff  ][2*j-yoff+1][2*k-zoff  ] +
	  0.031250 * fine[2*i-xoff  ][2*j-yoff+1][2*k-zoff+1] +


	  0.015625 * fine[2*i-xoff+1][2*j-yoff-1][2*k-zoff-1] +
	  0.031250 * fine[2*i-xoff+1][2*j-yoff-1][2*k-zoff  ] + 
	  0.015625 * fine[2*i-xoff+1][2*j-yoff-1][2*k-zoff+1] +

	  0.031250 * fine[2*i-xoff+1][2*j-yoff  ][2*k-zoff-1] +
	  0.062500 * fine[2*i-xoff+1][2*j-yoff  ][2*k-zoff  ] +
	  0.031250 * fine[2*i-xoff+1][2*j-yoff  ][2*k-zoff+1] +

	  0.015625 * fine[2*i-xoff+1][2*j-yoff+1][2*k-zoff-1] +
	  0.031250 * fine[2*i-xoff+1][2*j-yoff+1][2*k-zoff  ] +
	  0.015625 * fine[2*i-xoff+1][2*j-yoff+1][2*k-zoff+1];
      }
    }
  }
}

void old_restrict_fw(double ***fine, double ***coarse, int m, int n, int o, int xoff, int yoff, int zoff)
{
  int i, j, k;

  for (i=2-xoff;i<m-1;i+=2) {
    for (j=2-yoff;j<n-1;j+=2) {
      for (k=2-zoff;k<o-1;k+=2) {
	coarse[i/2][j/2][k/2] =
	  0.015625 * fine[i-1][j-1][k-1] + 0.03125 * fine[i-1][j  ][k-1] + 0.015625 * fine[i-1][j+1][k-1] +
	  0.03125  * fine[i  ][j-1][k-1] + 0.0625  * fine[i  ][j  ][k-1] + 0.03125  * fine[i  ][j+1][k-1] +
	  0.015625 * fine[i+1][j-1][k-1] + 0.03125 * fine[i+1][j  ][k-1] + 0.015625 * fine[i+1][j+1][k-1] +

	  0.03125  * fine[i-1][j-1][k  ] + 0.0625  * fine[i-1][j  ][k  ] + 0.03125  * fine[i-1][j+1][k  ] +
	  0.0625   * fine[i  ][j-1][k  ] + 0.125   * fine[i  ][j  ][k  ] + 0.0625   * fine[i  ][j+1][k  ] +
	  0.03125  * fine[i+1][j-1][k  ] + 0.0625  * fine[i+1][j  ][k  ] + 0.03125  * fine[i+1][j+1][k  ] +

	  0.015625 * fine[i-1][j-1][k+1] + 0.03125 * fine[i-1][j  ][k+1] + 0.015625 * fine[i-1][j+1][k+1] +
	  0.03125  * fine[i  ][j-1][k+1] + 0.0625  * fine[i  ][j  ][k+1] + 0.03125  * fine[i  ][j+1][k+1] +
	  0.015625 * fine[i+1][j-1][k+1] + 0.03125 * fine[i+1][j  ][k+1] + 0.015625 * fine[i+1][j+1][k+1];
      }
    }
  }

  return;
}
