/*
 *  interpolate.c
 *
 *  Copyright 2006 Matthias Bolten. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "interpolate.h"

void interpolate_prepare( double ***fine, double ***coarse, mg_data* data, int level )
{
  int i, j, k;

  for (i = 1-data[level].x_off;i < data[level].m_l-1; i+=2) {
    for (j = 1-data[level].y_off;j<data[level].n_l-1;j+=2) {
      for (k=1-data[level].z_off;k<data[level].o_l-1;k+=2) {
        fine[i+1][j+1][k+1] = 0.125 * coarse[i/2+1][j/2+1][k/2+1];
      }
    }
  }

  return;
}

void interpolate_finish(double ***fine, mg_data* data, int level )
{
  int i, j, k;
  int xoff = data[level].x_off, yoff = data[level].y_off, zoff = data[level].z_off;
  int m = data[level].m_l, n = data[level].n_l, o = data[level].o_l;

  for (i=1+xoff;i<m-1;i+=2) {
    for (j=1+yoff;j<n-1;j+=2) {
      for (k=1+zoff;k<o-1;k+=2) {
        fine[i][j][k] = 0.125 * (fine[i-1][j-1][k-1] + fine[i-1][j-1][k+1] + fine[i-1][j+1][k-1] + fine[i-1][j+1][k+1] + fine[i+1][j-1][k-1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
      }
    }
  }

 for (i=1+xoff;i<m-1;i+=2) {
    for (j=1+yoff;j<n-1;j+=2) {
      for (k=1-zoff;k<o-1;k+=2) {
        fine[i][j][k+1] = 0.25 * (fine[i-1][j-1][k+1] + fine[i-1][j+1][k+1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k+1]);
      }
    }
  }

 for (i=1+xoff;i<m-1;i+=2) {
    for (j=1-yoff;j<n-1;j+=2) {
      for (k=1+zoff;k<o-1;k+=2) {
        fine[i][j+1][k] = 0.25 * (fine[i-1][j+1][k-1] + fine[i-1][j+1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
      }
    }
  }

 for (i=1+xoff;i<m-1;i+=2) {
    for (j=1-yoff;j<n-1;j+=2) {
      for (k=1-zoff;k<o-1;k+=2) {
        fine[i][j+1][k+1] = 0.5 * (fine[i-1][j+1][k+1] + fine[i+1][j+1][k+1]);
      }
    }
  }

 for (i=1-xoff;i<m-1;i+=2) {
    for (j=1+yoff;j<n-1;j+=2) {
      for (k=1+zoff;k<o-1;k+=2) {
        fine[i+1][j][k] = 0.25 * (fine[i+1][j-1][k-1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
      }
    }
  }

 for (i=1-xoff;i<m-1;i+=2) {
    for (j=1+yoff;j<n-1;j+=2) {
      for (k=1-zoff;k<o-1;k+=2) {
        fine[i+1][j][k+1] = 0.5 * (fine[i+1][j-1][k+1] + fine[i+1][j+1][k+1]);
      }
    }
  }

 for (i=1-xoff;i<m-1;i+=2) {
    for (j=1-yoff;j<n-1;j+=2) {
      for (k=1+zoff;k<o-1;k+=2) {
        fine[i+1][j+1][k] = 0.5 * (fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
      }
    }
  }

  return;
}

