#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <stdbool.h>
#include "restrict.h"

void restrict_inj(double ***fine, double ***coarse, mg_data* data, int level )
{
  int i, j, k;

  for (i = data[level].x_ghosts-1+data[level].x_off; i < data[level+1].m_l-data[level].x_ghosts; i++) {
    for (j = data[level].y_ghosts-1+data[level].y_off; j < data[level+1].n_l-data[level].y_ghosts; j++) {
      for (k = data[level].z_ghosts-1+data[level].z_off; k < data[level+1].o_l-data[level].z_ghosts; k++) {
        coarse[i][j][k] = fine[2*i-data[level].x_off]
	  [2*j-data[level].y_off]
	  [2*k-data[level].z_off];
      }
    }
  }

  return;
}

void restrict_fw(double ***fine, double ***coarse, mg_data* data, int level )
{
  int i, j, k;
  int xoff , yoff, zoff;
  xoff = data[level].x_off;
  yoff = data[level].y_off;
  zoff = data[level].z_off;

  /* WARNING: currently only zeros at (0.0,0.0,0.0) or at (\pi,\pi,\pi) are supported! */
  if (data[level].zero_at_pi3==false) {
    for ( i = data[level+1].x_ghosts; i < data[level+1].m_l-data[level+1].x_ghosts; i++ ) {
      for ( j = data[level+1].y_ghosts; j < data[level+1].n_l-data[level+1].y_ghosts; j++ ) {
	for ( k = data[level+1].z_ghosts; k < data[level+1].o_l-data[level+1].z_ghosts; k++ ) {
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
  } else {
    /* zero at (\pi,\pi,\pi) */
    for ( i = data[level+1].x_ghosts; i < data[level+1].m_l-data[level+1].x_ghosts; i++ ) {
      for ( j = data[level+1].y_ghosts; j < data[level+1].n_l-data[level+1].y_ghosts; j++ ) {
	for ( k = data[level+1].z_ghosts; k < data[level+1].o_l-data[level+1].z_ghosts; k++ ) {
	  coarse[i][j][k] = 
	    -0.015625 * fine[2*i-xoff-1][2*j-yoff-1][2*k-zoff-1] +
	    0.031250 * fine[2*i-xoff-1][2*j-yoff-1][2*k-zoff  ] + 
	    -0.015625 * fine[2*i-xoff-1][2*j-yoff-1][2*k-zoff+1] +
	    
	    0.031250 * fine[2*i-xoff-1][2*j-yoff  ][2*k-zoff-1] +
	    -0.062500 * fine[2*i-xoff-1][2*j-yoff  ][2*k-zoff  ] +
	    0.031250 * fine[2*i-xoff-1][2*j-yoff  ][2*k-zoff+1] +
	    
	    -0.015625 * fine[2*i-xoff-1][2*j-yoff+1][2*k-zoff-1] +
	    0.031250 * fine[2*i-xoff-1][2*j-yoff+1][2*k-zoff  ] +
	    -0.015625 * fine[2*i-xoff-1][2*j-yoff+1][2*k-zoff+1] +
	    
	    
	    0.031250 * fine[2*i-xoff  ][2*j-yoff-1][2*k-zoff-1] +
	    -0.062500 * fine[2*i-xoff  ][2*j-yoff-1][2*k-zoff  ] + 
	    0.031250 * fine[2*i-xoff  ][2*j-yoff-1][2*k-zoff+1] +
	    
	    -0.062500 * fine[2*i-xoff  ][2*j-yoff  ][2*k-zoff-1] +
	    0.125000 * fine[2*i-xoff  ][2*j-yoff  ][2*k-zoff  ] +
	    -0.062500 * fine[2*i-xoff  ][2*j-yoff  ][2*k-zoff+1] +
	    
	    0.031250 * fine[2*i-xoff  ][2*j-yoff+1][2*k-zoff-1] +
	    -0.062500 * fine[2*i-xoff  ][2*j-yoff+1][2*k-zoff  ] +
	    0.031250 * fine[2*i-xoff  ][2*j-yoff+1][2*k-zoff+1] +
	    
	    
	    -0.015625 * fine[2*i-xoff+1][2*j-yoff-1][2*k-zoff-1] +
	    0.031250 * fine[2*i-xoff+1][2*j-yoff-1][2*k-zoff  ] + 
	    -0.015625 * fine[2*i-xoff+1][2*j-yoff-1][2*k-zoff+1] +
	    
	    0.031250 * fine[2*i-xoff+1][2*j-yoff  ][2*k-zoff-1] +
	    -0.062500 * fine[2*i-xoff+1][2*j-yoff  ][2*k-zoff  ] +
	    0.031250 * fine[2*i-xoff+1][2*j-yoff  ][2*k-zoff+1] +
	    
	    -0.015625 * fine[2*i-xoff+1][2*j-yoff+1][2*k-zoff-1] +
	    0.031250 * fine[2*i-xoff+1][2*j-yoff+1][2*k-zoff  ] +
	    -0.015625 * fine[2*i-xoff+1][2*j-yoff+1][2*k-zoff+1];
	}
      }
    }
  }

  return;
}
