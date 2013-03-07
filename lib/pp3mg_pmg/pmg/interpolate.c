#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <stdbool.h>
#include "interpolate.h"

void interpolate_prepare( double ***fine, double ***coarse, mg_data* data, int level )
{
  int i, j, k;

  for (i = data[level].x_ghosts-data[level].x_off; i < data[level].m_l-data[level].x_ghosts; i+=2) {
    for (j = data[level].y_ghosts-data[level].y_off; j < data[level].n_l-data[level].y_ghosts; j+=2) {
      for (k = data[level].z_ghosts-data[level].z_off; k < data[level].o_l-data[level].z_ghosts; k+=2) {
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

  /* WARNING: currently only zeros at (0.0,0.0,0.0) or (M_PI,M_PI,M_PI) are supported! */
  if (data[level].zero_at_pi3==false) {
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j][k] = 0.125 * (fine[i-1][j-1][k-1] + fine[i-1][j-1][k+1] + fine[i-1][j+1][k-1] + fine[i-1][j+1][k+1] + fine[i+1][j-1][k-1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts-zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j][k+1] = 0.25 * (fine[i-1][j-1][k+1] + fine[i-1][j+1][k+1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts-yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j+1][k] = 0.25 * (fine[i-1][j+1][k-1] + fine[i-1][j+1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts-yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts-zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j+1][k+1] = 0.5 * (fine[i-1][j+1][k+1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts-xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i+1][j][k] = 0.25 * (fine[i+1][j-1][k-1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts-xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts-zoff; k<o-data[level].x_ghosts; k+=2) {
        fine[i+1][j][k+1] = 0.5 * (fine[i+1][j-1][k+1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts-xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts-yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i+1][j+1][k] = 0.5 * (fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
  } else {
    /* zero at (\pi,\pi,\pi) */
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j][k] = -0.125 * (fine[i-1][j-1][k-1] + fine[i-1][j-1][k+1] + fine[i-1][j+1][k-1] + fine[i-1][j+1][k+1] + fine[i+1][j-1][k-1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts-zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j][k+1] = 0.25 * (fine[i-1][j-1][k+1] + fine[i-1][j+1][k+1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts-yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j+1][k] = 0.25 * (fine[i-1][j+1][k-1] + fine[i-1][j+1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts+xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts-yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts-zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i][j+1][k+1] = -0.5 * (fine[i-1][j+1][k+1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts-xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i+1][j][k] = 0.25 * (fine[i+1][j-1][k-1] + fine[i+1][j-1][k+1] + fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts-xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts+yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts-zoff; k<o-data[level].x_ghosts; k+=2) {
        fine[i+1][j][k+1] = -0.5 * (fine[i+1][j-1][k+1] + fine[i+1][j+1][k+1]);
	}
      }
    }
    
    for (i = data[level].x_ghosts-xoff; i<m-data[level].x_ghosts; i+=2) {
      for (j = data[level].y_ghosts-yoff; j<n-data[level].x_ghosts; j+=2) {
	for (k = data[level].z_ghosts+zoff; k<o-data[level].x_ghosts; k+=2) {
	  fine[i+1][j+1][k] = -0.5 * (fine[i+1][j+1][k-1] + fine[i+1][j+1][k+1]);
	}
      }
    }
  }
  return;
}

