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

  for (i=data[level].x_ghosts;i<data[level].m_l-data[level].x_ghosts;i++) {
    for (j=data[level].y_ghosts;j<data[level].n_l-data[level].y_ghosts;j++) {
      for (k=data[level].z_ghosts;k<data[level].o_l-data[level].z_ghosts;k++) {
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

  for (i=0;i<data[level].size;i++)
    if (data[level].x_offsets[i] == 0 && data[level].y_offsets[i] == 0 && data[level].z_offsets[i] == 0)
      alpha = 1.0 / data[level].values[i];

  for (i=data[level].x_ghosts;i<data[level].m_l-data[level].x_ghosts;i++) {
    for (j=data[level].y_ghosts;j<data[level].n_l-data[level].y_ghosts;j++) {
      for (k=data[level].z_ghosts;k<data[level].o_l-data[level].z_ghosts;k++) {
	r[i][j][k] = alpha * r[i][j][k];
      }
    }
  }

  return;
}
